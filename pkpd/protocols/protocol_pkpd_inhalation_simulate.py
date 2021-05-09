# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (info@kinestat.com)
# *
# * Kinestat Pharma
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'info@kinestat.com'
# *
# **************************************************************************

import math
import numpy as np

import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from .protocol_pkpd import ProtPKPD

from pkpd.objects import PKDepositionParameters, PKSubstanceLungParameters, PKPhysiologyLungParameters, PKLung,\
                         PKPDExperiment, PKPDVariable, PKPDSample
from pkpd.pkpd_units import createUnit
from pkpd.inhalation import diam2vol, saturable_2D_upwind_IE

# Tested in test_workflow_inhalation1

class ProtPKPDInhSimulate(ProtPKPD):
    """ Simulate inhalation PK\n
        See Hartung2020.
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'simulate inhalation'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('ptrDeposition', params.PointerParam, pointerClass='PKDepositionParameters', label="Deposition")
        form.addParam('doseMultiplier', params.FloatParam, default=1.0, label='Dose multiplier',
                      help="Multiply the dose in the deposition file by this factor")
        form.addParam('ptrPK', params.PointerParam, pointerClass='PKPDExperiment', label="PK parameters")
        form.addParam('simulationTime', params.FloatParam, label="Simulation time (min)", default=10*24*60)
        form.addParam('deltaT', params.FloatParam, label='Time step (min)', default=1, expertLevel=LEVEL_ADVANCED)
        form.addParam('diameters', params.StringParam, label='Diameters (um)',
                      default="0.1,1.1,0.1; 1.2,9.2,0.2",
                      help='Diameters to analyze. Syntax: start1, stop1, step1; start2, stop2, step2; ... They will be '
                           'used as arguments of numpy.arange',
                      expertLevel=LEVEL_ADVANCED)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulation')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runSimulation(self):
        self.deposition = PKDepositionParameters()
        self.deposition.setFiles(self.ptrDeposition.get().fnSubstance.get(),
                                 self.ptrDeposition.get().fnLung.get(),
                                 self.ptrDeposition.get().fnDeposition.get())
        self.deposition.doseMultiplier = self.doseMultiplier.get()
        self.deposition.read()

        substanceParams = PKSubstanceLungParameters()
        substanceParams.read(self.ptrDeposition.get().fnSubstance.get())

        lungParams = PKPhysiologyLungParameters()
        lungParams.read(self.ptrDeposition.get().fnLung.get())

        pkParams = PKPDExperiment()
        pkParams.load(self.ptrPK.get().fnPKPD)

        pkLungParams = PKLung()
        pkLungParams.prepare(substanceParams, lungParams, pkParams)

        # diameters = np.concatenate((np.arange(0.1,1.1,0.1),np.arange(1.2,9.2,0.2))) # [um]
        evalStr = "np.concatenate(("+",".join(["np.arange("+x.strip()+")" for x in self.diameters.get().split(";")])+"))"
        diameters = eval(evalStr, {'np': np})
        Sbnd = diam2vol(diameters)

        tt=np.arange(0,self.simulationTime.get()+self.deltaT.get(),self.deltaT.get())
        sol=saturable_2D_upwind_IE(lungParams, pkLungParams, self.deposition, tt, Sbnd)

        # Postprocessing
        depositionData = self.deposition.getData()
        alvDose = np.sum(depositionData['alveolar'])
        bronchDose = np.sum(depositionData['bronchial'])
        lungDose = alvDose + bronchDose
        AsysGut = sol['A']['sys']['gut']
        Acleared = sol['A']['sys']['clear'] + AsysGut - AsysGut[0]
        lungRetention = 100*(lungDose - Acleared)/lungDose;  # in percent of lung dose

        Csysnmol = sol['C']['sys']['ctr']
        Csys = Csysnmol * substanceParams.getData()['MW']

        # Create output
        self.experimentLungRetention = PKPDExperiment()
        self.experimentLungRetention.general["title"]="Inhalation simulate"

        tvar = PKPDVariable()
        tvar.varName = "t"
        tvar.varType = PKPDVariable.TYPE_NUMERIC
        tvar.role = PKPDVariable.ROLE_TIME
        tvar.units = createUnit("min")

        Rvar = PKPDVariable()
        Rvar.varName = "Retention"
        Rvar.varType = PKPDVariable.TYPE_NUMERIC
        Rvar.role = PKPDVariable.ROLE_MEASUREMENT
        Rvar.units = createUnit("none")
        Rvar.comment = "Lung retention (% lung dose)"

        Cnmolvar = PKPDVariable()
        Cnmolvar.varName = "Cnmol"
        Cnmolvar.varType = PKPDVariable.TYPE_NUMERIC
        Cnmolvar.role = PKPDVariable.ROLE_MEASUREMENT
        Cnmolvar.units = createUnit("nmol/mL")
        Cnmolvar.comment = "Central compartment concentration"

        Cvar = PKPDVariable()
        Cvar.varName = "C"
        Cvar.varType = PKPDVariable.TYPE_NUMERIC
        Cvar.role = PKPDVariable.ROLE_MEASUREMENT
        Cvar.units = createUnit("g/mL")
        Cvar.comment = "Central compartment concentration"

        self.experimentLungRetention.variables["t"] = tvar
        self.experimentLungRetention.variables["Retention"] = Rvar
        self.experimentLungRetention.variables["Cnmol"] = Cnmolvar
        self.experimentLungRetention.variables["C"] = Cvar

        # Samples
        simulationSample = PKPDSample()
        simulationSample.sampleName = "simulation"
        simulationSample.addMeasurementColumn("t", tt)
        simulationSample.addMeasurementColumn("Retention", lungRetention)
        simulationSample.addMeasurementColumn("Cnmol", Csysnmol)
        simulationSample.addMeasurementColumn("C", Csys)
        self.experimentLungRetention.samples["simulmvarNameation"] = simulationSample

        self.experimentLungRetention.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experimentLungRetention)
        self._defineSourceRelation(self.ptrDeposition.get(), self.experimentLungRetention)
        self._defineSourceRelation(self.ptrPK.get(), self.experimentLungRetention)

    def _citations(self):
        return ['Hartung2020']