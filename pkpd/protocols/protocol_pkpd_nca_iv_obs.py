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

import pyworkflow.protocol.params as params
from .protocol_pkpd_sa_base import ProtPKPDSABase
from pkpd.models.sa_models import NCAObsIVModel
from pkpd.pkpd_units import PKPDUnit
from pyworkflow.protocol.constants import LEVEL_ADVANCED

# TESTED in test_workflow_gabrielsson_pk01.py
# TESTED in test_workflow_gabrielsson_pk06.py
# TESTED in test_workflow_gabrielsson_pk07.py

class ProtPKPDNCAIVObs(ProtPKPDSABase):
    """ Non-compartmental analysis based on observations.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'nca iv observations'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        ProtPKPDSABase._defineParams1(self,form,False)
        form.addParam('protElimination', params.PointerParam, label="Elimination rate",
                      pointerClass='ProtPKPDEliminationRate',
                      help='Select an execution of a protocol estimating the elimination rate')
        form.addParam('areaCalc', params.EnumParam, choices=["Trapezoidal","Log-Trapezoidal"],
                      label="Method for AUC, AUMC calculation", default=1, expertLevel=LEVEL_ADVANCED,
                      help='See explanation at http://learnpkpd.com/2011/04/02/calculating-auc-linear-and-log-linear\n')
        form.addParam("absorptionF", params.FloatParam, label="Absorption fraction (bioavailability)", default=1,
                      expertLevel=LEVEL_ADVANCED, help="Between 0 (=no absorption) and 1 (=full absorption)")

    def getListOfFormDependencies(self):
        return [self.protElimination.get().getObjId()]

    #--------------------------- STEPS functions --------------------------------------------
    def setupFromFormParameters(self):
        self.fitting = self.readFitting(self.protElimination.get().outputFitting.fnFitting.get())

    def getXYvars(self):
        self.varNameX = self.protElimination.get().predictor.get()
        self.varNameY = self.protElimination.get().predicted.get()

    def createAnalysis(self):
        self.analysis = NCAObsIVModel()
        self.analysis.setExperiment(self.experiment)
        self.analysis.setXVar(self.varNameX)
        self.analysis.setYVar(self.varNameY)
        self.analysis.F = self.absorptionF.get()
        if self.areaCalc == 0:
            self.analysis.areaCalc = "Trapezoidal"
        else:
            self.analysis.areaCalc = "Log-Trapezoidal"

    def prepareForSampleAnalysis(self, sampleName):
        sampleFit = self.fitting.getSampleFit(sampleName)
        sample = self.experiment.samples[sampleName]
        sample.interpretDose()
        self.analysis.D = sample.getDoseAt(0.0)
        if sampleFit == None:
            print("  Cannot process %s because its elimination rate cannot be found\n\n"%sampleName)
            return False
        self.analysis.lambdaz = sampleFit.parameters[1]
        self.analysis.lambdazUnits = PKPDUnit()
        self.analysis.lambdazUnits.unit = self.fitting.modelParameterUnits[1]
        print("Elimination rate = %f [%s]"%(self.analysis.lambdaz,self.analysis.lambdazUnits._toString()))
        return True

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        msg.append("Non-compartmental analysis for the observations of the variable %s"%self.protElimination.get().predicted.get())
        return msg

    def _warnings(self):
        experiment = self.readExperiment(self.getInputExperiment().fnPKPD,show=False)
        incorrectList = experiment.getNonBolusDoses()
        if len(incorrectList)==0:
            return []
        else:
            return ["This protocol is meant only for intravenous bolus regimens. Check the doses for %s"%(','.join(incorrectList))]

