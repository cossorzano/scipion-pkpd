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
                      expertLevel=LEVEL_ADVANCED,
                      help="Multiply the dose in the deposition file by this factor")
        form.addParam('substanceMultiplier', params.StringParam, default="1 1 1 1 1 1 1 1", label='Substance multiplier',
                      expertLevel=LEVEL_ADVANCED,
                      help='Multiply a substance related parameter by this factor. This is a vector with the order: \n'
                           '1. Solubility (br) \n'                 
                           '2. Solubility (alv) \n'
                           '3. Dissolution rate (br) \n'
                           '4. Dissolution rate (alv) \n'
                           '5. Permeability (br) \n'
                           '6. Permeability (alv) \n'
                           '7. Partition coefficient (br) \n'
                           '8. Partition coefficient (alv)')
        form.addParam('physiologyMultiplier', params.StringParam, default="1 1 1 1 1 1 1 1 1", label='Physiology multiplier',
                      expertLevel=LEVEL_ADVANCED,
                      help='Multiply a physiology related parameter by this factor. This is a vector with the order: \n'
                           '1. Mucociliary clearance\n'           
                           '2. Fluid volume (br)\n'
                           '3. Fluid volume (alv)\n'
                           '4. Tissue volume (br)\n'
                           '5. Tissue volume (alv)\n'
                           '6. Surface area (br)\n'
                           '7. Surface area (alv)\n'
                           '8. Perfusion (br)\n'
                           '9. Perfusion (alv)\n')
        form.addParam('ptrPK', params.PointerParam, pointerClass='PKPDExperiment', label="PK parameters")
        form.addParam('pkMultiplier', params.StringParam, default="1 1 1 1 1 1", label='PK multiplier',
                      expertLevel=LEVEL_ADVANCED,
                      help='Multiply a substance related parameter by this factor. This is a vector with the order: \n'
                           '1. Systemic clearance (Cl)\n'
                           '2. Systemic apparent volume (V)\n'
                           '3. Flow between compartments (Q)\n'
                           '4. Peripheral apparent volume (Vp)\n'
                           '5. Absorption (k01)\n'
                           '6. Oral bioavailability (F)\n'
                      )
        form.addParam('ciliarySpeedType', params.EnumParam, choices=['Exponential', 'Fit', 'Interpolated'], default=2,
                      label='Ciliary speed', expertLevel=LEVEL_ADVANCED)

        form.addParam('simulationTime', params.FloatParam, label="Simulation time (min)", default=10*24*60)
        form.addParam('deltaT', params.FloatParam, label='Time step (min)', default=1, expertLevel=LEVEL_ADVANCED)
        form.addParam('diameters', params.StringParam, label='Diameters (um)',
                      default="0.1,1.1,0.1; 1.2,9.2,0.2",
                      help='Diameters to analyze. Syntax: start1, stop1, step1; start2, stop2, step2; ... They will be '
                           'used as arguments of numpy.arange',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('volMultiplier', params.FloatParam, label='Particle volume multiplier', default=1,
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
        substanceParams.multiplier = [float(x) for x in self.substanceMultiplier.get().split()]
        substanceParams.read(self.ptrDeposition.get().fnSubstance.get())

        lungParams = PKPhysiologyLungParameters()
        lungParams.multiplier = [float(x) for x in self.physiologyMultiplier.get().split()]
        lungParams.read(self.ptrDeposition.get().fnLung.get())

        pkParams = PKPDExperiment()
        pkParams.load(self.ptrPK.get().fnPKPD)

        pkLungParams = PKLung()
        pkLungParams.prepare(substanceParams, lungParams, pkParams,
                             [float(x) for x in self.pkMultiplier.get().split()],
                             self.ciliarySpeedType.get())

        # diameters = np.concatenate((np.arange(0.1,1.1,0.1),np.arange(1.2,9.2,0.2))) # [um]
        evalStr = "np.concatenate(("+",".join(["np.arange("+x.strip()+")" for x in self.diameters.get().split(";")])+"))"
        diameters = eval(evalStr, {'np': np})
        Sbnd = self.volMultiplier.get() * diam2vol(diameters)

        tt=np.arange(0,self.simulationTime.get()+self.deltaT.get(),self.deltaT.get())
        sol=saturable_2D_upwind_IE(lungParams, pkLungParams, self.deposition, tt, Sbnd)

        # Postprocessing
        depositionData = self.deposition.getData()
        alvDose = np.sum(depositionData['alveolar'])
        bronchDose = np.sum(depositionData['bronchial'])
        lungDose = alvDose + bronchDose
        Aalvsolid = sol['A']['alv']['solid']
        Aalvfluid = sol['A']['alv']['fluid']
        Aalvtissue = sol['A']['alv']['tissue']
        AsysGut = sol['A']['sys']['gut']
        AsysPer = sol['A']['sys']['per']
        AsysCtr = sol['A']['sys']['ctr']
        Atisbr = sol['A']['br']['tissue']
        Abrtissue = Atisbr
        Vtisbr = lungParams.getBronchial()['fVol'] * lungParams.getSystemic()['OWlung']
        Cavgbr = Atisbr / Vtisbr;
        Abrcleared = sol['A']['br']['clear']
        Abrcleared = np.reshape(Abrcleared,Abrcleared.size)
        Abrsolid = sol['A']['br']['solid']
        Abrfluid = sol['A']['br']['fluid']
        Acleared = sol['A']['sys']['clear'] + AsysGut - AsysGut[0]
        lungRetention = 100*(lungDose - Acleared)/lungDose;  # in percent of lung dose

        Csysnmol = sol['C']['sys']['ctr']
        Csys = Csysnmol * substanceParams.getData()['MW']

        CsysPer = AsysPer * substanceParams.getData()['MW'] / pkLungParams.pkData['Vp']

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

        alvTissueVar = PKPDVariable()
        alvTissueVar.varName = "alvTissue"
        alvTissueVar.varType = PKPDVariable.TYPE_NUMERIC
        alvTissueVar.role = PKPDVariable.ROLE_MEASUREMENT
        alvTissueVar.units = createUnit("nmol")
        alvTissueVar.comment = "Amount in alveoli"

        alvSolidVar = PKPDVariable()
        alvSolidVar.varName = "alvSolid"
        alvSolidVar.varType = PKPDVariable.TYPE_NUMERIC
        alvSolidVar.role = PKPDVariable.ROLE_MEASUREMENT
        alvSolidVar.units = createUnit("nmol")
        alvSolidVar.comment = "Amount undissolved in alveoli"

        alvFluidVar = PKPDVariable()
        alvFluidVar.varName = "alvFluid"
        alvFluidVar.varType = PKPDVariable.TYPE_NUMERIC
        alvFluidVar.role = PKPDVariable.ROLE_MEASUREMENT
        alvFluidVar.units = createUnit("nmol")
        alvFluidVar.comment = "Amount in alveolar lining fluid"

        brTissueVar = PKPDVariable()
        brTissueVar.varName = "brTissue"
        brTissueVar.varType = PKPDVariable.TYPE_NUMERIC
        brTissueVar.role = PKPDVariable.ROLE_MEASUREMENT
        brTissueVar.units = createUnit("nmol")
        brTissueVar.comment = "Amount in bronchii"

        brSolidVar = PKPDVariable()
        brSolidVar.varName = "brSolid"
        brSolidVar.varType = PKPDVariable.TYPE_NUMERIC
        brSolidVar.role = PKPDVariable.ROLE_MEASUREMENT
        brSolidVar.units = createUnit("nmol")
        brSolidVar.comment = "Amount undissolved in bronchii"

        brFluidVar = PKPDVariable()
        brFluidVar.varName = "brFluid"
        brFluidVar.varType = PKPDVariable.TYPE_NUMERIC
        brFluidVar.role = PKPDVariable.ROLE_MEASUREMENT
        brFluidVar.units = createUnit("nmol")
        brFluidVar.comment = "Amount in bronchial lining fluid"

        brClvar = PKPDVariable()
        brClvar.varName = "brClear"
        brClvar.varType = PKPDVariable.TYPE_NUMERIC
        brClvar.role = PKPDVariable.ROLE_MEASUREMENT
        brClvar.units = createUnit("nmol")
        brClvar.comment = "Cumulative amount cleared by mucociliary elevator"

        brCTisvar = PKPDVariable()
        brCTisvar.varName = "CbrTis"
        brCTisvar.varType = PKPDVariable.TYPE_NUMERIC
        brCTisvar.role = PKPDVariable.ROLE_MEASUREMENT
        brCTisvar.units = createUnit("nmol/mL")
        brCTisvar.comment = "Concentration in bronchial tissue"

        sysGutVar = PKPDVariable()
        sysGutVar.varName = "sysAbsorption"
        sysGutVar.varType = PKPDVariable.TYPE_NUMERIC
        sysGutVar.role = PKPDVariable.ROLE_MEASUREMENT
        sysGutVar.units = createUnit("nmol")
        sysGutVar.comment = "Amount in absorption compartment"

        sysCtrVar = PKPDVariable()
        sysCtrVar.varName = "sysCentral"
        sysCtrVar.varType = PKPDVariable.TYPE_NUMERIC
        sysCtrVar.role = PKPDVariable.ROLE_MEASUREMENT
        sysCtrVar.units = createUnit("nmol")
        sysCtrVar.comment = "Amount in central compartment"

        sysPerVar = PKPDVariable()
        sysPerVar.varName = "sysPeripheral"
        sysPerVar.varType = PKPDVariable.TYPE_NUMERIC
        sysPerVar.role = PKPDVariable.ROLE_MEASUREMENT
        sysPerVar.units = createUnit("nmol")
        sysPerVar.comment = "Amount in peripheral compartment"

        CsysPerVar = PKPDVariable()
        CsysPerVar.varName = "Cp"
        CsysPerVar.varType = PKPDVariable.TYPE_NUMERIC
        CsysPerVar.role = PKPDVariable.ROLE_MEASUREMENT
        CsysPerVar.units = createUnit("g/mL")
        CsysPerVar.comment = "Concentration in peripheral compartment"

        doseNmolVar = PKPDVariable()
        doseNmolVar.varName = "dose_nmol"
        doseNmolVar.varType = PKPDVariable.TYPE_NUMERIC
        doseNmolVar.role = PKPDVariable.ROLE_LABEL
        doseNmolVar.units = createUnit("nmol")
        doseNmolVar.comment = "Input dose in nmol"

        doseThroatVar = PKPDVariable()
        doseThroatVar.varName = "throat_dose_nmol"
        doseThroatVar.varType = PKPDVariable.TYPE_NUMERIC
        doseThroatVar.role = PKPDVariable.ROLE_LABEL
        doseThroatVar.units = createUnit("nmol")
        doseThroatVar.comment = "Throat dose in nmol"

        doseLungVar = PKPDVariable()
        doseLungVar.varName = "lung_dose_nmol"
        doseLungVar.varType = PKPDVariable.TYPE_NUMERIC
        doseLungVar.role = PKPDVariable.ROLE_LABEL
        doseLungVar.units = createUnit("nmol")
        doseLungVar.comment = "Lung dose in nmol"

        doseBronchialVar = PKPDVariable()
        doseBronchialVar.varName = "bronchial_dose_nmol"
        doseBronchialVar.varType = PKPDVariable.TYPE_NUMERIC
        doseBronchialVar.role = PKPDVariable.ROLE_LABEL
        doseBronchialVar.units = createUnit("nmol")
        doseBronchialVar.comment = "Bronchial dose in nmol"

        doseAlveolarVar = PKPDVariable()
        doseAlveolarVar.varName = "alveolar_dose_nmol"
        doseAlveolarVar.varType = PKPDVariable.TYPE_NUMERIC
        doseAlveolarVar.role = PKPDVariable.ROLE_LABEL
        doseAlveolarVar.units = createUnit("nmol")
        doseAlveolarVar.comment = "Alveolar dose in nmol"

        mccClearedLungDoseFractionVar = PKPDVariable()
        mccClearedLungDoseFractionVar.varName = "mcc_cleared_lung_dose_fraction"
        mccClearedLungDoseFractionVar.varType = PKPDVariable.TYPE_NUMERIC
        mccClearedLungDoseFractionVar.role = PKPDVariable.ROLE_LABEL
        mccClearedLungDoseFractionVar.units = createUnit("None")
        mccClearedLungDoseFractionVar.comment = "MCC cleared lung dose fraction"

        self.experimentLungRetention.variables["t"] = tvar
        self.experimentLungRetention.variables["Retention"] = Rvar
        self.experimentLungRetention.variables["Cnmol"] = Cnmolvar
        self.experimentLungRetention.variables["C"] = Cvar
        self.experimentLungRetention.variables["alvTissue"] = alvTissueVar
        self.experimentLungRetention.variables["alvFluid"] = alvFluidVar
        self.experimentLungRetention.variables["alvSolid"] = alvSolidVar
        self.experimentLungRetention.variables["brTissue"] = brTissueVar
        self.experimentLungRetention.variables["brFluid"] = brFluidVar
        self.experimentLungRetention.variables["brSolid"] = brSolidVar
        self.experimentLungRetention.variables["brClear"] = brClvar
        self.experimentLungRetention.variables["CbrTis"] = brCTisvar
        self.experimentLungRetention.variables["sysCentral"] = sysCtrVar
        self.experimentLungRetention.variables["sysAbsoprtion"] = sysGutVar
        self.experimentLungRetention.variables["sysPeripheral"] = sysPerVar
        self.experimentLungRetention.variables["Cp"] = CsysPerVar
        self.experimentLungRetention.variables["dose_nmol"] = doseNmolVar
        self.experimentLungRetention.variables["throat_dose_nmol"] = doseThroatVar
        self.experimentLungRetention.variables["lung_dose_nmol"] = doseLungVar
        self.experimentLungRetention.variables["bronchial_dose_nmol"] = doseBronchialVar
        self.experimentLungRetention.variables["alveolar_dose_nmol"] = doseAlveolarVar
        self.experimentLungRetention.variables["mcc_cleared_lung_dose_fraction"] = mccClearedLungDoseFractionVar

        # Samples
        simulationSample = PKPDSample()
        simulationSample.sampleName = "simulation"
        simulationSample.addMeasurementColumn("t", tt)
        simulationSample.addMeasurementColumn("Retention", lungRetention)
        simulationSample.addMeasurementColumn("Cnmol", Csysnmol)
        simulationSample.addMeasurementColumn("C", Csys)
        simulationSample.addMeasurementColumn("alvFluid", Aalvfluid)
        simulationSample.addMeasurementColumn("alvTissue", Aalvtissue)
        simulationSample.addMeasurementColumn("alvSolid", Aalvsolid)
        simulationSample.addMeasurementColumn("brFluid", Abrfluid)
        simulationSample.addMeasurementColumn("brTissue", Abrtissue)
        simulationSample.addMeasurementColumn("brSolid", Abrsolid)
        simulationSample.addMeasurementColumn("brClear", Abrcleared)
        simulationSample.addMeasurementColumn("CbrTis", Cavgbr)
        simulationSample.addMeasurementColumn("sysAbsorption", AsysGut)
        simulationSample.addMeasurementColumn("sysCentral", AsysCtr)
        simulationSample.addMeasurementColumn("sysPeripheral", AsysPer)
        simulationSample.addMeasurementColumn("Cp", CsysPer)
        simulationSample.setDescriptorValue("dose_nmol",depositionData['dose_nmol'])
        simulationSample.setDescriptorValue("throat_dose_nmol",depositionData['throat'])
        lungDose = depositionData['dose_nmol']-depositionData['throat']
        simulationSample.setDescriptorValue("lung_dose_nmol",lungDose)
        simulationSample.setDescriptorValue("bronchial_dose_nmol",bronchDose)
        simulationSample.setDescriptorValue("alveolar_dose_nmol",alvDose)
        simulationSample.setDescriptorValue("mcc_cleared_lung_dose_fraction",Abrcleared[-1]/lungDose)
        self.experimentLungRetention.samples["simulmvarNameation"] = simulationSample

        self.experimentLungRetention.write(self._getPath("experiment.pkpd"))

        # Plots
        import matplotlib.pyplot as plt
        lungData = lungParams.getBronchial()
        Xbnd = np.sort([0] + lungData['end_cm'].tolist() + lungData['pos'].tolist())
        Xctr = Xbnd[:-1] + np.diff(Xbnd) / 2

        tvec = tt / 60;
        T = np.max(tvec);

        Cflu = sol['C']['br']['fluid'];
        Cs = substanceParams.getData()['Cs_br'];

        plt.figure(figsize=(15, 9))
        plt.title('Concentration in bronchial fluid (Cs=%f [uM])'%Cs)
        plt.imshow(Cflu, interpolation='bilinear', aspect='auto', extent=[np.min(Xctr),np.max(Xctr),T,0])
        plt.clim(0,Cs)
        plt.xlim(np.min(Xctr),np.max(Xctr))
        plt.ylim(0,T)
        plt.colorbar()
        plt.ylabel('Time [h]')
        plt.xlabel('Distance from throat [cm]')
        plt.savefig(self._getPath('concentrationBronchialFluid.png'))

        Ctis = sol['C']['br']['fluid'];
        plt.figure(figsize=(15, 9))
        plt.title('Concentration in bronchial tissue')
        plt.imshow(Ctis, interpolation='bilinear', aspect='auto', extent=[np.min(Xctr),np.max(Xctr),T,0])
        plt.xlim(np.min(Xctr),np.max(Xctr))
        plt.ylim(0,T)
        plt.colorbar()
        plt.ylabel('Time [h]')
        plt.xlabel('Distance from throat [cm]')
        plt.savefig(self._getPath('concentrationBronchialTissue.png'))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experimentLungRetention)
        self._defineSourceRelation(self.ptrDeposition.get(), self.experimentLungRetention)
        self._defineSourceRelation(self.ptrPK.get(), self.experimentLungRetention)

    def _citations(self):
        return ['Hartung2020']