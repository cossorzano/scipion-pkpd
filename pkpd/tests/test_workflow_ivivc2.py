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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import *
from pkpd.protocols import *
from pkpd.objects import PKPDDataSet
from .test_workflow import TestWorkflow
import copy

class TestIVIVCWorkflow2(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)

    def testIVIVCWorkflow(self):
        # Check that Simulation is working
        prot1 = self.newProtocol(ProtPKPDDissolutionSimulate,
                                objLabel='pkpd - simulate dissolution 0th order',
                                modelType=0,
                                parameters="1",
                                timeUnits=0,
                                resampleT=1,
                                resampleTF=100)
        self.launchProtocol(prot1)
        self.assertIsNotNone(prot1.outputExperiment.fnPKPD, "There was a problem")

        # Fit a 0th order dissolution
        print("Fitting 0th order ...")
        prot2 = self.newProtocol(ProtPKPDDissolutionFit,
                                 objLabel='pkpd - fit dissolution 0th order',
                                 predicted='A',
                                 globalSearch=True,
                                 modelType=0)
        prot2.inputExperiment.set(prot1.outputExperiment)
        self.launchProtocol(prot2)
        self.assertIsNotNone(prot2.outputExperiment.fnPKPD, "There was a problem")
        self.assertIsNotNone(prot2.outputFitting.fnFitting, "There was a problem")

        experiment = PKPDExperiment()
        experiment.load(prot2.outputExperiment.fnPKPD)
        K = float(experiment.samples['simulatedProfile'].descriptors['K'])
        self.assertTrue(K>0.99 and K<1.01)

        fitting = PKPDFitting()
        fitting.load(prot2.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Simulate PK
        print("Simulate PK IV ...")
        prot3 = self.newProtocol(ProtPKPDODESimulate,
                                objLabel='PK simulate IV',
                                odeSource=1,
                                viaType=0,
                                pkType=2,
                                prmUser='50, 2000, 1, 350, 10, 250', # Vmax, Km, Cl, V, Clp, Vp
                                doses='Bolus ; via=Intravenous; bolus; t=0 h; d=750000 ug',
                                tF=84*60,
                                sampling=2,
                                addIndividuals=True)
        self.launchProtocol(prot3)
        self.assertIsNotNone(prot3.outputExperiment.fnPKPD, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot3.outputExperiment.fnPKPD)
        AUC0t = float(experiment.samples['Simulation_0'].descriptors['AUC0t'])
        self.assertTrue(AUC0t > 735100 and AUC0t < 735200)
        Cmax = float(experiment.samples['Simulation_0'].descriptors['Cmax'])
        self.assertTrue(Cmax > 2125 and Cmax < 2135)

        # Fit a monocompartmental to IV
        print("Fitting model to IV ...")
        prot4 = self.newProtocol(ProtPKPDTwoCompartmentsClintCl,
                                 objLabel='pkpd - iv 2 compartments Cl+Clint',
                                 predicted='C',
                                 bounds='(0.0, 100.0); (0.0, 5000.0); (0.0, 2.0); (0.0, 500.0); (0.0, 20.0); (0.0, 500.0)',
                                 globalSearch=False)
        prot4.inputExperiment.set(prot3.outputExperiment)
        self.launchProtocol(prot4)
        self.assertIsNotNone(prot4.outputExperiment.fnPKPD, "There was a problem")
        self.assertIsNotNone(prot4.outputFitting.fnFitting, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot4.outputExperiment.fnPKPD)
        Vmax = float(experiment.samples['Simulation_0'].descriptors['Vmax'])
        Km = float(experiment.samples['Simulation_0'].descriptors['Km'])
        Cl = float(experiment.samples['Simulation_0'].descriptors['Cl'])
        V = float(experiment.samples['Simulation_0'].descriptors['V'])
        Clp = float(experiment.samples['Simulation_0'].descriptors['Clp'])
        Vp = float(experiment.samples['Simulation_0'].descriptors['Vp'])
        self.assertTrue(Vmax>49 and Vmax<51)
        self.assertTrue(Km>1995 and Km<2005)
        self.assertTrue(Cl>0.99 and Cl<1.01)
        self.assertTrue(V>340 and V<360)
        self.assertTrue(Clp>9.9 and Clp<10.1)
        self.assertTrue(Vp>240 and Vp<260)
        fitting = PKPDFitting()
        fitting.load(prot4.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Simulate PK F=1
        print("Simulate PK F=1 ...")
        prot3F = self.newProtocol(ProtPKPDODESimulate,
                                  objLabel='PK simulate Oral F=1',
                                  odeSource=1,
                                  viaType=1,
                                  viaPrm='1000, 0, 1', # Rin, tlag, F
                                  pkType=2,
                                  prmUser='50, 2000, 1, 350, 10, 250',
                                  doses='Bolus ; via=Oral; bolus; t=0 h; d=750000 ug',
                                  tF=84 * 60,
                                  sampling=2,
                                  addIndividuals=True)
        self.launchProtocol(prot3F)
        self.assertIsNotNone(prot3F.outputExperiment.fnPKPD, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot3F.outputExperiment.fnPKPD)
        AUC0t = float(experiment.samples['Simulation_0'].descriptors['AUC0t'])
        self.assertTrue(AUC0t > 734350 and AUC0t < 734450)
        Cmax = float(experiment.samples['Simulation_0'].descriptors['Cmax'])
        self.assertTrue(Cmax > 700 and Cmax < 710)

        # Fit a monocompartmental to F=1
        print("Fitting model to F=1 ...")
        prot4F = self.newProtocol(ProtPKPDTwoCompartmentsClintCl,
                                 objLabel='pkpd - ev0 2 compartments Cl+Clint F=1',
                                 predicted='C',
                                 bounds='(975,1025); (47.5, 52.5); (1975, 2025.0); (0.975, 1.025); (345, 355.0); (9.75, 10.25); (247.5, 252.5)',
                                 globalSearch=False)
        prot4F.inputExperiment.set(prot3F.outputExperiment)
        self.launchProtocol(prot4F)
        self.assertIsNotNone(prot4F.outputExperiment.fnPKPD, "There was a problem")
        self.assertIsNotNone(prot4F.outputFitting.fnFitting, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot4F.outputExperiment.fnPKPD)
        Rin = float(experiment.samples['Simulation_0'].descriptors['Oral_Rin'])
        Vmax = float(experiment.samples['Simulation_0'].descriptors['Vmax'])
        Km = float(experiment.samples['Simulation_0'].descriptors['Km'])
        Cl = float(experiment.samples['Simulation_0'].descriptors['Cl'])
        V = float(experiment.samples['Simulation_0'].descriptors['V'])
        Clp = float(experiment.samples['Simulation_0'].descriptors['Clp'])
        Vp = float(experiment.samples['Simulation_0'].descriptors['Vp'])
        self.assertTrue(Rin>990 and Rin<1010)
        self.assertTrue(Vmax>48 and Vmax<52) # 100
        self.assertTrue(Km>1950 and Km<2050) #
        self.assertTrue(Cl>0.99 and Cl<1.01)
        self.assertTrue(V>340 and V<360)
        self.assertTrue(Clp>9.7 and Clp<10.3)
        self.assertTrue(Vp>240 and Vp<260)
        fitting = PKPDFitting()
        fitting.load(prot4F.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Simulate PK F=0.5
        print("Simulate PK F=0.5 ...")
        prot3F2 = self.newProtocol(ProtPKPDODESimulate,
                                  objLabel='PK simulate Oral F=0.5',
                                  odeSource=1,
                                  viaType=1,
                                  viaPrm='1000, 0, 0.5', # Rin, tlag, F
                                  pkType=2,
                                  prmUser='50, 2000, 1, 350, 10, 250',
                                  doses='Bolus ; via=Oral; bolus; t=0 h; d=750000 ug',
                                  tF=84 * 60,
                                  sampling=2,
                                  addIndividuals=True)
        self.launchProtocol(prot3F2)
        self.assertIsNotNone(prot3F2.outputExperiment.fnPKPD, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot3F2.outputExperiment.fnPKPD)
        AUC0t = float(experiment.samples['Simulation_0'].descriptors['AUC0t'])
        self.assertTrue(AUC0t > 366200 and AUC0t < 367200)
        Cmax = float(experiment.samples['Simulation_0'].descriptors['Cmax'])
        self.assertTrue(Cmax > 460 and Cmax < 470)

        # Fit a monocompartmental to F=0.5
        print("Fitting model to F=0.5 ...")
        prot4F2 = self.newProtocol(ProtPKPDTwoCompartmentsClintCl,
                                  objLabel='pkpd - ev0 2 compartments Cl+Clint F=0.5',
                                  predicted='C',
                                  bounds='(975,1025); (47.5, 52.5); (1975, 2025.0); (0.975, 1.025); (345, 355.0); (9.75, 10.25); (247.5, 252.5)',
                                  globalSearch=False)
        prot4F2.inputExperiment.set(prot3F2.outputExperiment)
        self.launchProtocol(prot4F2)
        self.assertIsNotNone(prot4F2.outputExperiment.fnPKPD, "There was a problem")
        self.assertIsNotNone(prot4F2.outputFitting.fnFitting, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot4F2.outputExperiment.fnPKPD)
        Rin = float(experiment.samples['Simulation_0'].descriptors['Oral_Rin'])
        Vmax = float(experiment.samples['Simulation_0'].descriptors['Vmax'])
        Km = float(experiment.samples['Simulation_0'].descriptors['Km'])
        Cl = float(experiment.samples['Simulation_0'].descriptors['Cl'])
        V = float(experiment.samples['Simulation_0'].descriptors['V'])
        Clp = float(experiment.samples['Simulation_0'].descriptors['Clp'])
        Vp = float(experiment.samples['Simulation_0'].descriptors['Vp'])
        self.assertTrue(Rin>990 and Rin<1010)
        self.assertTrue(Vmax>48 and Vmax<52)
        self.assertTrue(Km>1950 and Km<2050)
        self.assertTrue(Cl>0.99 and Cl<1.01)
        self.assertTrue(V>340 and V<360)
        self.assertTrue(Clp>9.7 and Clp<10.3)
        self.assertTrue(Vp>240 and Vp<260)
        fitting = PKPDFitting()
        fitting.load(prot4F2.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        print("Change via to F=1 ...")
        prot3F3 = self.newProtocol(ProtPKPDChangeVia,
                                  objLabel='change via - F=1',
                                  viaName='Oral',
                                  newViaType="ev0",
                                  tlag=0,
                                  bioavailability=1)
        prot3F3.inputExperiment.set(prot3F2.outputExperiment)
        self.launchProtocol(prot3F3)
        self.assertIsNotNone(prot3F3.outputExperiment.fnPKPD, "There was a problem")

        # Fit a monocompartmental to F=0.5
        print("Fitting model to F=0.5 as if F=1...")
        prot4F3 = self.newProtocol(ProtPKPDTwoCompartmentsClintCl,
                                  objLabel='pkpd - ev0 2 compartments Cl+Clint F=0.5',
                                  predicted='C',
                                  bounds='(1980.0, 2020.0); (90.0, 110.0); (1900.0, 2100.0); (1.9, 2.1); (700.0, 800.0); (19.0, 20.0); (450.0, 550.0)',
                                  globalSearch=False)
        prot4F3.inputExperiment.set(prot3F3.outputExperiment)
        self.launchProtocol(prot4F3)
        self.assertIsNotNone(prot4F3.outputExperiment.fnPKPD, "There was a problem")
        self.assertIsNotNone(prot4F3.outputFitting.fnFitting, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot4F3.outputExperiment.fnPKPD)
        Rin = float(experiment.samples['Simulation_0'].descriptors['Oral_Rin'])
        Vmax = float(experiment.samples['Simulation_0'].descriptors['Vmax'])
        Km = float(experiment.samples['Simulation_0'].descriptors['Km'])
        Cl = float(experiment.samples['Simulation_0'].descriptors['Cl'])
        V = float(experiment.samples['Simulation_0'].descriptors['V'])
        Clp = float(experiment.samples['Simulation_0'].descriptors['Clp'])
        Vp = float(experiment.samples['Simulation_0'].descriptors['Vp'])
        self.assertTrue(Rin>1980 and Rin<2020)
        self.assertTrue(Vmax>90 and Vmax<110)
        self.assertTrue(Km>1950 and Km<2050)
        self.assertTrue(Cl>1.9 and Cl<2.1)
        self.assertTrue(V>730 and V<770)
        self.assertTrue(Clp>18.5 and Clp<21.5)
        self.assertTrue(Vp>450 and Vp<550)
        fitting = PKPDFitting()
        fitting.load(prot4F2.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        print("Change via to F=unk ...")
        prot3Fu = self.newProtocol(ProtPKPDChangeVia,
                                   objLabel='change via - F=unk',
                                   viaName='Oral',
                                   newViaType="ev0",
                                   tlag=0)
        prot3Fu.inputExperiment.set(prot3F2.outputExperiment)
        self.launchProtocol(prot3Fu)
        self.assertIsNotNone(prot3Fu.outputExperiment.fnPKPD, "There was a problem")

        # Fit a monocompartmental to F=unk
        print("Fitting model to F=unk ...")
        prot4Fu = self.newProtocol(ProtPKPDTwoCompartmentsClintCl,
                                   objLabel='pkpd - ev0 2 compartments Cl+Clint F=unk',
                                   predicted='C',
                                   bounds='(0.0, 1.0); (975.0, 1025.0); (47.5, 52.5); (1975.0, 2025.0); (0.975, 1.025); (345.0, 355.0); (9.75, 10.25); (247.5, 252.5)',
                                   globalSearch=False)
        prot4Fu.inputExperiment.set(prot3Fu.outputExperiment)
        self.launchProtocol(prot4Fu)
        self.assertIsNotNone(prot4Fu.outputExperiment.fnPKPD, "There was a problem")
        self.assertIsNotNone(prot4Fu.outputFitting.fnFitting, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot4Fu.outputExperiment.fnPKPD)
        F = float(experiment.samples['Simulation_0'].descriptors['Oral_bioavailability'])
        Rin = float(experiment.samples['Simulation_0'].descriptors['Oral_Rin'])
        Vmax = float(experiment.samples['Simulation_0'].descriptors['Vmax'])
        Km = float(experiment.samples['Simulation_0'].descriptors['Km'])
        Cl = float(experiment.samples['Simulation_0'].descriptors['Cl'])
        V = float(experiment.samples['Simulation_0'].descriptors['V'])
        Clp = float(experiment.samples['Simulation_0'].descriptors['Clp'])
        Vp = float(experiment.samples['Simulation_0'].descriptors['Vp'])
        self.assertTrue(F > 0.49 and F < 0.51)
        self.assertTrue(Rin>980 and Rin<1020)
        self.assertTrue(Vmax>48 and Vmax<52)
        self.assertTrue(Km>1950 and Km<2050)
        self.assertTrue(Cl>0.99 and Cl<1.01)
        self.assertTrue(V>340 and V<360)
        self.assertTrue(Clp>9.6 and Clp<10.4)
        self.assertTrue(Vp>230 and Vp<270)
        fitting = PKPDFitting()
        fitting.load(prot4Fu.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2 > 0.99)

        print("Deconvolving Fourier ...")
        protDeconvFu = self.newProtocol(ProtPKPDDeconvolveFourier,
                                       objLabel='pkpd - deconvolve Fourier',
                                       externalIV=1,
                                       normalize=True,
                                       saturate=True,
                                       considerBioaval=1)
        protDeconvFu.inputODE.set(prot4Fu)
        protDeconvFu.externalIVODE.set(prot4)
        self.launchProtocol(protDeconvFu)
        self.assertIsNotNone(protDeconvFu.outputExperiment.fnPKPD, "There was a problem")

        print("Deconvolving Fourier amount...")
        protDeconvFu2 = self.newProtocol(ProtPKPDDeconvolveFourier,
                                       objLabel='pkpd - deconvolve Fourier amount',
                                       externalIV=1,
                                       normalize=False,
                                       considerBioaval=2)
        protDeconvFu2.inputODE.set(prot4Fu)
        protDeconvFu2.externalIVODE.set(prot4)
        self.launchProtocol(protDeconvFu2)
        self.assertIsNotNone(protDeconvFu2.outputExperiment.fnPKPD, "There was a problem")

        print("IVIVC F=unk...")
        protIVIVC = self.newProtocol(ProtPKPDDissolutionIVIVC,
                                     objLabel='pkpd - ivivc unk',
                                     timeScale=2,
                                     responseScale=1,
                                     ABounds="[0.1,10]")
        protIVIVC.inputInVitro.set(prot2)
        protIVIVC.inputInVivo.set(protDeconvFu)
        self.launchProtocol(protIVIVC)
        self.assertIsNotNone(protIVIVC.outputExperimentFabs.fnPKPD, "There was a problem")

        experiment = PKPDExperiment()
        experiment.load(protIVIVC.outputExperimentFabs.fnPKPD)
        R = float(experiment.samples['ivivc_0001'].descriptors['R'])
        self.assertTrue(R > 0.99)

        fhSummary = open(protIVIVC._getPath("summary.txt"))
        i=0
        for line in fhSummary.readlines():
            token = (line.split()[1]).split('=')[1]
            token = token[:-1]
            value = float(token)
            if i==0:
                k=value
                self.assertTrue(k>0.26 and k<0.28)
            elif i==1:
                A=value
                self.assertTrue(A>0.49 and A<0.51)
                break
            i+=1
        fhSummary.close()

        # Simulate PK with IVIVC
        print("Simulate PK")
        protSimulatePK = self.newProtocol(ProtPKPDDissolutionPKSimulation,
                                          objLabel='pkpd - simulate PK',
                                          ignorePKbioavailability=True,
                                          inputDose=750000,
                                          inputN=1,
                                          tF=84,
                                          addIndividuals=True)
        protSimulatePK.inputInVitro.set(prot2.outputFitting)
        protSimulatePK.inputPK.set(prot4Fu.outputFitting)
        protSimulatePK.inputIvIvC.set(protIVIVC.outputExperimentFabs)
        self.launchProtocol(protSimulatePK)
        self.assertIsNotNone(protSimulatePK.outputExperiment.fnPKPD, "There was a problem")

        # Internal validity
        print("Internal validity ...")
        protInternal = self.newProtocol(ProtPKPDIVIVCInternalValidity,
                                        objLabel='pkpd - internal validity')
        protInternal.inputExperiment.set(prot3F2.outputExperiment)
        protInternal.inputSimulated.set(protSimulatePK.outputExperiment)
        self.launchProtocol(protInternal)

        fnSummary = protInternal._getPath('summary.txt')
        self.assertTrue(os.path.exists(fnSummary))
        fh = open(fnSummary, 'r')
        for line in fh.readlines():
            tokens = line.split('=')
            self.assertTrue(abs(float(tokens[-1])) < 5)
        fh.close()

if __name__ == "__main__":
    unittest.main()
