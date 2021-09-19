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

class TestIVIVCWorkflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)

    def testIVIVCWorkflow(self):
        # Check that Simulation is working
        prot1 = self.newProtocol(ProtPKPDDissolutionSimulate,
                                objLabel='pkpd - simulate dissolution 1st order',
                                modelType=1,
                                parameters="100;0.01",
                                timeUnits=0,
                                resampleT=10,
                                resampleTF=800)
        self.launchProtocol(prot1)
        self.assertIsNotNone(prot1.outputExperiment.fnPKPD, "There was a problem")

        # Fit a first order dissolution
        print("Fitting 1st order ...")
        prot2 = self.newProtocol(ProtPKPDDissolutionFit,
                                 objLabel='pkpd - fit dissolution 1st order',
                                 predicted='A',
                                 globalSearch=True,
                                 modelType=1)
        prot2.inputExperiment.set(prot1.outputExperiment)
        self.launchProtocol(prot2)
        self.assertIsNotNone(prot2.outputExperiment.fnPKPD, "There was a problem")
        self.assertIsNotNone(prot2.outputFitting.fnFitting, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot2.outputExperiment.fnPKPD)
        Vmax = float(experiment.samples['simulatedProfile'].descriptors['Vmax'])
        self.assertTrue(Vmax>99.9)
        beta = float(experiment.samples['simulatedProfile'].descriptors['beta'])
        self.assertTrue(beta>0.0099 and beta<0.0101)

        fitting = PKPDFitting()
        fitting.load(prot2.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Simulate PK
        print("Simulate PK IV ...")
        prot3 = self.newProtocol(ProtPKPDODESimulate,
                                objLabel='PK simulate IV',
                                odeSource=1,
                                viaType=0,
                                pkType=0,
                                prmUser='0.0025, 1',
                                doses='Bolus ; via=Intravenous; bolus; t=0 h; d=1000 ug',
                                tF=36*60,
                                sampling=60,
                                addIndividuals=True)
        self.launchProtocol(prot3)
        self.assertIsNotNone(prot3.outputExperiment.fnPKPD, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot3.outputExperiment.fnPKPD)
        AUC0t = float(experiment.samples['Simulation_0'].descriptors['AUC0t'])
        self.assertTrue(AUC0t > 397980 and AUC0t < 397984)
        Cavg = float(experiment.samples['Simulation_0'].descriptors['Cavg'])
        self.assertTrue(Cavg > 184 and Cavg < 185)

        # Fit a monocompartmental to IV
        print("Fitting monocompartmental model to IV ...")
        prot4 = self.newProtocol(ProtPKPDMonoCompartment,
                                 objLabel='pkpd - iv monocompartment',
                                 predicted='C',
                                 bounds='(0,0.1);(0,5)')
        prot4.inputExperiment.set(prot3.outputExperiment)
        self.launchProtocol(prot4)
        self.assertIsNotNone(prot4.outputExperiment.fnPKPD, "There was a problem")
        self.assertIsNotNone(prot4.outputFitting.fnFitting, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot4.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Simulation_0'].descriptors['Cl'])
        V = float(experiment.samples['Simulation_0'].descriptors['V'])
        self.assertTrue(Cl>0.00249 and Cl<0.00251)
        self.assertTrue(V>0.99 and V<1.01)
        fitting = PKPDFitting()
        fitting.load(prot4.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Simulate PK Fast
        print("Simulate PK Fast ...")
        prot3F = self.newProtocol(ProtPKPDODESimulate,
                                  objLabel='PK simulate Oral fast',
                                  odeSource=1,
                                  viaType=2,
                                  viaPrm='0.005, 0, 1',
                                  pkType=0,
                                  prmUser='0.0025, 1',
                                  doses='Bolus ; via=Oral; bolus; t=0 h; d=1000 ug',
                                  tF=36 * 60,
                                  sampling=5,
                                  addIndividuals=True)
        self.launchProtocol(prot3F)
        self.assertIsNotNone(prot3F.outputExperiment.fnPKPD, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot3F.outputExperiment.fnPKPD)
        AUC0t = float(experiment.samples['Simulation_0'].descriptors['AUC0t'])
        self.assertTrue(AUC0t > 396380 and AUC0t < 396400)
        Cavg = float(experiment.samples['Simulation_0'].descriptors['Cavg'])
        self.assertTrue(Cavg > 183 and Cavg < 184)

        # Fit a monocompartmental to Fast
        print("Fitting monocompartmental model to Fast ...")
        prot4F = self.newProtocol(ProtPKPDMonoCompartment,
                                 objLabel='pkpd - fast monocompartment',
                                 predicted='C',
                                 bounds='(0,0.1);(0,0.1);(0,5)')
        prot4F.inputExperiment.set(prot3F.outputExperiment)
        self.launchProtocol(prot4F)
        self.assertIsNotNone(prot4F.outputExperiment.fnPKPD, "There was a problem")
        self.assertIsNotNone(prot4F.outputFitting.fnFitting, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot4F.outputExperiment.fnPKPD)
        Ka = float(experiment.samples['Simulation_0'].descriptors['Oral_Ka'])
        Cl = float(experiment.samples['Simulation_0'].descriptors['Cl'])
        V = float(experiment.samples['Simulation_0'].descriptors['V'])
        self.assertTrue(Ka > 0.0048 and Ka < 0.0052)
        self.assertTrue(Cl > 0.00248 and Cl < 0.00252)
        self.assertTrue(V > 0.99 and V < 1.01)
        fitting = PKPDFitting()
        fitting.load(prot4F.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2 > 0.99)

        # Simulate PK Slow
        print("Simulate PK Slow ...")
        prot3S = self.newProtocol(ProtPKPDODESimulate,
                                  objLabel='PK simulate Oral slow',
                                  odeSource=1,
                                  viaType=2,
                                  viaPrm='0.00125, 0, 1',
                                  pkType=0,
                                  prmUser='0.0025, 1',
                                  doses='Bolus ; via=Oral; bolus; t=0 h; d=1000 ug',
                                  tF=84 * 60,
                                  sampling=5,
                                  addIndividuals=True)
        self.launchProtocol(prot3S)
        self.assertIsNotNone(prot3S.outputExperiment.fnPKPD, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot3S.outputExperiment.fnPKPD)
        AUC0t = float(experiment.samples['Simulation_0'].descriptors['AUC0t'])
        self.assertTrue(AUC0t > 398520 and AUC0t < 398540)
        Cavg = float(experiment.samples['Simulation_0'].descriptors['Cavg'])
        self.assertTrue(Cavg > 78.5 and Cavg < 79.5)

        # Fit a monocompartmental to Slow
        print("Fitting monocompartmental model to Slow ...")
        prot4S = self.newProtocol(ProtPKPDMonoCompartment,
                                  objLabel='pkpd - slow monocompartment',
                                  predicted='C',
                                  bounds='(0,0.1);(0,0.1);(0.95,1.05)')
        prot4S.inputExperiment.set(prot3S.outputExperiment)
        self.launchProtocol(prot4S)
        self.assertIsNotNone(prot4S.outputExperiment.fnPKPD, "There was a problem")
        self.assertIsNotNone(prot4S.outputFitting.fnFitting, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot4S.outputExperiment.fnPKPD)
        Ka = float(experiment.samples['Simulation_0'].descriptors['Oral_Ka'])
        Cl = float(experiment.samples['Simulation_0'].descriptors['Cl'])
        V = float(experiment.samples['Simulation_0'].descriptors['V'])
        self.assertTrue(Ka > 0.00115 and Ka < 0.00135)
        self.assertTrue(Cl > 0.00248 and Cl < 0.00252)
        self.assertTrue(V > 0.95 and V < 1.05)
        fitting = PKPDFitting()
        fitting.load(prot4S.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2 > 0.99)

        # Simulate PK Bioavail
        print("Simulate PK Bioavail ...")
        prot3B = self.newProtocol(ProtPKPDODESimulate,
                                  objLabel='PK simulate Oral bioavail',
                                  odeSource=1,
                                  viaType=2,
                                  viaPrm='0.005, 0, 0.5',
                                  pkType=0,
                                  prmUser='0.0025, 1',
                                  doses='Bolus ; via=Oral; bolus; t=0 h; d=1000 ug',
                                  tF=36 * 60,
                                  sampling=5,
                                  addIndividuals=True)
        self.launchProtocol(prot3B)
        self.assertIsNotNone(prot3B.outputExperiment.fnPKPD, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot3B.outputExperiment.fnPKPD)
        AUC0t = float(experiment.samples['Simulation_0'].descriptors['AUC0t'])
        self.assertTrue(AUC0t > 198180 and AUC0t < 198210)
        Cavg = float(experiment.samples['Simulation_0'].descriptors['Cavg'])
        self.assertTrue(Cavg > 91.3 and Cavg < 92.3)

        print("Change via to unknown bioavailability ...")
        prot3B2 = self.newProtocol(ProtPKPDChangeVia,
                                  objLabel='change via - unknown bioavail',
                                  viaName='Oral',
                                  newViaType="ev1",
                                  tlag=0)
        prot3B2.inputExperiment.set(prot3B.outputExperiment)
        self.launchProtocol(prot3B2)
        self.assertIsNotNone(prot3B2.outputExperiment.fnPKPD, "There was a problem")

        # Fit a monocompartmental to Bioavail
        print("Fitting monocompartmental model to Bioavail ...")
        prot4B = self.newProtocol(ProtPKPDMonoCompartment,
                                  objLabel='pkpd - bioavail monocompartment',
                                  predicted='C',
                                  bounds='(0,1);(0,0.1);(0,0.003);(0.95,1.05)')
        prot4B.inputExperiment.set(prot3B2.outputExperiment)
        self.launchProtocol(prot4B)
        self.assertIsNotNone(prot4B.outputExperiment.fnPKPD, "There was a problem")
        self.assertIsNotNone(prot4B.outputFitting.fnFitting, "There was a problem")
        experiment = PKPDExperiment()
        experiment.load(prot4B.outputExperiment.fnPKPD)
        Ka = float(experiment.samples['Simulation_0'].descriptors['Oral_Ka'])
        F = float(experiment.samples['Simulation_0'].descriptors['Oral_bioavailability'])
        Cl = float(experiment.samples['Simulation_0'].descriptors['Cl'])
        V = float(experiment.samples['Simulation_0'].descriptors['V'])
        self.assertTrue(Ka > 0.0048 and Ka < 0.0052)
        self.assertTrue(F > 0.47 and F < 0.53)
        self.assertTrue(Cl > 0.0023 and Cl < 0.0027)
        self.assertTrue(V > 0.95 and V < 1.05)
        fitting = PKPDFitting()
        fitting.load(prot4B.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2 > 0.99)

        def checkDeconvolution(msg, protDeconv, refKa, refAmax, tol):
            print("Fitting 1st order to deconvolve %s ..."%msg)
            protFit = self.newProtocol(ProtPKPDDissolutionFit,
                                     objLabel='pkpd - fit dissolution 1st order %s'%msg,
                                     predicted='A',
                                     globalSearch=True,
                                     modelType=1)
            protFit.inputExperiment.set(protDeconv.outputExperiment)
            self.launchProtocol(protFit)
            self.assertIsNotNone(protFit.outputExperiment.fnPKPD, "There was a problem")
            self.assertIsNotNone(protFit.outputFitting.fnFitting, "There was a problem")
            experiment = PKPDExperiment()
            experiment.load(protFit.outputExperiment.fnPKPD)
            Vmax = float(experiment.samples['Simulation_0'].descriptors['Vmax'])
            self.assertTrue(Vmax>(1-tol)*refAmax and Vmax<(1+tol)*refAmax)
            beta = float(experiment.samples['Simulation_0'].descriptors['beta'])
            self.assertTrue(beta>(1-tol)*refKa and beta<(1+tol)*refKa)

            fitting = PKPDFitting()
            fitting.load(protFit.outputFitting.fnFitting)
            self.assertTrue(fitting.sampleFits[0].R2>0.99)

        def ivivc(msg, protVitro, protDeconv):
            print("IVIVC %s ..."%msg)
            protIVIVC = self.newProtocol(ProtPKPDDissolutionIVIVC,
                                         objLabel='pkpd - ivivc %s'%msg,
                                         timeScale=2)
            protIVIVC.inputInVitro.set(protVitro)
            protIVIVC.inputInVivo.set(protDeconv)
            self.launchProtocol(protIVIVC)
            self.assertIsNotNone(protIVIVC.outputExperimentFabs.fnPKPD, "There was a problem")
            return protIVIVC

        def deconvolve(msg, protVitro, protPK, protPKIV, refKa, refAmax, tol):
            print("Deconvolving %s ..."%msg)
            protDeconv = self.newProtocol(ProtPKPDDeconvolve,
                                      objLabel='pkpd - deconvolve %s'%msg,
                                      saturate=False)
            protDeconv.inputODE.set(protPK)
            self.launchProtocol(protDeconv)
            self.assertIsNotNone(protDeconv.outputExperiment.fnPKPD, "There was a problem")
            checkDeconvolution(msg, protDeconv, refKa, refAmax, tol)
            ivivc(msg, protVitro, protDeconv)

            print("Deconvolving Fourier %s ..." % msg)
            protDeconvF = self.newProtocol(ProtPKPDDeconvolveFourier,
                                          objLabel='pkpd - deconvolve Fourier %s' % msg,
                                          saturate=False)
            protDeconvF.inputODE.set(protPK)
            self.launchProtocol(protDeconvF)
            self.assertIsNotNone(protDeconvF.outputExperiment.fnPKPD, "There was a problem")
            checkDeconvolution("Fourier "+msg, protDeconvF, refKa, refAmax, tol)
            ivivc("Fourier "+msg, protVitro, protDeconvF)

            print("Deconvolving Fourier IV %s ..." % msg)
            protDeconvFIV = self.newProtocol(ProtPKPDDeconvolveFourier,
                                          objLabel='pkpd - deconvolve Fourier IV %s' % msg,
                                          externalIV=1,
                                          saturate=False)
            protDeconvFIV.inputODE.set(protPK)
            protDeconvFIV.externalIVODE.set(protPKIV)
            self.launchProtocol(protDeconvFIV)
            self.assertIsNotNone(protDeconvFIV.outputExperiment.fnPKPD, "There was a problem")
            checkDeconvolution("Fourier IV "+msg, protDeconvFIV, refKa, refAmax, tol)
            ivivc("Fourier IV "+msg, protVitro, protDeconvFIV)

            print("Deconvolving Wagner Nelson %s ..." % msg)
            protDeconvW = self.newProtocol(ProtPKPDDeconvolutionWagnerNelson,
                                           objLabel='pkpd - deconvolve Wagner %s' % msg,
                                           concVar='C',
                                           saturate=False)
            protDeconvW.inputExperiment.set(protPK.outputExperiment)
            self.launchProtocol(protDeconvW)
            self.assertIsNotNone(protDeconvW.outputExperiment.fnPKPD, "There was a problem")
            checkDeconvolution("Wagner " + msg, protDeconvW, refKa, refAmax, tol)
            ivivc("Wagner "+msg, protVitro, protDeconvW)

            print("Deconvolving Wagner Nelson IV %s ..." % msg)
            protDeconvWIV = self.newProtocol(ProtPKPDDeconvolutionWagnerNelson,
                                           objLabel='pkpd - deconvolve Wagner IV %s' % msg,
                                           externalIV=1,
                                           concVar='C',
                                           saturate=False)
            protDeconvWIV.inputExperiment.set(protPK.outputExperiment)
            protDeconvWIV.externalIVODE.set(protPKIV)
            self.launchProtocol(protDeconvWIV)
            self.assertIsNotNone(protDeconvWIV.outputExperiment.fnPKPD, "There was a problem")
            checkDeconvolution("WagnerIV " + msg, protDeconvWIV, refKa, refAmax, tol)
            ivivc("WagnerIV "+msg, protVitro, protDeconvWIV)

            return [protDeconv, protDeconvF, protDeconvW, protDeconvWIV]

        [protDeconvF, protDeconvFF, protDeconvWF, protDeconvWIVF] = deconvolve("fast", prot2, prot4F, prot4, 0.005, 100, 0.03)
        [protDeconvS, protDeconvFS, protDeconvWS, protDeconvWIVS] = deconvolve("slow", prot2, prot4S, prot4, 0.00125, 100, 0.03)
        [protDeconvB, protDeconvFB, protDeconvWB, protDeconvWIVB] = deconvolve("bioaval", prot2, prot4B, prot4, 0.005, 50, 0.06)

if __name__ == "__main__":
    unittest.main()
