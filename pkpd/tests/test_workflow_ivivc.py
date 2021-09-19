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
                                 objLabel='pkpd - iv monocompartment',
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
if __name__ == "__main__":
    unittest.main()
