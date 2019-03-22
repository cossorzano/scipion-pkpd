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


import unittest, sys
from pyworkflow.em import *
from pyworkflow.tests import *
from pkpd.protocols import *
from test_workflow import TestWorkflow
from pkpd.objects import PKPDDataSet


class TestGabrielssonPK01Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK01')
        cls.exptFn = cls.dataset.getFile('experiment')

    
    def testGabrielssonPK01Workflow(self):
        #First, import an experiment

        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)      

        # Change the concentration unit to mg/L
        print "Change Units"
        protChangeUnits = self.newProtocol(ProtPKPDChangeUnits,
                                           objLabel='pkpd - change units',
                                           labelToChange='Cp', newUnitsCategory=4,
                                           newUnitsCategoryConc=1)
        protChangeUnits.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protChangeUnits)
        self.assertIsNotNone(protChangeUnits.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeUnits)

        # Fit a single exponential to the input data_test
        print "Compute elimination rate..."
        protEliminationRate = self.newProtocol(ProtPKPDEliminationRate,
                                               objLabel='pkpd - elimination rate',
                                               predictor='t', predicted='Cp')
        protEliminationRate.inputExperiment.set(protChangeUnits.outputExperiment)
        self.launchProtocol(protEliminationRate)
        self.assertIsNotNone(protEliminationRate.outputExperiment.fnPKPD, "There was a problem with the exponential fitting")
        self.assertIsNotNone(protEliminationRate.outputFitting.fnFitting, "There was a problem with the exponential fitting")
        self.validateFiles('protEliminationRate', protEliminationRate)
        experiment = PKPDExperiment()
        experiment.load(protEliminationRate.outputExperiment.fnPKPD)
        self.assertAlmostEqual(float(experiment.samples['Individual1'].descriptors['c1']),1.011,3)
        self.assertAlmostEqual(float(experiment.samples['Individual1'].descriptors['lambda1']),0.0104,3)
        fitting = PKPDFitting()
        fitting.load(protEliminationRate.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.988)
        self.assertTrue(fitting.sampleFits[0].AIC<-45.8)

        # Non-compartmental analysis
        print "Performing Non-compartmental analysis..."
        protNCAIVObs = self.newProtocol(ProtPKPDNCAIVObs,
                                        objLabel='pkpd - nca iv observations')
        protNCAIVObs.inputExperiment.set(protChangeUnits.outputExperiment)
        protNCAIVObs.protElimination.set(protEliminationRate)
        self.launchProtocol(protNCAIVObs)
        self.assertIsNotNone(protNCAIVObs.outputExperiment.fnPKPD, "There was a problem with the Non-compartmental analysis ")
        self.assertIsNotNone(protNCAIVObs.outputAnalysis.fnAnalysis, "There was a problem with the Non-compartmental analysis ")
        self.validateFiles('protNCAIVObs', protNCAIVObs)

        experiment = PKPDExperiment()
        experiment.load(protNCAIVObs.outputExperiment.fnPKPD)
        AUC_0inf = float(experiment.samples['Individual1'].descriptors['AUC_0inf'])
        AUC_0t = float(experiment.samples['Individual1'].descriptors['AUC_0t'])
        AUMC_0inf = float(experiment.samples['Individual1'].descriptors['AUMC_0inf'])
        AUMC_0t = float(experiment.samples['Individual1'].descriptors['AUMC_0t'])
        CL_0inf = float(experiment.samples['Individual1'].descriptors['CL_0inf'])
        CL_0t = float(experiment.samples['Individual1'].descriptors['CL_0t'])
        MRT = float(experiment.samples['Individual1'].descriptors['MRT'])
        Vd_0inf = float(experiment.samples['Individual1'].descriptors['Vd_0inf'])
        Vd_0t = float(experiment.samples['Individual1'].descriptors['Vd_0t'])
        Vss = float(experiment.samples['Individual1'].descriptors['Vss'])
        thalf = float(experiment.samples['Individual1'].descriptors['thalf'])
        self.assertAlmostEqual(AUC_0inf,86.4555,3) # Gabrielsson p 495: 97.7
        self.assertAlmostEqual(AUC_0t,67.3003,3)
        self.assertAlmostEqual(AUMC_0inf,9013.1394,3)
        self.assertAlmostEqual(AUMC_0t,4305.2328,3)
        self.assertAlmostEqual(CL_0inf,0.1156,3)
        self.assertAlmostEqual(CL_0t,0.1485,3)
        self.assertAlmostEqual(MRT,104.2516,3) # Gabrielsson p 495: 97.5
        self.assertAlmostEqual(Vd_0inf,11.078,3)
        self.assertAlmostEqual(Vd_0t,14.2311,3)
        self.assertAlmostEqual(Vss,12.0584,3)
        self.assertAlmostEqual(thalf,66.387,3) # Gabrielsson p 495: 67.6

        # Fit a monocompartmental model
        print "Fitting monocompartmental model..."
        protIVMonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                 objLabel='pkpd - iv monocompartment',
                                                 bounds='(0.0, 0.2); (0.0, 20.0)')
        protIVMonoCompartment.inputExperiment.set(protNCAIVObs.outputExperiment)
        self.launchProtocol(protIVMonoCompartment)
        self.assertIsNotNone(protIVMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protIVMonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('ProtIVMonoCompartment', protIVMonoCompartment)

        experiment = PKPDExperiment()
        experiment.load(protIVMonoCompartment.outputExperiment.fnPKPD)
        Cl=float(experiment.samples['Individual1'].descriptors['Cl'])
        V=float(experiment.samples['Individual1'].descriptors['V'])
        self.assertTrue(Cl>0.09 and Cl<0.11) # Gabrielsson, p 495: 0.10
        self.assertTrue(V>9.7 and V<10) # Gabrielsson p 495: 9.98
        fitting = PKPDFitting()
        fitting.load(protIVMonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.9887)
        self.assertTrue(fitting.sampleFits[0].AIC<-45.8)

if __name__ == "__main__":
    unittest.main()
