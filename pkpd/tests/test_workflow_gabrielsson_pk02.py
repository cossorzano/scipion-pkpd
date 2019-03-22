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


class TestGabrielssonPK02Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK02')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK02Workflow(self):
        #First, import an experiment

        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)      

        # Filter time variable
        print "Filter time"
        protFilterTime = self.newProtocol(ProtPKPDFilterMeasurements,
                                                  objLabel='pkpd - filter measurements t>50',
                                                  filterType=1, condition='$(t)>50 and $(t)<=300')
        protFilterTime.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protFilterTime)
        self.assertIsNotNone(protFilterTime.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('protFilterTime', protFilterTime)

        # Filter concentration variable
        print "Filter Cp"
        protFilterCp = self.newProtocol(ProtPKPDFilterMeasurements,
                                        objLabel='pkpd - filter measurements Cp>0',
                                        filterType=0, condition='$(Cp)<=0')
        protFilterCp.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protFilterCp)
        self.assertIsNotNone(protFilterCp.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('ProtFilterCp', protFilterCp)

        # Fit a single exponential to the input data_test
        print "Fitting an exponential..."
        protEliminationRate = self.newProtocol(ProtPKPDEliminationRate,
                                               objLabel='pkpd - elimination rate',
                                               predictor='t', predicted='Cp')
        protEliminationRate.inputExperiment.set(protFilterTime.outputExperiment)
        self.launchProtocol(protEliminationRate)
        self.assertIsNotNone(protEliminationRate.outputExperiment.fnPKPD, "There was a problem with the exponential fitting")
        self.assertIsNotNone(protEliminationRate.outputFitting.fnFitting, "There was a problem with the exponential fitting")
        self.validateFiles('protEliminationRate', protEliminationRate)
        experiment = PKPDExperiment()
        experiment.load(protEliminationRate.outputExperiment.fnPKPD)
        c1 = float(experiment.samples['Individual1'].descriptors['c1'])
        lambda1 = float(experiment.samples['Individual1'].descriptors['lambda1'])
        self.assertAlmostEqual(c1,4.14,2)
        self.assertAlmostEqual(lambda1,0.00888,3) # Gabrielsson p 509: K=0.01

        fitting = PKPDFitting()
        fitting.load(protEliminationRate.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.93)
        self.assertTrue(fitting.sampleFits[0].AIC<-11)

        # Estimate absorption rate
        print "Estimation of the absorption rate..."
        protAbsorptionRate = self.newProtocol(ProtPKPDAbsorptionRate,
                                              objLabel='pkpd - absorption rate',
                                              bounds='(0,0.1);(10,60);(0,25)')
        protAbsorptionRate.inputExperiment.set(protFilterCp.outputExperiment)
        protAbsorptionRate.protElimination.set(protEliminationRate)
        self.launchProtocol(protAbsorptionRate)
        self.assertIsNotNone(protAbsorptionRate.outputExperiment.fnPKPD, "There was a problem with the absorption rate estimation ")
        self.assertIsNotNone(protAbsorptionRate.outputFitting.fnFitting, "There was a problem with the absorption rate estimation ")
        self.validateFiles('protAbsorptionRate', protAbsorptionRate)

        experiment = PKPDExperiment()
        experiment.load(protAbsorptionRate.outputExperiment.fnPKPD)
        Cmax = float(experiment.samples['Individual1'].descriptors['Cmax'])
        Ka = float(experiment.samples['Individual1'].descriptors['Ka'])
        Ke = float(experiment.samples['Individual1'].descriptors['Ke'])
        Vd = float(experiment.samples['Individual1'].descriptors['Vd'])
        tlag = float(experiment.samples['Individual1'].descriptors['tlag'])
        tmax = float(experiment.samples['Individual1'].descriptors['tmax'])
        self.assertTrue(Cmax>1.7)
        self.assertTrue(Ka>0.03 and Ka<0.06) # Gabrielsson p 514, Solution II: K01: 0.0428
        self.assertTrue(Ke<0.01)
        self.assertTrue(Vd>20)
        self.assertTrue(tlag>10 and tlag<23)
        self.assertTrue(tmax>30 and tmax<70)

        fitting = PKPDFitting()
        fitting.load(protAbsorptionRate.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.96)
        self.assertTrue(fitting.sampleFits[0].AIC<-17)

        # Fit a monocompartmental model with first order absorption
        print "Fitting monocompartmental model..."
        protEV1MonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                  objLabel='pkpd - ev1 monocompartment',
                                                  bounds='(0.0, 30.0); (0.0, 0.2); (0.0, 1.0); (0.0, 100.0)')
        protEV1MonoCompartment.inputExperiment.set(protAbsorptionRate.outputExperiment)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)

        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual1'].descriptors['Cl'])
        V = float(experiment.samples['Individual1'].descriptors['V'])
        Ka = float(experiment.samples['Individual1'].descriptors['Oral_Ka'])
        tlag = float(experiment.samples['Individual1'].descriptors['Oral_tlag'])
        self.assertTrue(Cl>0.27 and Cl<0.31) # Gabrielsson p 515, Solution II: CL/F=0.2819
        self.assertTrue(V>10 and V<40) # Gabrielsson p 515, Solution II: V/F=32.05 -------------- Mine: 27.5
        self.assertTrue(Ka>0.025 and Ka<0.06) # Gabrielsson p 511, Solution II: Ka=0.043 -------- Mine: 0.0264
        self.assertTrue(tlag>10 and tlag<25) # Gabrielsson p 511, Solution II: tlag=16

        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.97)
        self.assertTrue(fitting.sampleFits[0].AIC<-15)

        # Population bootstrap
        print "Bootstrapping model..."
        protBootstrap = self.newProtocol(ProtPKPDODEBootstrap,
                                         objLabel='pkpd - bootstrap')
        protBootstrap.inputODE.set(protEV1MonoCompartment)
        self.launchProtocol(protBootstrap)
        self.assertIsNotNone(protBootstrap.outputPopulation.fnFitting, "There was a problem with the bootstrap")
        self.validateFiles('protBootstrap', protBootstrap)

        # Filter population
        print "Filtering bootstrap..."
        protFilterBootstrap = self.newProtocol(ProtPKPDFilterPopulation,
                                               objLabel='pkpd - filter population',
                                               condition="$(R2)<0.95")
        protFilterBootstrap.inputPopulation.set(protBootstrap.outputPopulation)
        self.launchProtocol(protFilterBootstrap)
        self.assertIsNotNone(protFilterBootstrap.outputPopulation.fnFitting, "There was a problem with the population filter")
        self.validateFiles('protFilterBootstrap', protFilterBootstrap)

if __name__ == "__main__":
    unittest.main()
