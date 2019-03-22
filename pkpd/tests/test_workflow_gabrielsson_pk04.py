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


class TestGabrielssonPK04Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK04')
        cls.exptFn = cls.dataset.getFile('experiment')

    
    def testGabrielssonPK04Workflow(self):
        #First, import an experiment
        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Change the time unit to minute
        print "Change Units"
        protChangeTimeUnit = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (t to min)',
                                              labelToChange='t', newUnitsCategory=0, newUnitsCategoryTime=1)
        protChangeTimeUnit.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protChangeTimeUnit)
        self.assertIsNotNone(protChangeTimeUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeTimeUnit)

        # Filter time variable
        print "Filter time"
        protFilterTime = self.newProtocol(ProtPKPDFilterMeasurements,
                                          objLabel='pkpd - filter measurements (t<24h)',
                                          filterType=1, condition='$(t)<1440')
        protFilterTime.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protFilterTime)
        self.assertIsNotNone(protFilterTime.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('protFilterTime', protFilterTime)

        # Fit a model to the first dose
        print "Fitting monocompartmental model with 1st order to 1st dose..."
        protEV1MonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                  objLabel='pkpd - ev1 monocompartment',
                                                  bounds='(25, 60.0); (0.0, 0.01); (0.08, 0.2); (40.0, 70)')
        protEV1MonoCompartment.inputExperiment.set(protFilterTime.outputExperiment)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment.outputExperiment.fnPKPD)
        tlag = float(experiment.samples['Individual'].descriptors['Oral_tlag'])
        Ka = float(experiment.samples['Individual'].descriptors['Oral_Ka'])
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        self.assertTrue(tlag>35 and tlag<50)
        self.assertTrue(Ka>0.002 and Ka<0.003)
        self.assertTrue(Cl>0.13 and Cl<0.14)
        self.assertTrue(V>30 and V<60)
        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.98)
        self.assertTrue(fitting.sampleFits[0].AIC<-60)

        # Fit a model to the all doses
        print "Fitting monocompartmental model with 1st order to all doses..."
        protEV1MonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                  objLabel='pkpd - ev1 monocompartment',
                                                  bounds='(25, 60.0); (0.001, 0.005); (0.1, 0.15); (20.0, 70)',
                                                  globalSearch=False)
        protEV1MonoCompartment.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment.outputExperiment.fnPKPD)
        tlag = float(experiment.samples['Individual'].descriptors['Oral_tlag'])
        Ka = float(experiment.samples['Individual'].descriptors['Oral_Ka'])
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        self.assertTrue(tlag>35 and tlag<50) # Gabrielsson p 531 tlag=0.697h=41min ------- Mine: 40.66
        self.assertTrue(Ka>0.002 and Ka<0.004) # Gabrielsson p 531 k01=0.1964 1/h=0.0033 1/min ------- Mine: 0.0032
        self.assertTrue(Cl>0.13 and Cl<0.14)
        self.assertTrue(V>30 and V<70) # Gabrielsson p 531 Volume/F=94.946 ------- Mine: 63.73
        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.98)
        self.assertTrue(fitting.sampleFits[0].AIC<-90)

if __name__ == "__main__":
    unittest.main()
