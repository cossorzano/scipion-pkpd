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
from pkpd.objects import PKPDDataSet
from test_workflow import TestWorkflow


class TestGabrielssonPK17Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK17')
        cls.exptFn = cls.dataset.getFile('experiment')
        cls.exptFn2 = cls.dataset.getFile('experiment2')

    def testGabrielssonPK17Workflow(self):
        # Import an experiment (intravenous)

        print "Import Experiment (intravenous doses)"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Fit a mono-compartment model with intravenous absorption to a set of measurements
        print "Fitting a mono-compartment model ..."
        protPKPDPMonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                     objLabel='pkpd - iv mono-compartment',
                                                     globalSearch=False,fitType=0,
                                                     bounds='(0.0, 100.0); (0.0, 4000.0)')
        protPKPDPMonoCompartment.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPKPDPMonoCompartment)
        self.assertIsNotNone(protPKPDPMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the mono-compartmental model ")
        self.assertIsNotNone(protPKPDPMonoCompartment.outputFitting.fnFitting, "There was a problem with the mono-compartmental model ")
        self.validateFiles('protPKPDPMonoCompartment', protPKPDPMonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protPKPDPMonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        self.assertTrue(Cl>35.5 and Cl<36.5) # Gabrielsson p. 631: Cl=43.28
        self.assertTrue(V>1550 and V<1650) # Gabrielsson p. 631: V=1377.66
        fitting = PKPDFitting()
        fitting.load(protPKPDPMonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.8)

        # Fit a mono-compartment model with intravenous absorption to a set of measurements
        print "Fitting a mono-compartment model intrinsic ..."
        protPKPDPMonoCompartment = self.newProtocol(ProtPKPDMonoCompartmentClint,
                                                     objLabel='pkpd - iv mono-compartment intrinsic',
                                                     globalSearch=False,fitType=0,
                                                     bounds='(0.0, 200.0); (0.0, 2.0); (1000.0, 2000.0)')
        protPKPDPMonoCompartment.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPKPDPMonoCompartment)
        self.assertIsNotNone(protPKPDPMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the mono-compartmental model ")
        self.assertIsNotNone(protPKPDPMonoCompartment.outputFitting.fnFitting, "There was a problem with the mono-compartmental model ")
        self.validateFiles('protPKPDPMonoCompartment', protPKPDPMonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protPKPDPMonoCompartment.outputExperiment.fnPKPD)
        Vmax = float(experiment.samples['Individual'].descriptors['Vmax'])
        Km = float(experiment.samples['Individual'].descriptors['Km'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        self.assertTrue(Vmax>120 and Vmax<130) # Gabrielsson p. 632: Vm=107
        self.assertTrue(Km>1.2 and Km<1.25) # Gabrielsson p. 632: Km=0.566
        self.assertTrue(V>1500 and V<1600) # Gabrielsson p. 632: V=1454
        fitting = PKPDFitting()
        fitting.load(protPKPDPMonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.8)

        # This example is simulated data_test and it serves to verify that infusions are correctly simulated
        print "Import Experiment 2 (intravenous doses)"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment simulated',
                                      inputFile=self.exptFn2)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Fit a mono-compartment model with intravenous absorption to a set of measurements
        print "Fitting a mono-compartment model ..."
        protPKPDPMonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                     objLabel='pkpd - iv mono-compartment',
                                                     globalSearch=False,fitType=0,
                                                     bounds='(0.0, 1); (0.0, 100)')
        protPKPDPMonoCompartment.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPKPDPMonoCompartment)
        self.assertIsNotNone(protPKPDPMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the mono-compartmental model ")
        self.assertIsNotNone(protPKPDPMonoCompartment.outputFitting.fnFitting, "There was a problem with the mono-compartmental model ")
        self.validateFiles('protPKPDPMonoCompartment', protPKPDPMonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protPKPDPMonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        self.assertTrue(Cl>0.195 and Cl<0.205) # Domenech Vol 1, p 250: Cl=12 L/h=0.2 L/min
        self.assertTrue(V>69 and V<71) # Domenech Vol 1, p 250: V=70 L
        fitting = PKPDFitting()
        fitting.load(protPKPDPMonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

if __name__ == "__main__":
    unittest.main()
