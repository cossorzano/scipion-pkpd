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


class TestGabrielssonPK15Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK15')
        cls.exptFn = cls.dataset.getFile('experiment')
        cls.exptPDFn = cls.dataset.getFile('experimentPD')

    def testGabrielssonPK15Workflow(self):
        # Import an experiment
        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Fit a mono-compartment model with fractional oral absorption to a set of measurements
        print "Fitting a mono-compartment model ..."
        protPKPDPOMonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                     objLabel='pkpd - iv-ev1 mono-compartment',
                                                     globalSearch=True, fitType=2,
                                                     bounds='(0.0, 0.01); (0.9, 1.0); (0.2, 0.8); (3, 6)')
        protPKPDPOMonoCompartment.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPKPDPOMonoCompartment)
        self.assertIsNotNone(protPKPDPOMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the mono-compartmental model ")
        self.assertIsNotNone(protPKPDPOMonoCompartment.outputFitting.fnFitting, "There was a problem with the mono-compartmental model ")
        self.validateFiles('protPKPDPOMonoCompartments', protPKPDPOMonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protPKPDPOMonoCompartment.outputExperiment.fnPKPD)
        Ka = float(experiment.samples['Individual'].descriptors['IVEV1_Ka'])
        F = float(experiment.samples['Individual'].descriptors['IVEV1_F'])
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        self.assertTrue(Ka>0.005 and Ka<0.006)
        self.assertTrue(F>0.92 and F<0.97)
        self.assertTrue(Cl>0.45 and Cl<0.81)
        self.assertTrue(V>3.2 and V<5.7)
        fitting = PKPDFitting()
        fitting.load(protPKPDPOMonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.85)

        # Import an experiment
        print "Import Experiment PD"
        protImportPD = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment PD',
                                      inputFile=self.exptPDFn)
        self.launchProtocol(protImportPD)
        self.assertIsNotNone(protImportPD.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImportPD', protImportPD)

        # Fit a PD model
        print "Fitting a PD model ..."
        protPDfitting = self.newProtocol(ProtPKPDGenericFit,
                                         objLabel='pkpd - pd fitting', predicted='S',
                                         modelType=3, fitType=0,
                                         bounds='(60.0, 100.0); (-100.0, 0.0); (0.0, 1.0); (0.0, 10.0)')
        protPDfitting.inputExperiment.set(protImportPD.outputExperiment)
        self.launchProtocol(protPDfitting)
        self.assertIsNotNone(protPDfitting.outputExperiment.fnPKPD, "There was a problem with the mono-compartmental model ")
        self.assertIsNotNone(protPDfitting.outputFitting.fnFitting, "There was a problem with the mono-compartmental model ")
        self.validateFiles('protPDfitting', protPDfitting)
        experiment = PKPDExperiment()
        experiment.load(protPDfitting.outputExperiment.fnPKPD)
        e0 = float(experiment.samples['Population measures'].descriptors['e0'])
        ec50 = float(experiment.samples['Population measures'].descriptors['eC50'])
        emax = float(experiment.samples['Population measures'].descriptors['emax'])
        h = float(experiment.samples['Population measures'].descriptors['h'])
        self.assertTrue(e0>82 and e0<84) # Gabrielsson p 621: 82.25
        self.assertTrue(ec50>0.35 and ec50<0.39) # Gabrielsson p 621: 0.3488
        self.assertTrue(emax>-91 and emax<-0.88)
        self.assertTrue(h>5.3 and h<5.4) # Gabrielsson p 621: 6.69
        fitting = PKPDFitting()
        fitting.load(protPDfitting.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

if __name__ == "__main__":
    unittest.main()
