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


class TestGabrielssonPK10Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK10')
        cls.exptFnIV = cls.dataset.getFile('experimentIV')
        cls.exptFnPO = cls.dataset.getFile('experimentPO')

    def testGabrielssonPK10Workflow(self):
        # Import an experiment

        print "Import Experiment IV"
        protImportIV = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFnIV)
        self.launchProtocol(protImportIV)
        self.assertIsNotNone(protImportIV.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImportIV)

        # Fit a two-compartmentx model with intravenous absorption to a set of measurements
        print "Fitting a two-compartment model (intravenous)..."
        protPKPDIVTwoCompartments = self.newProtocol(ProtPKPDTwoCompartments,
                                                     objLabel='pkpd - iv two-compartments',
                                                     globalSearch=False,
                                                     bounds='(0.0, 2.0); (0.0, 100.0); (0.0, 3.0); (0.0, 100.0)')
        protPKPDIVTwoCompartments.inputExperiment.set(protImportIV.outputExperiment)
        self.launchProtocol(protPKPDIVTwoCompartments)
        self.assertIsNotNone(protPKPDIVTwoCompartments.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protPKPDIVTwoCompartments.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDIVTwoCompartments', protPKPDIVTwoCompartments)
        experiment = PKPDExperiment()
        experiment.load(protPKPDIVTwoCompartments.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        self.assertTrue(Cl>0.98 and Cl<1.06)
        self.assertTrue(Clp>2.38 and Clp<2.46)
        self.assertTrue(V>50 and V<55)
        self.assertTrue(Vp>42 and Vp<47)
        fitting = PKPDFitting()
        fitting.load(protPKPDIVTwoCompartments.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.98)
        self.assertTrue(fitting.sampleFits[0].AIC<-35)

        print "Import Experiment PO"
        protImportPO = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFnPO)
        self.launchProtocol(protImportPO)
        self.assertIsNotNone(protImportPO.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImportPO)

        # Fit a two-compartment model with oral absorption to a set of measurements
        print "Fitting a two-compartment model (oral)..."
        protPKPDPOTwoCompartments = self.newProtocol(ProtPKPDTwoCompartments,
                                                     objLabel='pkpd - ev1 two-compartments',
                                                     globalSearch=False,
                                                     bounds='(10.0, 20.0); (0.25, 0.45); (0.02, 0.06); (0.9, 1.15); (40.0, 70.0); (1.0, 3.0); (40.0, 70.0)')
        protPKPDPOTwoCompartments.inputExperiment.set(protImportPO.outputExperiment)
        self.launchProtocol(protPKPDPOTwoCompartments)
        self.assertIsNotNone(protPKPDPOTwoCompartments.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protPKPDPOTwoCompartments.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDPOTwoCompartments', protPKPDPOTwoCompartments)
        experiment = PKPDExperiment()
        experiment.load(protPKPDPOTwoCompartments.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        Ka = float(experiment.samples['Individual'].descriptors['Oral_Ka'])
        F = float(experiment.samples['Individual'].descriptors['Oral_bioavailability'])
        tlag = float(experiment.samples['Individual'].descriptors['Oral_tlag'])
        self.assertTrue(Cl>0.95 and Cl<1.06)
        self.assertTrue(Clp>1 and Clp<2.4)
        self.assertTrue(V>45)
        self.assertTrue(Vp>40 and Vp<70)
        self.assertTrue(Ka>0.018 and Ka<0.05)
        self.assertTrue(F>0.27 and F<0.38)
        self.assertTrue(tlag>9 and tlag<20)
        fitting = PKPDFitting()
        fitting.load(protPKPDPOTwoCompartments.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.6)

        # Fit a two-compartmentx model with oral absorption to a set of measurements
        print "Fitting IV and PO simultaneously ..."
        protTwoVias = self.newProtocol(ProtPKPDODETwoVias,
                                                     objLabel='pkpd - ode two vias',
                                                     globalSearch=False)
        protTwoVias.prot1ptr.set(protPKPDIVTwoCompartments)
        protTwoVias.prot2ptr.set(protPKPDPOTwoCompartments)
        self.launchProtocol(protTwoVias)
        self.assertIsNotNone(protTwoVias.outputExperiment1.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protTwoVias.outputFitting1.fnFitting, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protTwoVias.outputExperiment2.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protTwoVias.outputFitting2.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protTwoVias', protTwoVias)

        experiment = PKPDExperiment()
        experiment.load(protTwoVias.outputExperiment2.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        Ka = float(experiment.samples['Individual'].descriptors['Oral_Ka'])
        F = float(experiment.samples['Individual'].descriptors['Oral_bioavailability'])
        tlag = float(experiment.samples['Individual'].descriptors['Oral_tlag'])
        self.assertTrue(Cl>1.00 and Cl<1.07)
        self.assertTrue(Clp>1 and Clp<2.4)
        self.assertTrue(V>50 and V<57) # Gabrielsson, p. 583: V=59.9
        self.assertTrue(Vp>52 and Vp<60)
        self.assertTrue(Ka>0.032 and Ka<0.043) # Gabrielsson, p. 583: Ka=0.047
        self.assertTrue(F>0.32 and F<0.38) # Gabrielsson, p. 583: B10=0.3187
        self.assertTrue(tlag>12 and tlag<18) # Gabrielsson, p. 583: tlag=14.82
        fitting = PKPDFitting()
        fitting.load(protTwoVias.outputFitting1.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.98)
        self.assertTrue(fitting.sampleFits[0].AICc<-40)

if __name__ == "__main__":
    unittest.main()
