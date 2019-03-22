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


class TestGabrielssonPK09Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK09')
        cls.exptFn = cls.dataset.getFile('experiment')
    
    def testGabrielssonPK09Workflow(self):
        # Import an experiment

        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Fit a two-compartmentx model with intravenous absorption to a set of measurements
        print "Fitting a two-compartmentx model (intravenous)..."
        protPKPDIVTwoCompartments = self.newProtocol(ProtPKPDTwoCompartments,
                                                     objLabel='pkpd - iv two-compartments',
                                                     globalSearch=False,
                                                     bounds='(0.0, 0.8); (0.0, 115.0); (0.0, 2); (0.0, 115.0)')
        protPKPDIVTwoCompartments.inputExperiment.set(protImport.outputExperiment)
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
        self.assertTrue(Cl>0.35 and Cl<0.40)
        self.assertTrue(Clp>0.98 and Clp<1.05)
        self.assertTrue(V>50 and V<60)
        self.assertTrue(Vp>50 and Vp<60)
        fitting = PKPDFitting()
        fitting.load(protPKPDIVTwoCompartments.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-65)

if __name__ == "__main__":
    unittest.main()
