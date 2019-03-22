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


class TestGabrielssonPK18Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK18')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK18Workflow(self):
        # Import an experiment (intravenous)

        print "Import Experiment (intravenous doses)"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Fit a mono-compartment model with intravenous absorption to a set of measurements
        print "Fitting a two-compartments model intrinsic ..."
        protPKPDTwoCompartment = self.newProtocol(ProtPKPDTwoCompartmentsClint,
                                                     objLabel='pkpd - iv two-compartments intrinsic',
                                                     globalSearch=False,
                                                     bounds='(0.0, 0.2); (0.0, 0.1); (0.0, 15.0); (0.0, 3.0); (0.0, 50.0)')
        protPKPDTwoCompartment.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPKPDTwoCompartment)
        self.assertIsNotNone(protPKPDTwoCompartment.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protPKPDTwoCompartment.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDTwoCompartment', protPKPDTwoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protPKPDTwoCompartment.outputExperiment.fnPKPD)
        Vmax = float(experiment.samples['Individual'].descriptors['Vmax'])
        Km = float(experiment.samples['Individual'].descriptors['Km'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        self.assertTrue(Vmax>0.08 and Vmax<0.09) # Gabrielsson p. 644: Vmax=0.0822
        self.assertTrue(Km>0.018 and Km<0.028) # Gabrielsson p. 644: Km=0.0148
        self.assertTrue(V>9 and V<10) # Gabrielsson p. 644: V=8.7
        self.assertTrue(Clp>1.2 and Clp<1.3) # Gabrielsson p. 644: Cld=1.3
        self.assertTrue(Vp>28 and Vp<30) # Gabrielsson p. 644: Vt=31.17
        fitting = PKPDFitting()
        fitting.load(protPKPDTwoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

if __name__ == "__main__":
    unittest.main()
