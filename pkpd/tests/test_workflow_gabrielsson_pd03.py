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
import pkpd.protocols
from test_workflow import TestWorkflow
from pkpd.objects import PKPDDataSet


class TestGabrielssonPD03Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PD03')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPD03Workflow(self):
        # Import an experiment (intravenous)

        print "Import Experiment"
        protImport = self.newProtocol(pkpd.protocols.ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Fit a PD model
        print "Fitting a saturated PD model ..."
        protPDfitting = self.newProtocol(pkpd.protocols.ProtPKPDGenericFit,
                                         objLabel='pkpd - pd saturated',
                                         modelType=2, fitType=0,
                                         bounds='(150.0, 200.0); (-100.0, 0.0); (50.0, 400.0)')
        protPDfitting.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPDfitting)
        self.assertIsNotNone(protPDfitting.outputExperiment.fnPKPD, "There was a problem with the mono-compartmental model ")
        self.assertIsNotNone(protPDfitting.outputFitting.fnFitting, "There was a problem with the mono-compartmental model ")
        self.validateFiles('protPDfitting', protPDfitting)
        experiment = pkpd.protocols.PKPDExperiment()
        experiment.load(protPDfitting.outputExperiment.fnPKPD)
        e0 = float(experiment.samples['Population measures'].descriptors['e0'])
        ec50 = float(experiment.samples['Population measures'].descriptors['eC50'])
        emax = float(experiment.samples['Population measures'].descriptors['emax'])
        self.assertTrue(e0>175 and e0<177) # Gabrielsson p 905: 176
        self.assertTrue(ec50>233 and ec50<234) # Gabrielsson p 905: 231
        self.assertTrue(emax>-57 and emax<-0.55) # Gabrielsson p 905: 56.1
        fitting = pkpd.protocols.PKPDFitting()
        fitting.load(protPDfitting.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.96)

        # Fit a PD model
        print "Fitting a sigmoid PD model ..."
        protPDfitting = self.newProtocol(pkpd.protocols.ProtPKPDGenericFit,
                                         objLabel='pkpd - pd sigmoid',
                                         modelType=3, fitType=0,
                                         bounds='(150.0, 200.0); (-100.0, 0.0); (50.0, 400.0); (0,4)')
        protPDfitting.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPDfitting)
        self.assertIsNotNone(protPDfitting.outputExperiment.fnPKPD, "There was a problem with the mono-compartmental model ")
        self.assertIsNotNone(protPDfitting.outputFitting.fnFitting, "There was a problem with the mono-compartmental model ")
        self.validateFiles('protPDfitting', protPDfitting)
        experiment = pkpd.protocols.PKPDExperiment()
        experiment.load(protPDfitting.outputExperiment.fnPKPD)
        e0 = float(experiment.samples['Population measures'].descriptors['e0'])
        ec50 = float(experiment.samples['Population measures'].descriptors['eC50'])
        emax = float(experiment.samples['Population measures'].descriptors['emax'])
        h = float(experiment.samples['Population measures'].descriptors['h'])
        self.assertTrue(e0>171 and e0<172) # Gabrielsson p 905: 172
        self.assertTrue(ec50>143 and ec50<144) # Gabrielsson p 905: 143
        self.assertTrue(emax>-36 and emax<-0.35) # Gabrielsson p 905: 35.7
        self.assertTrue(h>1.9 and h<2) # Gabrielsson p 905: 1.94
        fitting = pkpd.protocols.PKPDFitting()
        fitting.load(protPDfitting.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.97)

if __name__ == "__main__":
    unittest.main()
