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
import copy

class TestImportExportCSVExcelWorkflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Dissolution')
        cls.exptFn = cls.dataset.getFile('excel')

    def testDissolutionWorkflow(self):
        print "Import Excel"
        protImport = self.newProtocol(ProtPKPDImportFromTable,
                                      objLabel='pkpd - import excel',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment_Sheet2.fnPKPD, "There was a problem with the import")
        self.assertIsNotNone(protImport.outputExperiment_experiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)
        experiment = PKPDExperiment()
        experiment.load(protImport.outputExperiment_Sheet2.fnPKPD)
        self.assertTrue(len(experiment.samples)==12)
        experiment = PKPDExperiment()
        experiment.load(protImport.outputExperiment_experiment.fnPKPD)
        self.assertTrue(len(experiment.samples)==12)

        print "Export by rows"
        protExport1 = self.newProtocol(ProtPKPDExportToCSV,
                                      objLabel='pkpd - export by rows',
                                      format=0)
        protExport1.inputExperiment.set(protImport.outputExperiment_Sheet2)
        self.launchProtocol(protExport1)


        print "Export as table"
        protExport2 = self.newProtocol(ProtPKPDExportToCSV,
                                       objLabel='pkpd - export by table',
                                       format=1)
        protExport2.inputExperiment.set(protImport.outputExperiment_Sheet2)
        self.launchProtocol(protExport2)


        print "Import from rows"
        protImport1 = self.newProtocol(ProtPKPDImportFromCSV,
                                       objLabel='pkpd - import from rows',
                                       inputFile=protExport1._getPath("experiment.csv"),
                                       variables="t ; h ; numeric ; time ; ;; C ; none ; numeric ; measurement ; dissolution")
        self.launchProtocol(protImport1)
        self.assertIsNotNone(protImport1.outputExperiment.fnPKPD, "There was a problem with the import")
        experiment = PKPDExperiment()
        experiment.load(protImport1.outputExperiment.fnPKPD)
        self.assertTrue(len(experiment.samples)==12)


        print "Import from table"
        protImport2 = self.newProtocol(ProtPKPDImportFromTable,
                                       objLabel='pkpd - import from table',
                                       inputFile=protExport2._getPath("experiment.csv"))
        self.launchProtocol(protImport2)
        self.assertIsNotNone(protImport2.outputExperiment.fnPKPD, "There was a problem with the import")
        experiment = PKPDExperiment()
        experiment.load(protImport2.outputExperiment.fnPKPD)
        self.assertTrue(len(experiment.samples) == 12)

if __name__ == "__main__":
    unittest.main()
