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

class TestDissolutionF2Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Dissolution')
        cls.exptTestFn = cls.dataset.getFile('experiment12Test')
        cls.exptRefFn = cls.dataset.getFile('experiment12Ref')

    def testDissolutionWorkflow(self):
        print "Import Experiment"
        protImportTest = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment test',
                                      inputFile=self.exptTestFn)
        self.launchProtocol(protImportTest)
        self.assertIsNotNone(protImportTest.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImportTest', protImportTest)

        protImportRef = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment ref',
                                      inputFile=self.exptRefFn)
        self.launchProtocol(protImportRef)
        self.assertIsNotNone(protImportRef.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImportRef', protImportRef)

        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2')
        prot.inputRef.set(protImportRef.outputExperiment)
        prot.inputTest.set(protImportTest.outputExperiment)
        self.launchProtocol(prot)
        self.assertIsNotNone(prot._getPath("summary.txt"), "There was a problem with the import")

if __name__ == "__main__":
    unittest.main()
