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

import os

from pyworkflow.tests import *
from pkpd.protocols import *
from pkpd.objects import PKPDDataSet
from test_workflow import TestWorkflow

class TestDissolutionF2Shah1998Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Dissolution')
        cls.expFn = cls.dataset.getFile('excelShah1998')
        cls.expAvgFn = cls.dataset.getFile('excelShah1998avg')

    def testDissolutionWorkflow(self):
        print "Import Experiment"
        protImport = self.newProtocol(ProtPKPDImportFromTable,
                                      objLabel='pkpd - import experiments',
                                      inputFile=self.expFn)
        self.launchProtocol(protImport)
        self.assertTrue(os.path.exists(protImport.outputExperiment_Reference.fnPKPD.get()), "There was a problem with the import")
        self.assertTrue(os.path.exists(protImport.outputExperiment_Test5.fnPKPD.get()), "There was a problem with the import")
        self.validateFiles('protImportTest', protImport)

        print "F2 Test1"
        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2 Test1')
        prot.inputRef.set(protImport.outputExperiment_Reference)
        prot.inputTest.set(protImport.outputExperiment_Test1)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")

        print "F2 Test2"
        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2 Test2')
        prot.inputRef.set(protImport.outputExperiment_Reference)
        prot.inputTest.set(protImport.outputExperiment_Test2)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")

        print "F2 Test3"
        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2 Test3')
        prot.inputRef.set(protImport.outputExperiment_Reference)
        prot.inputTest.set(protImport.outputExperiment_Test3)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")

        print "F2 Test4"
        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2 Test4')
        prot.inputRef.set(protImport.outputExperiment_Reference)
        prot.inputTest.set(protImport.outputExperiment_Test4)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")

        print "F2 Test5"
        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2 Test5')
        prot.inputRef.set(protImport.outputExperiment_Reference)
        prot.inputTest.set(protImport.outputExperiment_Test5)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")


        print "Import Experiment"
        protImport = self.newProtocol(ProtPKPDImportFromTable,
                                      objLabel='pkpd - import experiments avg',
                                      inputFile=self.expAvgFn)
        self.launchProtocol(protImport)
        self.assertTrue(os.path.exists(protImport.outputExperiment_Reference.fnPKPD.get()),
                        "There was a problem with the import")
        self.assertTrue(os.path.exists(protImport.outputExperiment_Test5.fnPKPD.get()), "There was a problem with the import")
        self.validateFiles('protImportTest', protImport)

        print "F2 Test1"
        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2 avg Test1')
        prot.inputRef.set(protImport.outputExperiment_Reference)
        prot.inputTest.set(protImport.outputExperiment_Test1)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")

        print "F2 Test2"
        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2 avg Test2')
        prot.inputRef.set(protImport.outputExperiment_Reference)
        prot.inputTest.set(protImport.outputExperiment_Test2)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")

        print "F2 Test3"
        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2 avg Test3')
        prot.inputRef.set(protImport.outputExperiment_Reference)
        prot.inputTest.set(protImport.outputExperiment_Test3)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")

        print "F2 Test4"
        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2 avg Test4')
        prot.inputRef.set(protImport.outputExperiment_Reference)
        prot.inputTest.set(protImport.outputExperiment_Test4)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")

        print "F2 Test5"
        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2 avg Test5')
        prot.inputRef.set(protImport.outputExperiment_Reference)
        prot.inputTest.set(protImport.outputExperiment_Test5)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")

if __name__ == "__main__":
    unittest.main()
