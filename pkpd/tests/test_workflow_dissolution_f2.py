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

import numpy as np

from pyworkflow.tests import *
from pkpd.protocols import *
from pkpd.objects import PKPDDataSet
from test_workflow import TestWorkflow

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
        self.assertTrue(os.path.exists(protImportTest.outputExperiment.fnPKPD.get()), "There was a problem with the import")
        self.validateFiles('protImportTest', protImportTest)

        protImportRef = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment ref',
                                      inputFile=self.exptRefFn)
        self.launchProtocol(protImportRef)
        self.assertTrue(os.path.exists(protImportRef.outputExperiment.fnPKPD.get()), "There was a problem with the import")
        self.validateFiles('protImportRef', protImportRef)

        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2')
        prot.inputRef.set(protImportRef.outputExperiment)
        prot.inputTest.set(protImportTest.outputExperiment)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")
        fh = open(prot._getPath("summary.txt"))
        for line in fh.readlines():
            if line.startswith("F2 mean+-std:"):
                tokens=line.split(":")
                tokens=tokens[1].split("+-")
                mean=float(tokens[0])
                self.assertTrue(mean>75 and mean<85)
            elif line.startswith("F1 mean+-std:"):
                tokens=line.split(":")
                tokens=tokens[1].split("+-")
                mean=float(tokens[0])
                self.assertTrue(mean>3 and mean<6)

        prot = self.newProtocol(ProtPKPDDissolutionF2,
                                objLabel='pkpd - dissol f1 & f2 resampling',
                                resampleT=0.25,
                                keepResample=True)
        prot.inputRef.set(protImportRef.outputExperiment)
        prot.inputTest.set(protImportTest.outputExperiment)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("summary.txt")), "There was a problem with the import")
        self.assertTrue(os.path.exists(prot._getPath("experiment.pkpd")), "There was a problem with the import")
        experiment = PKPDExperiment()
        experiment.load(prot._getPath("experiment.pkpd"))
        self.assertTrue(len(experiment.samples)==12)


        prot = self.newProtocol(ProtPKPDAverageSample,
                                objLabel='pkpd - sample average')
        prot.inputExperiment.set(protImportTest.outputExperiment)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("experiment.pkpd")), "There was a problem with the import")
        experiment = PKPDExperiment()
        experiment.load(prot._getPath("experiment.pkpd"))
        self.assertTrue(len(experiment.samples) == 1)
        self.assertTrue(np.abs(float(experiment.samples['AverageSample'].measurement_C[0])-0.3058)<1e-4)


        prot = self.newProtocol(ProtPKPDJoinSamples, prefix1="Test_", prefix2="Ref_",
                                objLabel='pkpd - join samples')
        prot.inputExperiment1.set(protImportTest.outputExperiment)
        prot.inputExperiment2.set(protImportRef.outputExperiment)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath("experiment.pkpd")), "There was a problem with the import")
        experiment = PKPDExperiment()
        experiment.load(prot._getPath("experiment.pkpd"))
        self.assertTrue(len(experiment.samples) == 24)

if __name__ == "__main__":
    unittest.main()
