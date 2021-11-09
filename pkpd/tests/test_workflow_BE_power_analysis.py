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
from .test_workflow import TestWorkflow

class TestBEPowerAnalysisWorkflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)

    def testDissolutionWorkflow(self):
        print("ABE")
        prot = self.newProtocol(ProtPKPDBEPowerAnalysis,
                                objLabel='pkpd - BE ABE',
                                method=0)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath('plot.png')))

        print("scABE")
        prot = self.newProtocol(ProtPKPDBEPowerAnalysis,
                                objLabel='pkpd - BE scABE',
                                method=1)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath('plot.png')))

        print("NTID")
        prot = self.newProtocol(ProtPKPDBEPowerAnalysis,
                                objLabel='pkpd - BE NTID',
                                method=2)
        self.launchProtocol(prot)
        self.assertTrue(os.path.exists(prot._getPath('plot.png')))

if __name__ == "__main__":
    unittest.main()
