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
from .test_workflow import TestWorkflow

class TestInhalation1Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Inhalation')

    def testDissolutionWorkflow(self):
        # Lung parameters
        print("Define lung parameters")
        protLung = self.newProtocol(ProtPKPDInhLungPhysiology,
                                      objLabel='pkpd - lung parameters')
        self.launchProtocol(protLung)
        self.assertIsNotNone(protLung.outputLungParameters.fnPhys, "There was a problem with the lung definition")

        # Substance parameters
        print("Substance (gold) ...")
        protGold= self.newProtocol(ProtPKPDInhSubstanceProperties,
                                   objLabel='pkpd - gold properties',
                                   name='gold',
                                   rho=19.3, MW=197,
                                   kdiss_alv=0, kdiss_br=0, kp_alv=0, kp_br=0, Cs_alv=1, Cs_br=1,
                                   Kpl_alv=1, Kpl_br=1, fu=1, R=0
                                   )
        self.launchProtocol(protGold)
        self.assertIsNotNone(protGold.outputSubstanceParameters.fnSubst, "There was a problem with the gold definition")

        # Deposition parameters
        print("Deposition ...")
        protDepo = self.newProtocol(ProtPKPDInhImportDepositionProperties,
                                    objLabel='pkpd - deposition',
                                    depositionFile=self.dataset.getFile('deposition1'))
        protDepo.substance.set(protGold.outputSubstanceParameters)
        protDepo.lungModel.set(protLung.outputLungParameters)
        self.launchProtocol(protDepo)
        self.assertIsNotNone(protDepo.outputDeposition.fnDeposition, "There was a problem with the deposition")

if __name__ == "__main__":
    unittest.main()
