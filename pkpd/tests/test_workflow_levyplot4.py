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
import os

from pyworkflow.tests import *
from pkpd.protocols import *
from pkpd.objects import PKPDDataSet
from test_workflow import TestWorkflow

class TestLevyPlotWorkflow4(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Dissolution')
        cls.fnInVivo = cls.dataset.getFile('invivo1')

    def testDissolutionWorkflow(self):
        print("Import Experiment in vivo")
        protImportInVivo = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import in vivo',
                                      inputFile=self.fnInVivo)
        self.launchProtocol(protImportInVivo)
        self.assertIsNotNone(protImportInVivo.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImportInVivo)

        # Fit bicompartment
        print("Fitting EV1-twocompartments model ...")
        protModelInVivo2 = self.newProtocol(ProtPKPDTwoCompartments,
                                       objLabel='pkpd - fit two compartments',
                                       bounds="(0.0, 8.0); (0.0, 0.2); (0.0, 10.0); (0.1, 50.0);  (0.0, 10.0); (0.1, 50.0)"
                                       )
        protModelInVivo2.inputExperiment.set(protImportInVivo.outputExperiment)
        self.launchProtocol(protModelInVivo2)
        self.assertIsNotNone(protModelInVivo2.outputExperiment.fnPKPD, "There was a problem with the PK model")
        self.assertIsNotNone(protModelInVivo2.outputFitting.fnFitting, "There was a problem with the PK model")
        self.validateFiles('ProtPKPDTwoCompartments', ProtPKPDTwoCompartments)

        # Deconvolve the in vivo
        print("Deconvolving in vivo Loo Riegelman ...")
        protDeconvLR = self.newProtocol(ProtPKPDDeconvolutionLooRiegelman,
                                        objLabel='pkpd - deconvolution Loo Riegelman'
                                       )
        protDeconvLR.inputExperiment.set(protModelInVivo2.outputExperiment)
        self.launchProtocol(protDeconvLR)
        self.assertIsNotNone(protDeconvLR.outputExperiment.fnPKPD, "There was a problem with the deconvolution Loo")
        self.validateFiles('ProtPKPDDeconvolutionLooRiegelman', ProtPKPDDeconvolutionLooRiegelman)

        # Change the time unit to minute
        print("Change Units")
        protChangeAunit = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (A to ug)',
                                              labelToChange='A', newUnitsCategory=2, newUnitsCategoryWeight=3)
        protChangeAunit.inputExperiment.set(protDeconvLR.outputExperiment)
        self.launchProtocol(protChangeAunit)
        self.assertIsNotNone(protChangeAunit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeAunit', protChangeAunit)

        # Inverse Loo-Riegelman
        print("Inverse Loo Riegelman ...")
        protInverseLR = self.newProtocol(ProtPKPDDeconvolutionLooRiegelmanInverse,
                                         objLabel='pkpd - inverse Loo Riegelman',
                                         k10=0.0185,
                                         k12=0.57,
                                         k21=2.54,
                                         V=16.7,
                                         Vp=3.76
                                        )
        protInverseLR.inputExperiment.set(protChangeAunit.outputExperiment)
        self.launchProtocol(protInverseLR)
        self.assertIsNotNone(protInverseLR.outputExperiment.fnPKPD, "There was a problem with the deconvolution Loo")
        experiment = PKPDExperiment()
        experiment.load(protInverseLR.outputExperiment.fnPKPD)
        C = np.asarray(experiment.samples['Individual1'].getValues('C'),dtype=np.float64)
        Cmax = np.max(C)
        self.assertTrue(Cmax>1.8 and Cmax<2.3)

if __name__ == "__main__":
    unittest.main()
