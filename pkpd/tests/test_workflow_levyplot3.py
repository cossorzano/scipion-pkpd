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

class TestLevyPlotWorkflow3(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Dissolution')
        cls.fnInVitro = cls.dataset.getFile('experiment12Test')
        cls.fnInVivo = cls.dataset.getFile('invivo12')

    def testDissolutionWorkflow(self):
        print("Import Experiment in vitro")
        protImportInVitro = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import in vitro',
                                      inputFile=self.fnInVitro)
        self.launchProtocol(protImportInVitro)
        self.assertIsNotNone(protImportInVitro.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImportInVitro)

        # Change the time unit to minute
        print("Change Units")
        protChangeTimeUnit = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (t to min)',
                                              labelToChange='t', newUnitsCategory=0, newUnitsCategoryTime=1)
        protChangeTimeUnit.inputExperiment.set(protImportInVitro.outputExperiment)
        self.launchProtocol(protChangeTimeUnit)
        self.assertIsNotNone(protChangeTimeUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeTimeUnit)

        print("Import Experiment in vivo")
        protImportInVivo = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import in vivo',
                                      inputFile=self.fnInVivo)
        self.launchProtocol(protImportInVivo)
        self.assertIsNotNone(protImportInVivo.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImportInVivo)

        # Fit a Weibull dissolution
        print("Fitting Weibull model ...")
        protWeibull = self.newProtocol(ProtPKPDDissolutionFit,
                                objLabel='pkpd - fit dissolution Weibull',
                                globalSearch=True, modelType=3)
        protWeibull.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protWeibull)
        self.assertIsNotNone(protWeibull.outputExperiment.fnPKPD, "There was a problem with the dissolution model")
        self.assertIsNotNone(protWeibull.outputFitting.fnFitting, "There was a problem with the dissolution model")
        self.validateFiles('ProtPKPDDissolutionFit', ProtPKPDDissolutionFit)

        # Fit Order 1
        print("Fitting EV1-monocompartment model ...")
        protModelInVivo = self.newProtocol(ProtPKPDMonoCompartment,
                                       objLabel='pkpd - fit monocompartment',
                                       bounds="(0.0, 8.0); (0.0, 0.2); (0.0, 10.0); (10.0, 50.0)"
                                       )
        protModelInVivo.inputExperiment.set(protImportInVivo.outputExperiment)
        self.launchProtocol(protModelInVivo)
        self.assertIsNotNone(protModelInVivo.outputExperiment.fnPKPD, "There was a problem with the PK model")
        self.assertIsNotNone(protModelInVivo.outputFitting.fnFitting, "There was a problem with the PK model")
        self.validateFiles('ProtPKPDMonoCompartment', ProtPKPDMonoCompartment)

        # Deconvolve the in vivo
        print("Deconvolving in vivo ...")
        protDeconv = self.newProtocol(ProtPKPDDeconvolve,
                                       objLabel='pkpd - deconvolution'
                                       )
        protDeconv.inputODE.set(protModelInVivo)
        self.launchProtocol(protDeconv)
        self.assertIsNotNone(protDeconv.outputExperiment.fnPKPD, "There was a problem with the deconvolution")
        self.validateFiles('ProtPKPDDeconvolve', ProtPKPDDeconvolve)

        # Levy plot
        print("Levy plot ...")
        protLevy = self.newProtocol(ProtPKPDDissolutionLevyPlot,
                                      objLabel='pkpd - levy plot'
                                      )
        protLevy.inputInVitro.set(protWeibull)
        protLevy.inputInVivo.set(protDeconv)
        self.launchProtocol(protLevy)
        self.assertIsNotNone(protLevy.outputExperiment.fnPKPD, "There was a problem with the Levy plot")
        self.validateFiles('ProtPKPDDissolutionLevyPlot', ProtPKPDDissolutionLevyPlot)

if __name__ == "__main__":
    unittest.main()
