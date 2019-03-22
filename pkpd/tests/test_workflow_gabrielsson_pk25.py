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


class TestGabrielssonPK25Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK25')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK25Workflow(self):
        # Import an experiment (intravenous)

        print "Import Experiment (intravenous doses)"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Change the time unit to minute
        print "Change Units"
        protChangeTimeUnit = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (t to min)',
                                              labelToChange='t', newUnitsCategory=0, newUnitsCategoryTime=1)
        protChangeTimeUnit.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protChangeTimeUnit)
        self.assertIsNotNone(protChangeTimeUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeTimeUnit)

        # # Fit a monocompartmental model to a set of measurements obtained by intravenous doses and urine
        print "Fitting a two-compartmental model (intravenous doses and urine) ..."
        protIVTwoCompartmentsUrine = self.newProtocol(ProtPKPDTwoCompartmentsUrine,
                                                      objLabel='pkpd - iv two-compartments urine',
                                                      globalSearch=False,
                                                      bounds='(2,4); (200,300); (0,1); (300,500); (0,0.1)')
        protIVTwoCompartmentsUrine.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protIVTwoCompartmentsUrine)
        self.assertIsNotNone(protIVTwoCompartmentsUrine.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protIVTwoCompartmentsUrine.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protIVTwoCompartmentsUrine', protIVTwoCompartmentsUrine)
        experiment = PKPDExperiment()
        experiment.load(protIVTwoCompartmentsUrine.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        fe = float(experiment.samples['Individual'].descriptors['fe'])
        self.assertTrue(Clp>0.65 and Clp<0.75)
        self.assertTrue(Cl>2.9 and Cl<3)
        self.assertTrue(V>240 and V<260)
        self.assertTrue(Vp>390 and Vp<410)
        self.assertTrue(fe>0.025 and fe<0.04)
        fitting = PKPDFitting()
        fitting.load(protIVTwoCompartmentsUrine.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-30)

if __name__ == "__main__":
    unittest.main()
