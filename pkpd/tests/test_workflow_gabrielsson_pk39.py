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


class TestGabrielssonPK39Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK39')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK39Workflow(self):
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

        # Fit a monocompartmental model to a set of measurements obtained by intravenous doses and urine
        print "Fitting a two-compartmental model (intravenous doses) ..."
        protIVTwoCompartments = self.newProtocol(ProtPKPDTwoCompartments,
                                                      objLabel='pkpd - iv two-compartments',
                                                      bounds='(0,0.03); (0,1); (0,0.03); (1,3)')
        protIVTwoCompartments.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protIVTwoCompartments)
        self.assertIsNotNone(protIVTwoCompartments.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protIVTwoCompartments.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protIVTwoCompartmentsUrine', protIVTwoCompartments)
        experiment = PKPDExperiment()
        experiment.load(protIVTwoCompartments.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        self.assertTrue(Clp>0.011 and Clp<0.014) # Gabrielsson p 801: 0.8985 h^-1=.014975 min^-1
        self.assertTrue(Cl>0.006 and Cl<0.008) # Gabrielsson p 801: 0.417 h^-1=.006950 min^-1
        self.assertTrue(V>0.35 and V<0.37) # Gabrielsson p 801: Vc=0.322
        self.assertTrue(Vp>2.1 and Vp<2.3) # Gabrielsson p 802: Vt=2.14
        fitting = PKPDFitting()
        fitting.load(protIVTwoCompartments.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-50)

if __name__ == "__main__":
    unittest.main()
