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


class TestGabrielssonPK43Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK43')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK43Workflow(self):
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
        print "Fitting a two-compartmental model ..."
        protMonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                      objLabel='pkpd - mono-compartment',
                                                      globalSearch=False,
                                                      bounds='(0.0, 0.2); (0.0, 0.05); (60.0, 200.0); (0.0, 1.0); (0.0, 0.05); (0.0, 30.0)')
        protMonoCompartment.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protMonoCompartment)
        self.assertIsNotNone(protMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protMonoCompartment.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protMonoCompartment', protMonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protMonoCompartment.outputExperiment.fnPKPD)
        Ka1 = float(experiment.samples['Individual'].descriptors['Oral_Ka1'])
        Ka2 = float(experiment.samples['Individual'].descriptors['Oral_Ka2'])
        tlag12 = float(experiment.samples['Individual'].descriptors['Oral_tlag12'])
        F1 = float(experiment.samples['Individual'].descriptors['Oral_F1'])
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        self.assertTrue(Ka1>0.11 and Ka1<0.12) # Gabrielsson p 826: Ka1=7.62 h^-1=0.127 min^-1
        self.assertTrue(Ka2>0.015 and Ka2<0.025) # Gabrielsson p 826: Ka2=1.07 h^-1=0.01783 min^-1
        self.assertTrue(tlag12>135 and tlag12<145) # Gabrielsson p 826: tlag=2.29 h= 137 min
        self.assertTrue(F1>0.5 and F1<0.55) # Gabrielsson p 826: frct=0.514
        self.assertTrue(Cl>0.025 and Cl<0.035) # Gabrielsson p 826: 0.0899 L/h^-1=.00149 L/min^-1
        self.assertTrue(V>20 and V<22) # Gabrielsson p 826: V=20.61
        fitting = PKPDFitting()
        fitting.load(protMonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-80)

if __name__ == "__main__":
    unittest.main()
