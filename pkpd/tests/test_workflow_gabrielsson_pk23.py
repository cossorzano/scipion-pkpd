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


class TestGabrielssonPK23Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK23')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK23Workflow(self):
        # Import an experiment (intravenous)

        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # # Change the time unit to minute
        # print "Change Units"
        # protChangeTimeUnit = self.newProtocol(ProtPKPDChangeUnits,
        #                                       objLabel='pkpd - change units (t to min)',
        #                                       labelToChange='t', newUnitsCategory=0, newUnitsCategoryTime=1)
        # protChangeTimeUnit.inputExperiment.set(protImport.outputExperiment)
        # self.launchProtocol(protChangeTimeUnit)
        # self.assertIsNotNone(protChangeTimeUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        # self.validateFiles('protChangeUnits', protChangeTimeUnit)
        #
        # # Change the time unit to minute
        # print "Change Units"
        # protChangeConcUnit = self.newProtocol(ProtPKPDChangeUnits,
        #                                       objLabel='pkpd - change units (Cp to mg/L)',
        #                                       labelToChange='Cp', newUnitsCategory=4, newUnitsCategoryConc=1)
        # protChangeConcUnit.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        # self.launchProtocol(protChangeConcUnit)
        # self.assertIsNotNone(protChangeConcUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        # self.validateFiles('protChangeUnits', protChangeConcUnit)
        #
        # # Fit a two-compartment model with intravenous absorption to a set of measurements
        # print "Fitting a two-compartment model autoinduction ..."
        # prot = self.newProtocol(ProtPKPDTwoCompartmentsAutoinduction,
        #                                              objLabel='pkpd - iv two-compartment autoinduction',
        #                                              globalSearch=False,fitType=1,
        #                                              bounds='(0, 10.0); (0.0, 0.1); (0.0, 0.002); (135.0, 175.0); (0.75, 1.15); (30.0, 70.0)')
        # prot.inputExperiment.set(protChangeConcUnit.outputExperiment)
        # self.launchProtocol(prot)
        # self.assertIsNotNone(prot.outputExperiment.fnPKPD, "There was a problem with the mono-compartmental model ")
        # self.assertIsNotNone(prot.outputFitting.fnFitting, "There was a problem with the mono-compartmental model ")
        # self.validateFiles('ProtPKPDTwoCompartmentsAutoinduction', ProtPKPDTwoCompartmentsAutoinduction)
        # experiment = PKPDExperiment()
        # experiment.load(prot.outputExperiment.fnPKPD)
        # E0 = float(experiment.samples['Individual'].descriptors['E0'])
        # a = float(experiment.samples['Individual'].descriptors['a'])
        # kout = float(experiment.samples['Individual'].descriptors['kout'])
        # V = float(experiment.samples['Individual'].descriptors['V'])
        # Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        # Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        # self.assertTrue(E0>3.8 and E0<3.9) # Gabrielsson p. 688, E0=138
        # self.assertTrue(a>0.05 and a<0.06) # Gabrielsson p. 688, a=0.041 L/h=0.00068 L/min
        # self.assertTrue(kout>0.0008 and kout<0.001) # Gabrielsson p. 688, kout=0.023
        # self.assertTrue(V>155 and V<160) # Gabrielsson p. 688: V=146L
        # self.assertTrue(Clp>0.9 and Clp<1) # Gabrielsson p. 688, Cld=120 L/h=2 L/min
        # self.assertTrue(Vp>48 and Vp<52) # Gabrielsson p. 688, Vt=58.4 L
        #
        # fitting = PKPDFitting()
        # fitting.load(prot.outputFitting.fnFitting)
        # self.assertTrue(fitting.sampleFits[0].R2>0.99)

if __name__ == "__main__":
    unittest.main()
