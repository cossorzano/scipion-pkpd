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


class TestGabrielssonPK14Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK14')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK14Workflow(self):
        # Import an experiment
        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # # Change the concentration to mg/L
        print "Change Units"
        protChangeTimeUnit = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (t to min)',
                                              labelToChange='t', newUnitsCategory=0, newUnitsCategoryConc=1)
        protChangeTimeUnit.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protChangeTimeUnit)
        self.assertIsNotNone(protChangeTimeUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeTimeUnit)

        # Fit a two-compartment model with oral absorption to a set of measurements
        print "Fitting a two-compartment model ..."
        protPKPDPOTwoCompartments = self.newProtocol(ProtPKPDTwoCompartments,
                                                     objLabel='pkpd - iv two-compartments',
                                                     globalSearch=True,
                                                     bounds='(0.0, 30.0); (0.1, 0.3); (0.8, 1.2); (60, 110); (0.15, 0.25); (100, 120)')
        protPKPDPOTwoCompartments.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protPKPDPOTwoCompartments)
        self.assertIsNotNone(protPKPDPOTwoCompartments.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protPKPDPOTwoCompartments.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDPOTwoCompartments', protPKPDPOTwoCompartments)
        experiment = PKPDExperiment()
        experiment.load(protPKPDPOTwoCompartments.outputExperiment.fnPKPD)
        Ka = float(experiment.samples['Individual'].descriptors['Oral_Ka'])
        tlag = float(experiment.samples['Individual'].descriptors['Oral_tlag'])
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        self.assertTrue(Ka>0.185 and Ka<0.29) # Gabrielsson p. 613: k01=10 h^-1=0.166 min^-1
        self.assertTrue(tlag>5 and tlag<22) # Gabrielsson p. 613: Tlag=0.078 h=4.7 min
        self.assertTrue(Clp>0.165 and Clp<0.22) # Gabrielsson p. 613: k12=0.13
        self.assertTrue(Cl>0.85 and Cl<1.1) # Gabrielsson p. 613: k10=0.66
        self.assertTrue(V>82 and V<93) # Gabrielsson p. 613: Vc=83
        self.assertTrue(Vp>95 and Vp<121)
        fitting = PKPDFitting()
        fitting.load(protPKPDPOTwoCompartments.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)


if __name__ == "__main__":
    unittest.main()
