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


class TestGabrielssonPK12Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK12')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK12Workflow(self):
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
        protChangeConcUnit = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (Cp to mg/L)',
                                              labelToChange='Cp', newUnitsCategory=4, newUnitsCategoryConc=1)
        protChangeConcUnit.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protChangeConcUnit)
        self.assertIsNotNone(protChangeConcUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeConcUnit)

        # Fit a two-compartment model with oral absorption to a set of measurements
        print "Fitting a two-compartment model ..."
        protPKPDPOTwoCompartments = self.newProtocol(ProtPKPDTwoCompartments,
                                                     objLabel='pkpd - ev1 two-compartments',
                                                     globalSearch=True,
                                                     bounds='(0.0, 5.0); (0.0, 0.15); (0.0, 0.2); (0.0, 0.02); (0.0, 1.0); (0.0, 0.02); (0.0, 1.0)')
        protPKPDPOTwoCompartments.inputExperiment.set(protChangeConcUnit.outputExperiment)
        self.launchProtocol(protPKPDPOTwoCompartments)
        self.assertIsNotNone(protPKPDPOTwoCompartments.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protPKPDPOTwoCompartments.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDPOTwoCompartments', protPKPDPOTwoCompartments)
        experiment = PKPDExperiment()
        experiment.load(protPKPDPOTwoCompartments.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        Ka = float(experiment.samples['Individual'].descriptors['Oral_Ka'])
        tlag = float(experiment.samples['Individual'].descriptors['Oral_tlag'])
        bioavailability = float(experiment.samples['Individual'].descriptors['Oral_bioavailability'])
        self.assertTrue(Cl>0.0125 and Cl<0.015) # Gabrielsson p. 598: Cl=0.0145
        self.assertTrue(Clp>0.008 and Clp<0.013) # Gabrielsson p. 598: Cld=0.0208
        self.assertTrue(V>0.19 and V<0.23) # Gabrielsson p. 598: Vc=0.120
        self.assertTrue(Vp>0.24 and Vp<0.32) # Gabrielsson p. 598: Vt=0.2759
        self.assertTrue(Ka>0.09 and Ka<0.17) # Gabrielsson p. 598: Ka=0.103
        self.assertTrue(tlag>0 and tlag<6) # Gabrielsson p. 598: tlag=4.67
        fitting = PKPDFitting()
        fitting.load(protPKPDPOTwoCompartments.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.98)


if __name__ == "__main__":
    unittest.main()
