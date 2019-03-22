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


class TestGabrielssonPK11Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK11')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK11Workflow(self):
        # Import an experiment
        print "Import Experiment"
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

        # Filter time variable
        print "Filter time"
        protFilterTime = self.newProtocol(ProtPKPDFilterMeasurements,
                                                  objLabel='pkpd - filter measurements t<=1440',
                                                  filterType=1, condition='$(t)<=1440')
        protFilterTime.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protFilterTime)
        self.assertIsNotNone(protFilterTime.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('protFilterTime', protFilterTime)

        # Fit a two-compartment model with oral absorption to a set of measurements
        print "Fitting a two-compartment model ..."
        protPKPDPOTwoCompartments = self.newProtocol(ProtPKPDTwoCompartments,
                                                     objLabel='pkpd - ev1 two-compartments',
                                                     globalSearch=False,
                                                     bounds='(13, 20.0); (0.0, 0.03); (0.05, 0.15); (0.5, 11); (0.01, 0.04); (9, 15)')
        protPKPDPOTwoCompartments.inputExperiment.set(protFilterTime.outputExperiment)
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
        self.assertTrue(Cl>0.085 and Cl<0.1)
        self.assertTrue(Clp>0.024 and Clp<0.028)
        self.assertTrue(V>8 and V<11)
        self.assertTrue(Vp>10 and Vp<14)
        self.assertTrue(Ka>0.008 and Ka<0.02)
        self.assertTrue(tlag>10 and tlag<17)
        fitting = PKPDFitting()
        fitting.load(protPKPDPOTwoCompartments.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Fit a two-compartment model with oral absorption to a set of measurements
        print "Fitting a two-compartment model ..."
        protPKPDPOTwoCompartments = self.newProtocol(ProtPKPDTwoCompartments,
                                                     objLabel='pkpd - ev1 two-compartments',
                                                     globalSearch=False,
                                                     bounds='(13, 20.0); (0.0, 0.03); (0.05, 0.15); (0.5, 11); (0.01, 0.04); (9, 18)')
        protPKPDPOTwoCompartments.inputExperiment.set(protChangeTimeUnit.outputExperiment)
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
        self.assertTrue(Cl>0.08 and Cl<0.09)
        self.assertTrue(Clp>0.01 and Clp<0.03)
        self.assertTrue(V>3.5 and V<4.5)
        self.assertTrue(Vp>10 and Vp<15)
        self.assertTrue(Ka>0.004 and Ka<0.03) # Gabrielsson p.590 K01=1.934 h^-1=0.032 min^-1
        self.assertTrue(tlag>10 and tlag<19) # Gabrielsson p.590, tlag=0.327 h=19.6 min
        fitting = PKPDFitting()
        fitting.load(protPKPDPOTwoCompartments.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.98)

        # Fit a two-compartment model with oral absorption to a set of measurements
        print "Fitting a PD model ..."
        protFitPD = self.newProtocol(ProtPKPDGenericFit,
                                     objLabel='pkpd - fit pd',
                                     modelType=2,
                                     bounds='(0.0, 20.0); (0.0, 100.0); (0.0, 10.0)')
        protFitPD.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protFitPD)
        self.assertIsNotNone(protFitPD.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protFitPD.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protFitPD', protFitPD)
        experiment = PKPDExperiment()
        experiment.load(protFitPD.outputExperiment.fnPKPD)
        e0 = float(experiment.samples['Individual'].descriptors['e0'])
        ec50 = float(experiment.samples['Individual'].descriptors['eC50'])
        emax = float(experiment.samples['Individual'].descriptors['emax'])
        self.assertTrue(e0>10 and e0<11)
        self.assertTrue(ec50>1.2 and ec50<1.4) # Gabrielsson p. 591: 0.78
        self.assertTrue(emax>88 and emax<92) # Gabrielsson p. 591: 97.2
        fitting = PKPDFitting()
        fitting.load(protFitPD.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.9)

if __name__ == "__main__":
    unittest.main()
