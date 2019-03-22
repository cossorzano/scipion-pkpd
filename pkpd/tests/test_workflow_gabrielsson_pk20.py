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


class TestGabrielssonPK20Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK20')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK20Workflow(self):
        # Import an experiment (intravenous)

        print "Import Experiment (intravenous doses) Group"
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

        # Change the time unit to minute
        print "Change Units"
        protChangeConcUnit = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (Cp to mg/L)',
                                              labelToChange='Cp', newUnitsCategory=4, newUnitsCategoryConc=1)
        protChangeConcUnit.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protChangeConcUnit)
        self.assertIsNotNone(protChangeConcUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeConcUnit)

        # Fit a mono-compartment model with intravenous absorption to a set of measurements
        print "Fitting a mono-compartment model intrinsic ..."
        protPKPDPMonoCompartment = self.newProtocol(ProtPKPDMonoCompartmentClint,
                                                     objLabel='pkpd - iv mono-compartment intrinsic',
                                                     globalSearch=True,fitType=1,
                                                     bounds='(0,1); (0,1); (0,100)')
        protPKPDPMonoCompartment.inputExperiment.set(protChangeConcUnit.outputExperiment)
        self.launchProtocol(protPKPDPMonoCompartment)
        self.assertIsNotNone(protPKPDPMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the mono-compartmental model ")
        self.assertIsNotNone(protPKPDPMonoCompartment.outputFitting.fnFitting, "There was a problem with the mono-compartmental model ")
        self.validateFiles('protPKPDPMonoCompartment', protPKPDPMonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protPKPDPMonoCompartment.outputExperiment.fnPKPD)
        Vmax = float(experiment.samples['Individual100'].descriptors['Vmax'])
        Km = float(experiment.samples['Individual100'].descriptors['Km'])
        V = float(experiment.samples['Individual100'].descriptors['V'])
        self.assertTrue(Vmax>0.55 and Vmax<0.60) # Gabrielsson p. 671: Vmax=37775 ug/h=0.629 mg/min
        self.assertTrue(Km>0.68 and Km<0.69) # Gabrielsson p. 671: Km2=812 ug/L=0.812 mg/L
        self.assertTrue(V>49 and V<50) # Gabrielsson p. 671: V=48.56 L

        Vmax = float(experiment.samples['Individual25'].descriptors['Vmax'])
        Km = float(experiment.samples['Individual25'].descriptors['Km'])
        V = float(experiment.samples['Individual25'].descriptors['V'])
        self.assertTrue(Vmax>0.85 and Vmax<0.90) # Gabrielsson p. 671: Vmax=37775 ug/h=0.629 mg/min
        self.assertTrue(Km>0.45 and Km<0.50) # Gabrielsson p. 671: Km1=282 ug/L=0.282 mg/L
        self.assertTrue(V>44 and V<46) # Gabrielsson p. 671: V=48.56 L

        fitting = PKPDFitting()
        fitting.load(protPKPDPMonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[1].R2>0.99)

if __name__ == "__main__":
    unittest.main()
