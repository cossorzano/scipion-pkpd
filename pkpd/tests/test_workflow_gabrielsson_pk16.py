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


class TestGabrielssonPK16Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK16')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPK16Workflow(self):
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

        # Fit a two-compartment model with oral absorption to a set of measurements
        print "Fitting a two-compartment model ..."
        protPKPDPOTwoCompartments = self.newProtocol(ProtPKPDTwoCompartments,
                                                     objLabel='pkpd - iv two-compartments',
                                                     globalSearch=False,
                                                     bounds='(0.002, 0.01); (1, 2); (0.0001, 0.001); (0, 0.4)')
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
        self.assertTrue(Clp>0.00045 and Clp<0.00055) # Gabrielsson p. 626: Cld=0.030407 h^-1=0.0005067 min^-1
        self.assertTrue(Cl>0.005 and Cl<0.007) # Gabrielsson p. 626: (ClR+Clm)=(0.3149+0.0540) h^-1=.006148 min^-1
        self.assertTrue(V>1.5 and V<1.7) # Gabrielsson p. 626: Vc=1.61
        self.assertTrue(Vp>0.15 and Vp<0.17) # Gabrielsson p. 626: Vt=0.1647
        fitting = PKPDFitting()
        fitting.load(protPKPDPOTwoCompartments.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Fit a monocompartmental model to a set of measurements obtained by intravenous doses and urine
        print "Fitting a two-compartmental model (intravenous doses and urine) ..."
        protIVTwoCompartmentsUrine = self.newProtocol(ProtPKPDTwoCompartmentsUrine,
                                                      objLabel='pkpd - iv two-compartments urine',
                                                      globalSearch=False,
                                                      bounds='(0.002, 0.01); (1, 2); (0.0001, 0.001); (0, 0.4); (0,1)')
        protIVTwoCompartmentsUrine.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protIVTwoCompartmentsUrine)
        self.assertIsNotNone(protIVTwoCompartmentsUrine.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protIVTwoCompartmentsUrine.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protIVTwoCompartmentsUrine', protIVTwoCompartmentsUrine)
        experiment = PKPDExperiment()
        experiment.load(protIVTwoCompartmentsUrine.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        fe = float(experiment.samples['Individual'].descriptors['fe'])
        self.assertTrue(Clp>0.00045 and Clp<0.00055) # Gabrielsson p. 626: Cld=0.030407 h^-1=0.0005067 min^-1
        self.assertTrue(Cl>0.005 and Cl<0.007) # Gabrielsson p. 626: (ClR+Clm)=(0.3149+0.0540) h^-1=.006148 min^-1
        self.assertTrue(V>1.5 and V<1.7) # Gabrielsson p. 626: Vc=1.61
        self.assertTrue(Vp>0.15 and Vp<0.17) # Gabrielsson p. 626: Vt=0.1647
        self.assertTrue(fe>0.84 and fe<0.88) # Gabrielsson p. 626: Clr=0.3149 h^-1, fe=ClR/(ClR+Clm)=0.3149/(0.3149+0.0540)=0.8536
        fitting = PKPDFitting()
        fitting.load(protIVTwoCompartmentsUrine.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-120)

if __name__ == "__main__":
    unittest.main()
