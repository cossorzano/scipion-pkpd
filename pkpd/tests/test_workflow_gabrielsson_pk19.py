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


class TestGabrielssonPK19Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK19')
        cls.exptFn = cls.dataset.getFile('experiment')
        cls.expIFn = cls.dataset.getFile('experimentIndividual')

    def testGabrielssonPK19Workflow(self):
        # Import an experiment (intravenous)

        print "Import Experiment (intravenous doses) Group"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment Group',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        print "Fitting a two-compartments model intrinsic with metabolite ..."
        protPKPDTwoCompartment = self.newProtocol(ProtPKPDTwoCompartmentsClintMetabolite,
                                                     objLabel='pkpd - iv two-compartments intrinsic, metabolite, group',
                                                     globalSearch=False,
                                                     bounds='(0.5, 3.0); (15.0, 200.0); (1.0, 5.0); (0.001, 0.4); (0.5, 5.0); (0.01, 0.5); (0.1, 4.0)')
        protPKPDTwoCompartment.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPKPDTwoCompartment)
        self.assertIsNotNone(protPKPDTwoCompartment.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protPKPDTwoCompartment.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDTwoCompartment', protPKPDTwoCompartment)
        # Individual10; Clm=0.0412476111466; Clp=0.0752579955459; Km=39.8444811617; V=1.46970385411; Vm=0.327649681816; Vmax=1.4155764387; Vp=2.24003743827
        experiment = PKPDExperiment()
        experiment.load(protPKPDTwoCompartment.outputExperiment.fnPKPD)
        Vmax = float(experiment.samples['Individual10'].descriptors['Vmax'])
        Clm = float(experiment.samples['Individual10'].descriptors['Clm'])
        Vm = float(experiment.samples['Individual10'].descriptors['Vm'])
        Km = float(experiment.samples['Individual10'].descriptors['Km'])
        V = float(experiment.samples['Individual10'].descriptors['V'])
        Clp = float(experiment.samples['Individual10'].descriptors['Clp'])
        Vp = float(experiment.samples['Individual10'].descriptors['Vp'])
        self.assertTrue(Vmax>1.3 and Vmax<1.5) # Gabrielsson p. 655: Vmax=1.69
        self.assertTrue(Clm>0.035 and Clm<0.045) # Gabrielsson p. 655: kme=0.14
        self.assertTrue(Vm>0.25 and Vm<0.37) # Gabrielsson p. 655: Vme=0.29
        self.assertTrue(Km>37 and Km<43) # Gabrielsson p. 655: Km=57
        self.assertTrue(V>1.4 and V<1.5) # Gabrielsson p. 655: V=1.07
        self.assertTrue(Clp>0.07 and Clp<0.08) # Gabrielsson p. 655: Cld=0.128
        self.assertTrue(Vp>2 and Vp<2.5) # Gabrielsson p. 655: Vt=1.99
        fitting = PKPDFitting()
        fitting.load(protPKPDTwoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.98)

        print "Fitting a two-compartments model intrinsic ..."
        protPKPDTwoCompartment = self.newProtocol(ProtPKPDTwoCompartmentsClint,
                                                     objLabel='pkpd - iv two-compartments intrinsic, group',
                                                     globalSearch=False,
                                                     bounds='(0.1, 1.5); (0.1, 30.0); (0.1, 3.0); (0.01, 0.15); (0.1, 5.0)')
        protPKPDTwoCompartment.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPKPDTwoCompartment)
        self.assertIsNotNone(protPKPDTwoCompartment.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protPKPDTwoCompartment.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDTwoCompartment', protPKPDTwoCompartment)
        # Individual10; Clp=0.0656654314355; Km=17.4883311291; V=1.49956827889; Vmax=0.75629408911; Vp=2.54888209969
        experiment = PKPDExperiment()
        experiment.load(protPKPDTwoCompartment.outputExperiment.fnPKPD)
        Vmax = float(experiment.samples['Individual10'].descriptors['Vmax'])
        Km = float(experiment.samples['Individual10'].descriptors['Km'])
        V = float(experiment.samples['Individual10'].descriptors['V'])
        Clp = float(experiment.samples['Individual10'].descriptors['Clp'])
        Vp = float(experiment.samples['Individual10'].descriptors['Vp'])
        self.assertTrue(Vmax>0.7 and Vmax<0.8) # Gabrielsson p. 655: Vmax=1.69
        self.assertTrue(Km>15 and Km<19) # Gabrielsson p. 655: Km=57
        self.assertTrue(V>1.4 and V<1.6) # Gabrielsson p. 655: V=1.07
        self.assertTrue(Clp>0.06 and Clp<0.07) # Gabrielsson p. 655: Cld=0.128
        self.assertTrue(Vp>2.25 and Vp<2.75) # Gabrielsson p. 655: Vt=1.99
        fitting = PKPDFitting()
        fitting.load(protPKPDTwoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Fit a mono-compartment model with intravenous absorption to a set of measurements
        print "Fitting a two-compartments model intrinsic with metabolite, group, linear ..."
        protPKPDTwoCompartment = self.newProtocol(ProtPKPDTwoCompartmentsClintMetabolite,
                                                     objLabel='pkpd - iv two-compartments intrinsic, metabolite, group, linear',
                                                     globalSearch=False,
                                                     fitType=0,
                                                     bounds='(0.5, 3.0); (15.0, 200.0); (1.0, 5.0); (0.001, 0.4); (0.5, 5.0); (0.01, 0.5); (0.1, 4.0)')
        protPKPDTwoCompartment.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPKPDTwoCompartment)
        self.assertIsNotNone(protPKPDTwoCompartment.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protPKPDTwoCompartment.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDTwoCompartment', protPKPDTwoCompartment)
        # Individual10; Clm=0.0432221623361; Clp=0.127294496616; Km=96.3645307085; V=0.975045098718; Vm=0.396948904704; Vmax=2.30837989545; Vp=1.9212940747
        experiment = PKPDExperiment()
        experiment.load(protPKPDTwoCompartment.outputExperiment.fnPKPD)
        Vmax = float(experiment.samples['Individual10'].descriptors['Vmax'])
        Clm = float(experiment.samples['Individual10'].descriptors['Clm'])
        Vm = float(experiment.samples['Individual10'].descriptors['Vm'])
        Km = float(experiment.samples['Individual10'].descriptors['Km'])
        V = float(experiment.samples['Individual10'].descriptors['V'])
        Clp = float(experiment.samples['Individual10'].descriptors['Clp'])
        Vp = float(experiment.samples['Individual10'].descriptors['Vp'])
        self.assertTrue(Vmax>2.1 and Vmax<2.4) # Gabrielsson p. 655: Vmax=1.69
        self.assertTrue(Clm>0.038 and Clm<0.048) # Gabrielsson p. 655: kme=0.14
        self.assertTrue(Vm>0.35 and Vm<0.45) # Gabrielsson p. 655: Vme=0.29
        self.assertTrue(Km>85 and Km<100) # Gabrielsson p. 655: Km=57
        self.assertTrue(V>0.9 and V<1.1) # Gabrielsson p. 655: V=1.07
        self.assertTrue(Clp>0.11 and Clp<0.15) # Gabrielsson p. 655: Cld=0.128
        self.assertTrue(Vp>1.8 and Vp<2) # Gabrielsson p. 655: Vt=1.99
        fitting = PKPDFitting()
        fitting.load(protPKPDTwoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        print "Import Experiment (intravenous doses) Individual"
        protImportIndividual = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment Individual',
                                      inputFile=self.expIFn)
        self.launchProtocol(protImportIndividual)
        self.assertIsNotNone(protImportIndividual.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImportIndividual', protImportIndividual)

        print "Fitting a two-compartments model intrinsic with metabolite, individual ..."
        protPKPDTwoCompartment = self.newProtocol(ProtPKPDTwoCompartmentsClintMetabolite,
                                                     objLabel='pkpd - iv two-compartments intrinsic, metabolite',
                                                     globalSearch=False,
                                                     bounds='(0.5, 3.0); (15.0, 200.0); (1.0, 5.0); (0.001, 0.4); (0.5, 5.0); (0.01, 0.5); (0.1, 4.0)')
        protPKPDTwoCompartment.inputExperiment.set(protImportIndividual.outputExperiment)
        self.launchProtocol(protPKPDTwoCompartment)
        self.assertIsNotNone(protPKPDTwoCompartment.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protPKPDTwoCompartment.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDTwoCompartment', protPKPDTwoCompartment)
        # Individual10; dose=Infusion10; group=__Individual10 ; Clm=0.0499683694388; Clp=0.238822441847; Km=24.3049234128; V=1.00000275169; Vm=0.640665157721; Vmax=1.28036830695; Vp=1.97543492947
        # Individual300; dose=Infusion300; group=__Individual300 ; Clm=0.0443708277703; Clp=0.117178521139; Km=154.25975031; V=1.02601824853; Vm=0.430637326901; Vmax=2.99941542931; Vp=1.86417420681
        # Individual50; dose=Infusion50; group=__Individual50 ; Clm=0.0402926886276; Clp=0.0668657643453; Km=16.9762033784; V=1.48256537935; Vm=0.195238425892; Vmax=0.706407025484; Vp=2.72796064002
        experiment = PKPDExperiment()
        experiment.load(protPKPDTwoCompartment.outputExperiment.fnPKPD)
        Vmax = float(experiment.samples['Individual10'].descriptors['Vmax'])
        Clm = float(experiment.samples['Individual10'].descriptors['Clm'])
        Vm = float(experiment.samples['Individual10'].descriptors['Vm'])
        Km = float(experiment.samples['Individual10'].descriptors['Km'])
        V = float(experiment.samples['Individual10'].descriptors['V'])
        Clp = float(experiment.samples['Individual10'].descriptors['Clp'])
        Vp = float(experiment.samples['Individual10'].descriptors['Vp'])
        self.assertTrue(Vmax>0.65 and Vmax<3.2) # Gabrielsson p. 655: Vmax=1.69
        self.assertTrue(Clm>0.035 and Clm<0.055) # Gabrielsson p. 655: kme=0.14
        self.assertTrue(Vm>0.17 and Vm<0.7) # Gabrielsson p. 655: Vme=0.29
        self.assertTrue(Km>10 and Km<160) # Gabrielsson p. 655: Km=57
        self.assertTrue(V>0.95 and V<1.6) # Gabrielsson p. 655: V=1.07
        self.assertTrue(Clp>0.05 and Clp<0.28) # Gabrielsson p. 655: Cld=0.128
        self.assertTrue(Vp>1.7 and Vp<2.9) # Gabrielsson p. 655: Vt=1.99

        Vmax = float(experiment.samples['Individual50'].descriptors['Vmax'])
        Clm = float(experiment.samples['Individual50'].descriptors['Clm'])
        Vm = float(experiment.samples['Individual50'].descriptors['Vm'])
        Km = float(experiment.samples['Individual50'].descriptors['Km'])
        V = float(experiment.samples['Individual50'].descriptors['V'])
        Clp = float(experiment.samples['Individual50'].descriptors['Clp'])
        Vp = float(experiment.samples['Individual50'].descriptors['Vp'])
        self.assertTrue(Vmax>0.65 and Vmax<3.2) # Gabrielsson p. 655: Vmax=1.69
        self.assertTrue(Clm>0.035 and Clm<0.055) # Gabrielsson p. 655: kme=0.14
        self.assertTrue(Vm>0.17 and Vm<0.7) # Gabrielsson p. 655: Vme=0.29
        self.assertTrue(Km>10 and Km<160) # Gabrielsson p. 655: Km=57
        self.assertTrue(V>0.95 and V<1.6) # Gabrielsson p. 655: V=1.07
        self.assertTrue(Clp>0.05 and Clp<0.28) # Gabrielsson p. 655: Cld=0.128
        self.assertTrue(Vp>1.7 and Vp<2.9) # Gabrielsson p. 655: Vt=1.99

        Vmax = float(experiment.samples['Individual300'].descriptors['Vmax'])
        Clm = float(experiment.samples['Individual300'].descriptors['Clm'])
        Vm = float(experiment.samples['Individual300'].descriptors['Vm'])
        Km = float(experiment.samples['Individual300'].descriptors['Km'])
        V = float(experiment.samples['Individual300'].descriptors['V'])
        Clp = float(experiment.samples['Individual300'].descriptors['Clp'])
        Vp = float(experiment.samples['Individual300'].descriptors['Vp'])
        self.assertTrue(Vmax>0.65 and Vmax<3.2) # Gabrielsson p. 655: Vmax=1.69
        self.assertTrue(Clm>0.035 and Clm<0.055) # Gabrielsson p. 655: kme=0.14
        self.assertTrue(Vm>0.17 and Vm<0.7) # Gabrielsson p. 655: Vme=0.29
        self.assertTrue(Km>10 and Km<160) # Gabrielsson p. 655: Km=57
        self.assertTrue(V>0.95 and V<1.6) # Gabrielsson p. 655: V=1.07
        self.assertTrue(Clp>0.05 and Clp<0.28) # Gabrielsson p. 655: Cld=0.128
        self.assertTrue(Vp>1.7 and Vp<2.9) # Gabrielsson p. 655: Vt=1.99

        fitting = PKPDFitting()
        fitting.load(protPKPDTwoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.97)
        self.assertTrue(fitting.sampleFits[1].R2>0.97)
        self.assertTrue(fitting.sampleFits[2].R2>0.97)

        print "Fitting a two-compartments model, individual ..."
        protPKPDTwoCompartment = self.newProtocol(ProtPKPDTwoCompartments,
                                                     objLabel='pkpd - iv two-compartments',
                                                     globalSearch=False,
                                                     bounds='(0.0, 0.1); (0.0, 3.0); (0.0, 0.25); (0.0, 4.0)')
        protPKPDTwoCompartment.inputExperiment.set(protImportIndividual.outputExperiment)
        self.launchProtocol(protPKPDTwoCompartment)
        self.assertIsNotNone(protPKPDTwoCompartment.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protPKPDTwoCompartment.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDTwoCompartment', protPKPDTwoCompartment)
        # Individual10; dose=Infusion10; group=__Individual10 ; Cl=0.04925267021; Clp=0.127978254678; V=1.48992920408; Vp=1.5007459044
        # Individual300; dose=Infusion300; group=__Individual300 ; Cl=0.0131117993155; Clp=0.107585836679; V=1.04464809241; Vp=1.80698471861
        # Individual50; dose=Infusion50; group=__Individual50 ; Cl=0.0293448221159; Clp=0.0896564329198; V=1.1715951472; Vp=2.15610019635
        experiment = PKPDExperiment()
        experiment.load(protPKPDTwoCompartment.outputExperiment.fnPKPD)
        V = float(experiment.samples['Individual10'].descriptors['V'])
        Cl = float(experiment.samples['Individual10'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual10'].descriptors['Clp'])
        Vp = float(experiment.samples['Individual10'].descriptors['Vp'])
        self.assertTrue(V>0.9 and V<1.6) # Gabrielsson p. 655: V=1.07
        self.assertTrue(Cl>0.009 and Cl<0.055)
        self.assertTrue(Clp>0.07 and Clp<0.15) # Gabrielsson p. 655: Cld=0.128
        self.assertTrue(Vp>1.4 and Vp<2.3) # Gabrielsson p. 655: Vt=1.99

        V = float(experiment.samples['Individual50'].descriptors['V'])
        Cl = float(experiment.samples['Individual50'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual50'].descriptors['Clp'])
        Vp = float(experiment.samples['Individual50'].descriptors['Vp'])
        self.assertTrue(V>0.9 and V<1.6) # Gabrielsson p. 655: V=1.07
        self.assertTrue(Cl>0.009 and Cl<0.055)
        self.assertTrue(Clp>0.07 and Clp<0.15) # Gabrielsson p. 655: Cld=0.128
        self.assertTrue(Vp>1.4 and Vp<2.3) # Gabrielsson p. 655: Vt=1.99

        V = float(experiment.samples['Individual300'].descriptors['V'])
        Cl = float(experiment.samples['Individual300'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual300'].descriptors['Clp'])
        Vp = float(experiment.samples['Individual300'].descriptors['Vp'])
        self.assertTrue(V>0.9 and V<1.6) # Gabrielsson p. 655: V=1.07
        self.assertTrue(Cl>0.009 and Cl<0.055)
        self.assertTrue(Clp>0.07 and Clp<0.15) # Gabrielsson p. 655: Cld=0.128
        self.assertTrue(Vp>1.4 and Vp<2.3) # Gabrielsson p. 655: Vt=1.99

        fitting = PKPDFitting()
        fitting.load(protPKPDTwoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.97)
        self.assertTrue(fitting.sampleFits[1].R2>0.97)
        self.assertTrue(fitting.sampleFits[2].R2>0.97)

if __name__ == "__main__":
    unittest.main()
