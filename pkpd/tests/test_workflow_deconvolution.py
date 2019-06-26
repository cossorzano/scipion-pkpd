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

import numpy as np
import unittest
from pyworkflow.em import *
from pyworkflow.tests import *
from pkpd.protocols import *
from test_workflow import TestWorkflow
from pkpd.objects import PKPDDataSet


class TestDeconvolutionWorkflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK02')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testDeconvolutionWorkflow(self):
        #First, import an experiment

        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)      

        # Fit a monocompartmental model with first order absorption
        print "Fitting monocompartmental model..."
        protEV1MonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                  objLabel='pkpd - ev1 monocompartment',
                                                  bounds='(0.0, 20.0); (0.0, 0.2); (0.0, 1.0); (0.0, 100.0)')
        protEV1MonoCompartment.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)

        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual1'].descriptors['Cl'])
        V = float(experiment.samples['Individual1'].descriptors['V'])
        Ka = float(experiment.samples['Individual1'].descriptors['Oral_Ka'])
        tlag = float(experiment.samples['Individual1'].descriptors['Oral_tlag'])
        self.assertTrue(Cl>0.27 and Cl<0.32) # Gabrielsson p 515, Solution II: CL/F=0.2819
        self.assertTrue(V>10 and V<40) # Gabrielsson p 515, Solution II: V/F=32.05 -------------- Mine: 27.5
        self.assertTrue(Ka>0.015 and Ka<0.06) # Gabrielsson p 511, Solution II: Ka=0.043 -------- Mine: 0.0264
        self.assertTrue(tlag>10 and tlag<25) # Gabrielsson p 511, Solution II: tlag=16

        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.93)
        self.assertTrue(fitting.sampleFits[0].AIC<-14)

        # Deconvolution
        print "Deconvolving ..."
        prot = self.newProtocol(ProtPKPDDeconvolve,
                                objLabel='dissol deconv')
        prot.inputODE.set(protEV1MonoCompartment)
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputExperiment.fnPKPD, "There was a problem with the deconvolution")
        self.validateFiles('prot', prot)
        experiment = PKPDExperiment()
        experiment.load(prot.outputExperiment.fnPKPD)
        A = np.asarray(experiment.samples['Individual1'].getValues('A'),dtype=np.float64)
        self.assertTrue(A[0]==0.0)
        self.assertTrue(A[64]>42 and A[64]<46)
        self.assertTrue(A[200]>90 and A[200]<95)


        # Simulate PK
        print "Simulate PK ..."
        prot = self.newProtocol(ProtPKPDODESimulate,
                                objLabel='PK simulate',
                                paramsSource=1,
                                prmUser='0.0, 0.035339, 0.2835, 28.7765',
                                doses='Bolus ; via=Oral; bolus; t=0 h; d=1 ug',
                                tF=24)
        prot.inputODE.set(protEV1MonoCompartment)
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputExperiment.fnPKPD, "There was a problem with the deconvolution")
        self.validateFiles('prot', prot)
        experiment = PKPDExperiment()
        experiment.load(prot.outputExperiment.fnPKPD)
        AUC0t = float(experiment.samples['Simulation_0'].descriptors['AUC0t'])
        self.assertTrue(AUC0t > 3.4 and AUC0t < 3.6)
        AUMC0t = float(experiment.samples['Simulation_0'].descriptors['AUMC0t'])
        self.assertTrue(AUMC0t > 455 and AUMC0t < 457)

        # Fit a monocompartmental model with first order absorption
        print "Fitting monocompartmental model..."
        protEV1MonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                  objLabel='pkpd - ev1 monocompartment simulated',
                                                  bounds='(0.0, 20.0); (0.0, 0.2); (0.0, 1.0); (0.0, 100.0)')
        protEV1MonoCompartment.inputExperiment.set(prot.outputExperiment)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD,
                             "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartment.outputFitting.fnFitting,
                             "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)

        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Simulation_0'].descriptors['Cl'])
        V = float(experiment.samples['Simulation_0'].descriptors['V'])
        Ka = float(experiment.samples['Simulation_0'].descriptors['Oral_Ka'])
        tlag = float(experiment.samples['Simulation_0'].descriptors['Oral_tlag'])
        self.assertTrue(Cl > 0.26 and Cl < 0.30)
        self.assertTrue(V > 26 and V < 30)
        self.assertTrue(tlag > 0 and tlag < 3)
        self.assertTrue(Ka > 0.03 and Ka < 0.042)

        # Deconvolution Fourier
        print "Deconvolving Fourier ..."
        prot = self.newProtocol(ProtPKPDDeconvolveFourier,
                                objLabel='dissol deconv Fourier')
        prot.inputODE.set(protEV1MonoCompartment)
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputExperiment.fnPKPD, "There was a problem with the deconvolution")
        self.validateFiles('prot', prot)
        experiment = PKPDExperiment()
        experiment.load(prot.outputExperiment.fnPKPD)
        A = np.asarray(experiment.samples['Simulation_0'].getValues('A'), dtype=np.float64)
        self.assertTrue(A[400] > 97 and A[400] <= 100)

        # Change via
        print "Change via ..."
        protVia = self.newProtocol(ProtPKPDChangeVia,
                                objLabel='change via - spline',
                                viaName='Oral',newViaType='spline2')
        protVia.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protVia)
        self.assertIsNotNone(protVia.outputExperiment.fnPKPD, "There was a problem with the change of via")
        self.validateFiles('protVia', protVia)


        # Fit a monocompartmental model with first order absorption
        print "Fitting monocompartmental model..."
        protEV1MonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                  objLabel='pkpd - ev1 monocompartment spline',
                                                  bounds='(0.0, 20.0); (0.9, 1.0); (0.0, 100.0); (0.0, 1.0); (0.0, 1.0); (0.2, 0.4); (10.0, 40.0)')
        protEV1MonoCompartment.inputExperiment.set(protVia.outputExperiment)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD,
                             "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartment.outputFitting.fnFitting,
                             "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)

        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual1'].descriptors['Cl'])
        V = float(experiment.samples['Individual1'].descriptors['V'])
        a0 = float(experiment.samples['Individual1'].descriptors['Oral_spline2_A0'])
        a1 = float(experiment.samples['Individual1'].descriptors['Oral_spline2_A1'])
        tlag = float(experiment.samples['Individual1'].descriptors['Oral_tlag'])
        self.assertTrue(Cl > 0.25 and Cl < 0.32)  # Gabrielsson p 515, Solution II: CL/F=0.2819
        self.assertTrue(V > 10 and V < 40)  # Gabrielsson p 515, Solution II: V/F=32.05 -------------- Mine: 27.5
        self.assertTrue(tlag > 0 and tlag < 25)  # Gabrielsson p 511, Solution II: tlag=16
        self.assertTrue(a0 > 0.63 and a0 < 0.86)
        self.assertTrue(a1 > 0.2 and a1 < 0.88)

        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2 > 0.94)
        self.assertTrue(fitting.sampleFits[0].AIC < -14)

        # Deconvolution
        print "Deconvolving ..."
        prot = self.newProtocol(ProtPKPDDeconvolve,
                                objLabel='dissol deconv spline')
        prot.inputODE.set(protEV1MonoCompartment)
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputExperiment.fnPKPD, "There was a problem with the deconvolution")
        self.validateFiles('prot', prot)
        experiment = PKPDExperiment()
        experiment.load(prot.outputExperiment.fnPKPD)
        A = np.asarray(experiment.samples['Individual1'].getValues('A'), dtype=np.float64)
        self.assertTrue(A[0] == 0.0)
        self.assertTrue(A[64] > 42 and A[64] < 55)
        self.assertTrue(A[200] > 90 and A[200] <=100)

        # NCA numeric
        print "NCA numeric ..."
        prot = self.newProtocol(ProtPKPDNCANumeric,
                                objLabel='nca numeric')
        prot.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputExperiment.fnPKPD, "There was a problem with the deconvolution")
        self.validateFiles('prot', prot)
        experiment = PKPDExperiment()
        experiment.load(prot.outputExperiment.fnPKPD)
        AUC0t = float(experiment.samples['Individual1'].descriptors['AUC0t'])
        self.assertTrue(AUC0t > 325.5 and AUC0t < 327.5)
        AUMC0t = float(experiment.samples['Individual1'].descriptors['AUMC0t'])
        self.assertTrue(AUMC0t > 42165 and AUMC0t < 42170)
        Cmax = float(experiment.samples['Individual1'].descriptors['Cmax'])
        self.assertTrue(Cmax > 1.9 and Cmax < 2.1)
        Tmax = float(experiment.samples['Individual1'].descriptors['Tmax'])
        self.assertTrue(Tmax > 39 and Tmax < 41)
        MRT = float(experiment.samples['Individual1'].descriptors['MRT'])
        self.assertTrue(MRT > 129 and MRT < 130)

if __name__ == "__main__":
    unittest.main()
