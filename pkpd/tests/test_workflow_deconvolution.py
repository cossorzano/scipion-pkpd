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

    def testGabrielssonPK02Workflow(self):
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
        self.assertTrue(a0 > 0.63 and a0 < 0.69)
        self.assertTrue(a1 > 0.81 and a1 < 0.88)

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
        self.assertTrue(A[64] > 42 and A[64] < 46)
        self.assertTrue(A[200] > 90 and A[200] < 95)

if __name__ == "__main__":
    unittest.main()
