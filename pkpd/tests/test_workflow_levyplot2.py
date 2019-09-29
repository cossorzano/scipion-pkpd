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
import os

from pyworkflow.tests import *
from pkpd.protocols import *
from pkpd.objects import PKPDDataSet, PKPDExperiment
from test_workflow import TestWorkflow

def unitResponse(D,V,Ka,Cl,t):
    Ke=Cl/V
    t=np.clip(t,0.0,None)
    C=D/V*Ka/(Ka-Ke)*(np.exp(-Ke*t)-np.exp(-Ka*t))
    return C

class TestLevyPlotWorkflow2(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)

    def testDissolutionWorkflow(self):
        # Create invivo data
        experimentStr = """
[EXPERIMENT] ===========================
comment = 
title = Dissolution

[VARIABLES] ============================
C ; none ; numeric[%f] ; measurement ; Concentration in solution (%)
t ; min ; numeric[%f] ; time ; Time in minutes since start

[VIAS] ================================

[DOSES] ================================

[GROUPS] ================================
__Profile

[SAMPLES] ================================
Profile; group=__Profile

[MEASUREMENTS] ===========================
Profile ; t; C
0 0 
2.5 1.7 
5 8.3 
7.5 13.3 
10 20.0 
20 44.0 
30 61.0 
40 70.7 
60 78.0 
80 79.7 
100 80.7 
120 80.0 
160 81.3 
200 82.0 
240 82.3 
"""
        fnExperiment = "experimentInVitro.pkpd"
        fhExperiment = open(fnExperiment, "w")
        fhExperiment.write(experimentStr)
        fhExperiment.close()

        print "Import Experiment in vitro"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment in vitro',
                                      inputFile=fnExperiment)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        os.remove(fnExperiment)

        # Fit a Weibull dissolution
        print "Fitting Weibull model ..."
        protWeibull = self.newProtocol(ProtPKPDDissolutionFit,
                                objLabel='pkpd - fit dissolution Weibull',
                                globalSearch=True, modelType=3)
        protWeibull.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protWeibull)
        self.assertIsNotNone(protWeibull.outputExperiment.fnPKPD, "There was a problem with the dissolution model ")
        self.assertIsNotNone(protWeibull.outputFitting.fnFitting, "There was a problem with the dissolution model ")
        self.validateFiles('ProtPKPDDissolutionFit', ProtPKPDDissolutionFit)
        experiment = PKPDExperiment()
        experiment.load(protWeibull.outputExperiment.fnPKPD)
        Vmax = float(experiment.samples['Profile'].descriptors['Vmax'])
        self.assertTrue(Vmax>80 and Vmax<82)
        lambdda = float(experiment.samples['Profile'].descriptors['lambda'])
        self.assertTrue(lambdda>0.009 and lambdda<0.011)
        b = float(experiment.samples['Profile'].descriptors['b'])
        self.assertTrue(b>1.4 and b<1.5)

        fitting = PKPDFitting()
        fitting.load(protWeibull.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.997)

        # Create invivo data
        experimentStr = """
[EXPERIMENT] ===========================
comment = Generated as C(t)=D0/V*Ka/(Ka-Ke)*)(exp(-Ke*t)-exp(-Ka*t))
title = My experiment

[VARIABLES] ============================
Cp ; ug/L ; numeric[%f] ; measurement ; Plasma concentration
t ; min ; numeric[%f] ; time ; 

[VIAS] ================================
Oral; splineXY5;  tlag min; bioavailability=1.000000

[DOSES] ================================
Bolus1; via=Oral; bolus; t=0.000000 h; d=200 ug

[GROUPS] ================================
__Individual1

[SAMPLES] ================================
Individual1; dose=Bolus1; group=__Individual1

[MEASUREMENTS] ===========================
Individual1 ; t; Cp
"""
        t = np.arange(0,1000,10)
        Cp = unitResponse(100,50,0.05,0.2,t-20)+unitResponse(100,50,0.05,0.2,t-120)
        for n in range(t.size):
            experimentStr+="%f %f\n"%(t[n],Cp[n])
        fnExperiment ="experimentInVivo.pkpd"
        fhExperiment = open(fnExperiment,"w")
        fhExperiment.write(experimentStr)
        fhExperiment.close()

        print "Import Experiment in vivo"
        protImportInVivo = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment in vivo',
                                      inputFile=fnExperiment)
        self.launchProtocol(protImportInVivo)
        self.assertIsNotNone(protImportInVivo.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImportInVivo)

        os.remove(fnExperiment)

        # NCA numeric
        print "NCA numeric ..."
        protNCA = self.newProtocol(ProtPKPDNCANumeric,
                                objLabel='pkpd - nca numeric')
        protNCA.inputExperiment.set(protImportInVivo.outputExperiment)
        self.launchProtocol(protNCA)
        self.assertIsNotNone(protNCA.outputExperiment.fnPKPD, "There was a problem with the NCA numeric")
        self.validateFiles('prot', protNCA)
        experiment = PKPDExperiment()
        experiment.load(protNCA.outputExperiment.fnPKPD)
        AUC0t = float(experiment.samples['Individual1'].descriptors['AUC0t'])
        self.assertTrue(AUC0t > 960 and AUC0t < 980)
        AUMC0t = float(experiment.samples['Individual1'].descriptors['AUMC0t'])
        self.assertTrue(AUMC0t > 305000 and AUMC0t < 306000)
        Cmax = float(experiment.samples['Individual1'].descriptors['Cmax'])
        self.assertTrue(Cmax > 2.7 and Cmax < 2.9)
        Tmax = float(experiment.samples['Individual1'].descriptors['Tmax'])
        self.assertTrue(Tmax > 155 and Tmax < 165)
        MRT = float(experiment.samples['Individual1'].descriptors['MRT'])
        self.assertTrue(MRT > 314 and MRT < 315)

        # Fit Order 1
        print "Fitting splines5-monocompartment model ..."
        protModelInVivo = self.newProtocol(ProtPKPDMonoCompartment,
                                       objLabel='pkpd - fit monocompartment',
                                       bounds="(15.0, 30.0); (0.0, 400.0); (0.0, 1.0); (0.0, 1.0); (0.0, 1.0); (0.0, 1.0); (0.0, 1.0); (0.0, 1.0); (0.0, 1.0); (0.0, 1.0); (0.0, 1.0); (0.0, 1.0); (0.15, 0.25); (47, 53)"
                                       )
        protModelInVivo.inputExperiment.set(protImportInVivo.outputExperiment)
        self.launchProtocol(protModelInVivo)
        self.assertIsNotNone(protModelInVivo.outputExperiment.fnPKPD, "There was a problem with the PK model")
        self.assertIsNotNone(protModelInVivo.outputFitting.fnFitting, "There was a problem with the PK model")
        self.validateFiles('ProtPKPDMonoCompartment', ProtPKPDMonoCompartment)

        experiment = PKPDExperiment()
        experiment.load(protModelInVivo.outputExperiment.fnPKPD)
        V = float(experiment.samples['Individual1'].descriptors['V'])
        self.assertTrue(V>48 and V<52)
        Cl = float(experiment.samples['Individual1'].descriptors['Cl'])
        self.assertTrue(Cl>0.19 and Cl<0.21)

        fitting = PKPDFitting()
        fitting.load(protModelInVivo.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.998)

        # Deconvolve the in vivo
        print "Deconvolving in vivo ..."
        protDeconv = self.newProtocol(ProtPKPDDeconvolve,
                                       objLabel='pkpd - deconvolution'
                                       )
        protDeconv.inputODE.set(protModelInVivo)
        self.launchProtocol(protDeconv)
        self.assertIsNotNone(protDeconv.outputExperiment.fnPKPD, "There was a problem with the deconvolution")
        self.validateFiles('ProtPKPDDeconvolve', ProtPKPDDeconvolve)

        # Levy plot
        print "Levy plot ..."
        protLevy = self.newProtocol(ProtPKPDDissolutionLevyPlot,
                                      objLabel='pkpd - levy plot'
                                      )
        protLevy.inputInVitro.set(protWeibull)
        protLevy.inputInVivo.set(protDeconv)
        self.launchProtocol(protLevy)
        self.assertIsNotNone(protLevy.outputExperiment.fnPKPD, "There was a problem with the Levy plot")
        self.validateFiles('ProtPKPDDissolutionLevyPlot', ProtPKPDDissolutionLevyPlot)

        # IVIVC
        print "In vitro-in vivo correlation ..."
        protIVIVC = self.newProtocol(ProtPKPDDissolutionIVIVC,
                                     timeScale=5,
                                     responseScale=1,
                                     objLabel='pkpd - ivivc'
                                    )
        protIVIVC.inputInVitro.set(protWeibull)
        protIVIVC.inputInVivo.set(protDeconv)
        self.launchProtocol(protIVIVC)
        self.assertIsNotNone(protIVIVC.outputExperimentFabs.fnPKPD, "There was a problem with the IVIVC")
        self.assertIsNotNone(protIVIVC.outputExperimentAdissol.fnPKPD, "There was a problem with the IVIVC")
        self.validateFiles('ProtPKPDDissolutionIVIVC', ProtPKPDDissolutionIVIVC)

        # IVIVC generic
        print "In vitro-in vivo generic ..."
        protIVIVCG = self.newProtocol(ProtPKPDDissolutionIVIVCGeneric,
                                      timeScale='$[k1]*$(t)+$[k2]*np.power($(t),2)+$[k3]*np.power($(t),3)',
                                      timeBounds='k1: [0,3]; k2: [-0.1,0.01];  k3: [0,1e-3]',
                                      responseScale='$[A]*$(Adissol)+$[B]+$[C]*np.power($(Adissol),2)',
                                      responseBounds='A: [0.01,1]; B: [-50,30]; C: [-0.05,0.05]',
                                      objLabel='pkpd - ivivc generic'
                                     )
        protIVIVCG.inputInVitro.set(protWeibull)
        protIVIVCG.inputInVivo.set(protDeconv)
        self.launchProtocol(protIVIVCG)
        self.assertIsNotNone(protIVIVCG.outputExperimentFabs.fnPKPD, "There was a problem with the IVIVC Generic")
        self.validateFiles('ProtPKPDDissolutionIVIVCG', ProtPKPDDissolutionIVIVCGeneric)

        # IVIVC splines
        print "In vitro-in vivo splies ..."
        protIVIVCS = self.newProtocol(ProtPKPDDissolutionIVIVCSplines,
                                      timeScale=1,
                                      responseScale=1,
                                      objLabel='pkpd - ivivc splines'
                                     )
        protIVIVCS.inputInVitro.set(protWeibull)
        protIVIVCS.inputInVivo.set(protDeconv)
        self.launchProtocol(protIVIVCS)
        self.assertIsNotNone(protIVIVCS.outputExperimentFabs.fnPKPD, "There was a problem with the IVIVC Splines")
        self.validateFiles('ProtPKPDDissolutionIVIVCS', ProtPKPDDissolutionIVIVCSplines)

        # Dissolution simulation
        print "IVIV+PK simulation ..."
        protIVIVPKL = self.newProtocol(ProtPKPDDissolutionPKSimulation,
                                      objLabel='pkpd - ivivc+pk',
                                      conversionType=1,
                                      inputN=1,
                                      tF=16.66,
                                      addIndividuals=True,
                                      inputDose=200
                                      )
        protIVIVPKL.inputInVitro.set(protWeibull.outputFitting)
        protIVIVPKL.inputPK.set(protModelInVivo.outputFitting)
        protIVIVPKL.inputLevy.set(protLevy.outputExperiment)
        self.launchProtocol(protIVIVPKL)
        self.assertIsNotNone(protIVIVPKL.outputExperiment.fnPKPD, "There was a problem with the simulation")
        self.validateFiles('ProtPKPDDissolutionPKSimulation', ProtPKPDDissolutionPKSimulation)

        # Dissolution simulation
        print "IVIV+PK simulation ..."
        protIVIVPKS = self.newProtocol(ProtPKPDDissolutionPKSimulation,
                                       objLabel='pkpd - ivivc+pk',
                                       inputN=1,
                                       tF=16.66,
                                       addIndividuals=True,
                                       inputDose=200
                                       )
        protIVIVPKS.inputInVitro.set(protWeibull.outputFitting)
        protIVIVPKS.inputPK.set(protModelInVivo.outputFitting)
        protIVIVPKS.inputIvIvC.set(protIVIVCS.outputExperimentFabs)
        self.launchProtocol(protIVIVPKS)
        self.assertIsNotNone(protIVIVPKS.outputExperiment.fnPKPD, "There was a problem with the simulation")
        self.validateFiles('ProtPKPDDissolutionPKSimulation', ProtPKPDDissolutionPKSimulation)

        # Internal validity
        print "Internal validity ..."
        protInternal = self.newProtocol(ProtPKPDIVIVCInternalValidity,
                                        objLabel='pkpd - internal validity')
        protInternal.inputExperiment.set(protNCA.outputExperiment)
        protInternal.inputSimulated.set(protIVIVPKL.outputExperiment)
        self.launchProtocol(protInternal)
        fnSummary = protInternal._getPath("summary.txt")
        self.assertTrue(os.path.exists(fnSummary))
        lineNo = 0
        for line in open(fnSummary).readlines():
            tokens = line.split('=')
            if lineNo == 0:
                AUCmean = np.abs(float(tokens[-1]))
                self.assertTrue(AUCmean < 20)
            elif lineNo == 1:
                Cmaxmean = np.abs(float(tokens[-1]))
                self.assertTrue(Cmaxmean < 10)
            lineNo += 1

        # Internal validity
        print "Internal validity ..."
        protInternal = self.newProtocol(ProtPKPDIVIVCInternalValidity,
                                    objLabel='pkpd - internal validity')
        protInternal.inputExperiment.set(protNCA.outputExperiment)
        protInternal.inputSimulated.set(protIVIVPKS.outputExperiment)
        self.launchProtocol(protInternal)
        fnSummary = protInternal._getPath("summary.txt")
        self.assertTrue(os.path.exists(fnSummary))
        lineNo = 0
        for line in open(fnSummary).readlines():
            tokens = line.split('=')
            if lineNo == 0:
                AUCmean = np.abs(float(tokens[-1]))
                self.assertTrue(AUCmean < 10)
            elif lineNo == 1:
                Cmaxmean = np.abs(float(tokens[-1]))
                self.assertTrue(Cmaxmean < 20)
            lineNo += 1

if __name__ == "__main__":
    unittest.main()
