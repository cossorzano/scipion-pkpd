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

from pyworkflow.tests import *
from pkpd.protocols import *
from pkpd.objects import PKPDDataSet
from test_workflow import TestWorkflow

class TestDeconvolution2Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Dissolution')
        cls.dissolutionFn = cls.dataset.getFile('exceldissolutionMean')
        cls.pkFn = cls.dataset.getFile('excelinvivoMean')

    def testDissolutionWorkflow(self):
        print "Import Dissolution"
        protImport = self.newProtocol(ProtPKPDImportFromTable,
                                      objLabel='pkpd - import dissolution',
                                      tVar='t; min; numeric; time; Time variable',
                                      inputFile=self.dissolutionFn)
        self.launchProtocol(protImport)
        self.assertTrue(os.path.exists(protImport.outputExperiment_Test.fnPKPD.get()), "There was a problem with the import")
        self.validateFiles('protImportTest', protImport)

        # Fit a Weibull dissolution
        print "Fitting Weibull model ..."
        protWeibull = self.newProtocol(ProtPKPDDissolutionFit,
                                objLabel='pkpd - fit dissolution Weibull',
                                allowTlag=True,bounds="(120,120);(80,100);(0,2);(0.1,4)",
                                globalSearch=True, modelType=3)
        protWeibull.inputExperiment.set(protImport.outputExperiment_Test)
        self.launchProtocol(protWeibull)
        self.assertIsNotNone(protWeibull.outputExperiment.fnPKPD, "There was a problem with the dissolution model ")
        self.assertIsNotNone(protWeibull.outputFitting.fnFitting, "There was a problem with the dissolution model ")
        self.validateFiles('ProtPKPDDissolutionFit', ProtPKPDDissolutionFit)
        experiment = PKPDExperiment()
        experiment.load(protWeibull.outputExperiment.fnPKPD)
        Vmax = float(experiment.samples['Mean'].descriptors['Vmax'])
        self.assertTrue(Vmax>97)
        lambdda = float(experiment.samples['Mean'].descriptors['lambda'])
        self.assertTrue(lambdda>0.04 and lambdda<0.08)
        b = float(experiment.samples['Mean'].descriptors['b'])
        self.assertTrue(b>0.56 and b<0.73)

        fitting = PKPDFitting()
        fitting.load(protWeibull.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.997)


        print "Import In vivo"
        protImport = self.newProtocol(ProtPKPDImportFromExcel,
                                      objLabel='pkpd - import invivo',
                                      variables='t;min;numeric;time; ;;Cp;ug/L;numeric;measurement;plasma concentration',
                                      vias='Oral;spline3',
                                      doses='Bolus;via=Oral;bolus;t=0 h;dose=500 ug',
                                      inputFile=self.pkFn)
        self.launchProtocol(protImport)
        self.assertTrue(os.path.exists(protImport.outputExperiment.fnPKPD.get()), "There was a problem with the import")
        self.validateFiles('protImportTest', protImport)

        # NCA numeric
        print "NCA numeric ..."
        protNCA = self.newProtocol(ProtPKPDNCANumeric,
                                objLabel='pkpd - nca numeric')
        protNCA.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protNCA)
        self.assertIsNotNone(protNCA.outputExperiment.fnPKPD, "There was a problem with the NCA numeric")
        self.validateFiles('prot', protNCA)
        experiment = PKPDExperiment()
        experiment.load(protNCA.outputExperiment.fnPKPD)
        AUC0t = float(experiment.samples['Mean'].descriptors['AUC0t'])
        self.assertTrue(AUC0t > 11600 and AUC0t < 11650)
        AUMC0t = float(experiment.samples['Mean'].descriptors['AUMC0t'])
        self.assertTrue(AUMC0t > 11e6 and AUMC0t < 11.5e6)
        Cmax = float(experiment.samples['Mean'].descriptors['Cmax'])
        self.assertTrue(Cmax > 13.5 and Cmax < 13.8)
        Tmax = float(experiment.samples['Mean'].descriptors['Tmax'])
        self.assertTrue(Tmax > 325 and Tmax < 335)
        MRT = float(experiment.samples['Mean'].descriptors['MRT'])
        self.assertTrue(MRT > 950 and MRT < 1000)

        # Fit a monocompartmental model with first order absorption
        print "Fitting two-compartments model..."
        protPK = self.newProtocol(ProtPKPDTwoCompartments,
                                                  objLabel='pkpd - twocompartments',
                                                  globalSearch=False,
                                                  bounds='(261.0, 350.0); (0.0, 0.2); (0.36, 0.56); (0.53, 0.73); (0.032, 0.052); (27.0, 29.0); (0.003, 0.023); (9.0, 11.0)')
        protPK.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPK)
        self.assertIsNotNone(protPK.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protPK.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protPK', protPK)

        experiment = PKPDExperiment()
        experiment.load(protPK.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Mean'].descriptors['Cl'])
        V = float(experiment.samples['Mean'].descriptors['V'])
        Clp = float(experiment.samples['Mean'].descriptors['Clp'])
        Vp = float(experiment.samples['Mean'].descriptors['Vp'])
        self.assertTrue(Cl>0.036 and Cl<0.046)
        self.assertTrue(V>27 and V<29)
        self.assertTrue(Clp>0.008 and Clp<0.018)
        self.assertTrue(Vp>9 and Vp<11)

        # Deconvolution
        print "Deconvolving ..."
        protDeconv = self.newProtocol(ProtPKPDDeconvolve,
                                objLabel='dissol deconv')
        protDeconv.inputODE.set(protPK)
        self.launchProtocol(protDeconv)
        self.assertIsNotNone(protDeconv.outputExperiment.fnPKPD, "There was a problem with the deconvolution")
        self.validateFiles('protDeconv', protDeconv)
        experiment = PKPDExperiment()
        experiment.load(protDeconv.outputExperiment.fnPKPD)
        A = np.asarray(experiment.samples['Mean'].getValues('A'),dtype=np.float64)
        self.assertTrue(A[0]<1)
        self.assertTrue(A[341]>49 and A[341]<51)
        self.assertTrue(A[606]>98)

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

        # Dissolution simulation
        print "IVIV+PK simulation ..."
        protIVIVPK = self.newProtocol(ProtPKPDDissolutionPKSimulation,
                                      objLabel='pkpd - ivivc+pk',
                                      conversionType=1,
                                      inputN=1,
                                      tF=72,
                                      addIndividuals=True,
                                      inputDose=500
                                      )
        protIVIVPK.inputInVitro.set(protWeibull.outputFitting)
        protIVIVPK.inputPK.set(protPK.outputFitting)
        protIVIVPK.inputLevy.set(protLevy.outputExperiment)
        self.launchProtocol(protIVIVPK)
        self.assertIsNotNone(protIVIVPK.outputExperiment.fnPKPD, "There was a problem with the simulation")
        self.validateFiles('ProtPKPDDissolutionPKSimulation', ProtPKPDDissolutionPKSimulation)

        # Internal validity
        print "Internal validity ..."
        protInternal = self.newProtocol(ProtPKPDIVIVCInternalValidity,
                                      objLabel='pkpd - internal validity'
                                      )
        protInternal.inputExperiment.set(protNCA.outputExperiment)
        protInternal.inputSimulated.set(protIVIVPK.outputExperiment)
        self.launchProtocol(protInternal)
        self.validateFiles('protInternal', ProtPKPDIVIVCInternalValidity)

        fnSummary=protInternal._getPath("summary.txt")
        self.assertTrue(os.path.exists(fnSummary), "There was a problem with the import")
        fh=open(fnSummary)
        lineNo=1
        for line in fh.readlines():
            mean=float(line.split('=')[-1])
            if lineNo==1:
                self.assertTrue(np.abs(mean)<1, "The AUC error is not correct")
            elif lineNo==2:
                self.assertTrue(np.abs(mean)<5, "The Cmax error is not correct")
            lineNo+=1

        fh.close()

if __name__ == "__main__":
    unittest.main()
