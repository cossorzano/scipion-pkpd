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

import os

from pyworkflow.tests import *
from pkpd.protocols import *
from pkpd.objects import PKPDDataSet
from test_workflow import TestWorkflow

class TestLevyPlotWorkflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Dissolution')
        cls.fnInVitro = cls.dataset.getFile('experiment12Test')
        cls.fnInVivo = cls.dataset.getFile('invivo12')

    def testDissolutionWorkflow(self):
        print "Import Experiment in vitro"
        protImportInVitro = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import in vitro',
                                      inputFile=self.fnInVitro)
        self.launchProtocol(protImportInVitro)
        self.assertIsNotNone(protImportInVitro.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImportInVitro)

        # Mahalanobis
        print "Mahalanobis ..."
        protMah1= self.newProtocol(ProtPKPDStatsMahalanobis,
                                   objLabel='pkpd - Mahalanobis measurements',
                                   labels='C'
                                   )
        protMah1.inputExperiment1.set(protImportInVitro.outputExperiment)
        self.launchProtocol(protMah1)
        fnSummary = protMah1._getPath("report.txt")
        self.assertTrue(os.path.exists(fnSummary))

        # Average sample
        print "Average sample ..."
        protAvg= self.newProtocol(ProtPKPDAverageSample,
                                   objLabel='pkpd - Average sample'
                                   )
        protAvg.inputExperiment.set(protImportInVitro.outputExperiment)
        self.launchProtocol(protAvg)
        self.assertIsNotNone(protAvg.outputExperiment.fnPKPD, "There was a problem with the import")

        # Change the time unit to minute
        print "Change Units"
        protChangeTimeUnit = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (t to min)',
                                              labelToChange='t', newUnitsCategory=0, newUnitsCategoryTime=1)
        protChangeTimeUnit.inputExperiment.set(protImportInVitro.outputExperiment)
        self.launchProtocol(protChangeTimeUnit)
        self.assertIsNotNone(protChangeTimeUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeTimeUnit)

        print "Import Experiment in vivo"
        protImportInVivo = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import in vivo',
                                      inputFile=self.fnInVivo)
        self.launchProtocol(protImportInVivo)
        self.assertIsNotNone(protImportInVivo.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImportInVivo)

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
        self.assertTrue(AUC0t > 325.5 and AUC0t < 327.5)
        AUMC0t = float(experiment.samples['Individual1'].descriptors['AUMC0t'])
        self.assertTrue(AUMC0t > 42165 and AUMC0t < 42170)
        Cmax = float(experiment.samples['Individual1'].descriptors['Cmax'])
        self.assertTrue(Cmax > 1.9 and Cmax < 2.1)
        Tmax = float(experiment.samples['Individual1'].descriptors['Tmax'])
        self.assertTrue(Tmax > 39 and Tmax < 41)
        MRT = float(experiment.samples['Individual1'].descriptors['MRT'])
        self.assertTrue(MRT > 129 and MRT < 130)

        # Fit a Weibull dissolution
        print "Fitting Weibull model ..."
        protWeibull = self.newProtocol(ProtPKPDDissolutionFit,
                                objLabel='pkpd - fit dissolution Weibull',
                                globalSearch=True, modelType=3)
        protWeibull.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protWeibull)
        self.assertIsNotNone(protWeibull.outputExperiment.fnPKPD, "There was a problem with the dissolution model")
        self.assertIsNotNone(protWeibull.outputFitting.fnFitting, "There was a problem with the dissolution model")
        self.validateFiles('ProtPKPDDissolutionFit', ProtPKPDDissolutionFit)

        # Fit Order 1
        print "Fitting EV1-monocompartment model ..."
        protModelInVivo = self.newProtocol(ProtPKPDMonoCompartment,
                                       objLabel='pkpd - fit monocompartment',
                                       bounds="(0.0, 8.0); (0.0, 0.2); (0.0, 10.0); (10.0, 50.0)"
                                       )
        protModelInVivo.inputExperiment.set(protImportInVivo.outputExperiment)
        self.launchProtocol(protModelInVivo)
        self.assertIsNotNone(protModelInVivo.outputExperiment.fnPKPD, "There was a problem with the PK model")
        self.assertIsNotNone(protModelInVivo.outputFitting.fnFitting, "There was a problem with the PK model")
        self.validateFiles('ProtPKPDMonoCompartment', ProtPKPDMonoCompartment)

        # Deconvolve the in vivo
        print "Deconvolving in vivo ..."
        protDeconv = self.newProtocol(ProtPKPDDeconvolve,
                                       objLabel='pkpd - deconvolution'
                                       )
        protDeconv.inputODE.set(protModelInVivo)
        self.launchProtocol(protDeconv)
        self.assertIsNotNone(protDeconv.outputExperiment.fnPKPD, "There was a problem with the deconvolution")
        self.validateFiles('ProtPKPDDeconvolve', ProtPKPDDeconvolve)

        # Deconvolve the in vivo
        print "Deconvolving in vivo Wagner Nelson..."
        protDeconvWN = self.newProtocol(ProtPKPDDeconvolutionWagnerNelson,
                                        objLabel='pkpd - deconvolution Wagner Nelson'
                                       )
        protDeconvWN.inputExperiment.set(protModelInVivo.outputExperiment)
        self.launchProtocol(protDeconvWN)
        self.assertIsNotNone(protDeconvWN.outputExperiment.fnPKPD, "There was a problem with the deconvolution Wagner")
        self.validateFiles('ProtPKPDDeconvolutionWagnerNelson', ProtPKPDDeconvolutionWagnerNelson)

        # Fit bicompartment
        print "Fitting EV1-twocompartments model ..."
        protModelInVivo2 = self.newProtocol(ProtPKPDTwoCompartments,
                                       objLabel='pkpd - fit two compartments',
                                       bounds="(0.0, 8.0); (0.0, 0.2); (0.0, 10.0); (0.1, 50.0);  (0.0, 10.0); (0.1, 50.0)"
                                       )
        protModelInVivo2.inputExperiment.set(protImportInVivo.outputExperiment)
        self.launchProtocol(protModelInVivo2)
        self.assertIsNotNone(protModelInVivo2.outputExperiment.fnPKPD, "There was a problem with the PK model")
        self.assertIsNotNone(protModelInVivo2.outputFitting.fnFitting, "There was a problem with the PK model")
        self.validateFiles('ProtPKPDTwoCompartments', ProtPKPDTwoCompartments)

        # Deconvolve the in vivo
        print "Deconvolving in vivo Loo Riegelman ..."
        protDeconvLR = self.newProtocol(ProtPKPDDeconvolutionLooRiegelman,
                                        objLabel='pkpd - deconvolution Loo Riegelman'
                                       )
        protDeconvLR.inputExperiment.set(protModelInVivo2.outputExperiment)
        self.launchProtocol(protDeconvLR)
        self.assertIsNotNone(protDeconvLR.outputExperiment.fnPKPD, "There was a problem with the deconvolution Loo")
        self.validateFiles('ProtPKPDDeconvolutionLooRiegelman', ProtPKPDDeconvolutionLooRiegelman)

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
                                      timeScale='$[a0]+$[a1]*$(t)+$[a2]*np.power($(t),2.0)',
                                      timeBounds='a0: [-100,100]; a1: [1e-4,1e3]; a2: [1e-8,1e3]',
                                      responseScale='$[k1]*$(Adissol)',
                                      responseBounds='k1: [1e-3,1e3]',
                                      objLabel='pkpd - ivivc generic'
                                     )
        protIVIVCG.inputInVitro.set(protWeibull)
        protIVIVCG.inputInVivo.set(protDeconv)
        self.launchProtocol(protIVIVCG)
        self.assertIsNotNone(protIVIVCG.outputExperimentFabs.fnPKPD, "There was a problem with the IVIVC Generic")
        self.assertIsNotNone(protIVIVCG.outputExperimentAdissol.fnPKPD, "There was a problem with the IVIVC Generic")
        self.validateFiles('ProtPKPDDissolutionIVIVCG', ProtPKPDDissolutionIVIVCGeneric)

        # IVIVC Wagner
        print "In vitro-in vivo correlation Wagner Nelson..."
        protIVIVCWN = self.newProtocol(ProtPKPDDissolutionIVIVC,
                                       responseScale=1,
                                       objLabel='pkpd - ivivc'
                                      )
        protIVIVCWN.inputInVitro.set(protWeibull)
        protIVIVCWN.inputInVivo.set(protDeconvWN)
        self.launchProtocol(protIVIVCWN)
        self.assertIsNotNone(protIVIVCWN.outputExperimentFabs.fnPKPD, "There was a problem with the IVIVC")
        self.validateFiles('ProtPKPDDissolutionIVIVC', ProtPKPDDissolutionIVIVC)

        # Dissolution simulation
        print "IVIV+PK simulation ..."
        protIVIVPK = self.newProtocol(ProtPKPDDissolutionPKSimulation,
                                      objLabel='pkpd - ivivc+pk',
                                      inputN=100,
                                      tF=15,
                                      addIndividuals=True,
                                      inputDose=100
                                      )
        protIVIVPK.inputInVitro.set(protWeibull.outputFitting)
        protIVIVPK.inputPK.set(protModelInVivo.outputFitting)
        protIVIVPK.inputIvIvC.set(protIVIVC.outputExperimentFabs)
        self.launchProtocol(protIVIVPK)
        self.assertIsNotNone(protIVIVPK.outputExperiment.fnPKPD, "There was a problem with the simulation")
        self.validateFiles('ProtPKPDDissolutionPKSimulation', ProtPKPDDissolutionPKSimulation)

        # Internal validity
        print "Internal validity ..."
        protInternal = self.newProtocol(ProtPKPDIVIVCInternalValidity,
                                       objLabel='pkpd - internal validity',
                                       )
        protInternal.inputExperiment.set(protNCA.outputExperiment)
        protInternal.inputSimulated.set(protIVIVPK.outputExperiment)
        self.launchProtocol(protInternal)
        fnSummary = protInternal._getPath("summary.txt")
        self.assertTrue(os.path.exists(fnSummary))
        lineNo = 0
        for line in open(fnSummary).readlines():
            tokens = line.split('=')
            if lineNo==0:
                AUCmean=float(tokens[-1])
                self.assertTrue(AUCmean>5 and AUCmean<60)
            elif lineNo==1:
                Cmaxmean = float(tokens[-1])
                self.assertTrue(Cmaxmean > 40 and Cmaxmean < 80)
            lineNo+=1

        # Bootstrap dissolution
        print "Dissolution bootstrap ..."
        protDissolBootstrap = self.newProtocol(ProtPKPDFitBootstrap,
                                      objLabel='pkpd - dissol bootstrap',
                                      Nbootstrap=10
                                      )
        protDissolBootstrap.inputFit.set(protWeibull)
        self.launchProtocol(protDissolBootstrap)
        self.assertIsNotNone(protDissolBootstrap.outputPopulation.fnFitting, "There was a problem with the dissolution bootstrap")
        self.validateFiles('ProtPKPDFitBootstrap', ProtPKPDFitBootstrap)

        # Bootstrap ODE
        print "ODE bootstrap ..."
        protODEBootstrap = self.newProtocol(ProtPKPDODEBootstrap,
                                               objLabel='pkpd - ode bootstrap',
                                               Nbootstrap=3
                                               )
        protODEBootstrap.inputODE.set(protModelInVivo)
        self.launchProtocol(protODEBootstrap)
        self.assertIsNotNone(protODEBootstrap.outputPopulation.fnFitting,
                             "There was a problem with the ODE bootstrap")
        self.validateFiles('ProtPKPDODEBootstrap', ProtPKPDODEBootstrap)

        # Dissolution simulation
        print "IVIV+PK simulation ..."
        protIVIVPKBoot = self.newProtocol(ProtPKPDDissolutionPKSimulation,
                                      objLabel='pkpd - ivivc+pk bootstrap',
                                      inputN=100,
                                      tF=15,
                                      addIndividuals=True
                                      )
        protIVIVPKBoot.inputInVitro.set(protDissolBootstrap.outputPopulation)
        protIVIVPKBoot.inputPK.set(protODEBootstrap.outputPopulation)
        protIVIVPKBoot.inputIvIvC.set(protIVIVC.outputExperimentFabs)
        self.launchProtocol(protIVIVPKBoot)
        self.assertIsNotNone(protIVIVPKBoot.outputExperiment.fnPKPD, "There was a problem with the dissolution model")
        self.validateFiles('ProtPKPDDissolutionPKSimulation', ProtPKPDDissolutionPKSimulation)

        # t Test
        print "T-test ..."
        protTtest= self.newProtocol(ProtPKPDStatsExp2Subgroups2Mean,
                                    objLabel='pkpd - t test',
                                    label1 = 'MRT',
                                    label2 = 'MRT'
                                    )
        protTtest.inputExperiment1.set(protIVIVPKBoot.outputExperiment)
        protTtest.inputExperiment2.set(protIVIVPK.outputExperiment)
        self.launchProtocol(protTtest)
        fnSummary = protTtest._getPath("report.txt")
        self.assertTrue(os.path.exists(fnSummary))
        # for line in open(fnSummary).readlines(): # In some cases, 5% it fails
        #     if '-statistic' in line:
        #         tokens = line.split('=')
        #         pval=float(tokens[-1])
        #         self.assertTrue(pval>0.1)

        # Kolmogorov test
        print "Kolmogorov test ..."
        protKtest= self.newProtocol(ProtPKPDStatsExp2Subgroups2Kolmogorov,
                                    objLabel='pkpd - Kolmogorov test',
                                    label1='MRT',
                                    label2='MRT'
                                    )
        protKtest.inputExperiment1.set(protIVIVPKBoot.outputExperiment)
        protKtest.inputExperiment2.set(protIVIVPK.outputExperiment)
        self.launchProtocol(protKtest)
        fnSummary = protKtest._getPath("report.txt")
        self.assertTrue(os.path.exists(fnSummary))
        # for line in open(fnSummary).readlines(): # In some cases, 5% it fails
        #     if '-statistic' in line:
        #         tokens = line.split('=')
        #         pval=float(tokens[-1])
        #         self.assertTrue(pval>0.1)


        # Mahalanobis
        print "Mahalanobis ..."
        protMah2 = self.newProtocol(ProtPKPDStatsMahalanobis,
                                    objLabel='pkpd - Mahalanobis labels',
                                    labels='AUC0t AUMC0t Cmax'
                                    )
        protMah2.inputExperiment1.set(protIVIVPKBoot.outputExperiment)
        protMah2.inputExperiment2.set(protIVIVPK.outputExperiment)
        self.launchProtocol(protMah2)
        fnSummary = protMah2._getPath("report.txt")
        self.assertTrue(os.path.exists(fnSummary))

if __name__ == "__main__":
    unittest.main()
