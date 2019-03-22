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
from test_workflow import TestWorkflow
from pkpd.objects import PKPDDataSet


class TestGabrielssonPK07Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK07')
        cls.exptFn = cls.dataset.getFile('experiment')
    
    def testGabrielssonPK07Workflow(self):
        # Import an experiment

        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Fit a two-compartmentx model with intravenous absorption to a set of measurements
        print "Fitting a two-compartmentx model (intravenous)..."
        protPKPDIVTwoCompartments = self.newProtocol(ProtPKPDTwoCompartments,
                                                     objLabel='pkpd - iv two-compartments',
                                                     globalSearch=False,
                                                     bounds='(0.0, 1.0); (0.0, 100.0); (0.0, 2.0); (0.0, 100.0)')
        protPKPDIVTwoCompartments.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPKPDIVTwoCompartments)
        self.assertIsNotNone(protPKPDIVTwoCompartments.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protPKPDIVTwoCompartments.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protPKPDIVTwoCompartments', protPKPDIVTwoCompartments)
        experiment = PKPDExperiment()
        experiment.load(protPKPDIVTwoCompartments.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        self.assertTrue(Cl>0.35 and Cl<0.38)
        self.assertTrue(Clp>1 and Clp<1.04)
        self.assertTrue(V>55 and V<60)
        self.assertTrue(Vp>55 and Vp<60)
        fitting = PKPDFitting()
        fitting.load(protPKPDIVTwoCompartments.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-65)

        # Fit a two-compartmentx model with intravenous absorption to a set of measurements
        print "Fitting a two-compartmentx model (intravenous)..."
        protODERefine = self.newProtocol(ProtPKPDODERefine,
                                                     objLabel='pkpd - ode refinement')
        protODERefine.inputODE.set(protPKPDIVTwoCompartments)
        self.launchProtocol(protODERefine)
        self.assertIsNotNone(protODERefine.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protODERefine.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protODERefine', protODERefine)
        experiment = PKPDExperiment()
        experiment.load(protODERefine.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        Clp = float(experiment.samples['Individual'].descriptors['Clp'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Vp = float(experiment.samples['Individual'].descriptors['Vp'])
        self.assertTrue(Cl>0.35 and Cl<0.38)
        self.assertTrue(Clp>1 and Clp<1.04)
        self.assertTrue(V>55 and V<60)
        self.assertTrue(Vp>55 and Vp<60)
        fitting = PKPDFitting()
        fitting.load(protODERefine.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-65)

        # Fit a two-compartmentx model with intravenous absorption to a set of measurements
        print "Fitting two exponentials ..."
        protExponentials = self.newProtocol(ProtPKPDExponentialFit,
                                            objLabel='pkpd - exponential fit',
                                            Nexp=2)
        protExponentials.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protExponentials)
        self.assertIsNotNone(protExponentials.outputExperiment.fnPKPD, "There was a problem with the two-compartmental model ")
        self.assertIsNotNone(protExponentials.outputFitting.fnFitting, "There was a problem with the two-compartmental model ")
        self.validateFiles('protExponentials', protExponentials)
        experiment = PKPDExperiment()
        experiment.load(protExponentials.outputExperiment.fnPKPD)
        c1 = float(experiment.samples['Individual'].descriptors['c1'])
        lambda1 = float(experiment.samples['Individual'].descriptors['lambda1'])
        c2 = float(experiment.samples['Individual'].descriptors['c2'])
        lambda2 = float(experiment.samples['Individual'].descriptors['lambda2'])
        self.assertTrue(c1>1.03 and c1<1.05) # Gabrielsson p. 555 A=1.06
        self.assertTrue(lambda1>0.038 and lambda1<0.040) # Gabrielsson p. 555 alpha=0.048
        self.assertTrue(c2>0.72 and c2<0.73) # Gabrielsson p. 555 B=0.79
        self.assertTrue(lambda2>0.0028 and lambda2<0.0030) # Gabrielsson p. 555 beta=0.0031
        fitting = PKPDFitting()
        fitting.load(protExponentials.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-65)

        # Filter time variable
        print "Filter time"
        protFilterTime = self.newProtocol(ProtPKPDFilterMeasurements,
                                                  objLabel='pkpd - filter measurements t>116',
                                                  filterType=1, condition='$(t)>116')
        protFilterTime.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protFilterTime)
        self.assertIsNotNone(protFilterTime.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('protFilterTime', protFilterTime)

        # Elimination rate of PO
        print "Elimination rate ..."
        protEliminationRate = self.newProtocol(ProtPKPDEliminationRate,
                                               objLabel='pkpd - elimination rate',
                                               predictor='t', predicted='Cp')
        protEliminationRate.inputExperiment.set(protFilterTime.outputExperiment)
        self.launchProtocol(protEliminationRate)
        self.assertIsNotNone(protEliminationRate.outputExperiment.fnPKPD, "There was a problem with the exponential fitting")
        self.assertIsNotNone(protEliminationRate.outputFitting.fnFitting, "There was a problem with the exponential fitting")
        self.validateFiles('protEliminationRate', protEliminationRate)
        experiment = PKPDExperiment()
        experiment.load(protEliminationRate.outputExperiment.fnPKPD)
        c1 = float(experiment.samples['Individual'].descriptors['c1'])
        lambda1 = float(experiment.samples['Individual'].descriptors['lambda1'])
        self.assertTrue(c1>0.71 and c1<0.74) # Gabrielsson p. 555 B=0.79
        self.assertTrue(lambda1>0.0028 and lambda1<0.0030) # Gabrielsson p. 555 beta=0.0031

        fitting = PKPDFitting()
        fitting.load(protEliminationRate.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.96)
        self.assertTrue(fitting.sampleFits[0].AIC<-20)

        # Non-compartmental analysis IV
        print "Performing Non-compartmental analysis of IV ..."
        protNCAIVObs = self.newProtocol(ProtPKPDNCAIVObs,
                                        objLabel='pkpd - nca iv observations')
        protNCAIVObs.inputExperiment.set(protImport.outputExperiment)
        protNCAIVObs.protElimination.set(protEliminationRate)
        self.launchProtocol(protNCAIVObs)
        self.assertIsNotNone(protNCAIVObs.outputExperiment.fnPKPD, "There was a problem with the Non-compartmental analysis ")
        self.assertIsNotNone(protNCAIVObs.outputAnalysis.fnAnalysis, "There was a problem with the Non-compartmental analysis ")
        self.validateFiles('protNCAIVObs', protNCAIVObs)

        experiment = PKPDExperiment()
        experiment.load(protNCAIVObs.outputExperiment.fnPKPD)
        AUC_0inf = float(experiment.samples['Individual'].descriptors['AUC_0inf'])
        CL_0inf = float(experiment.samples['Individual'].descriptors['CL_0inf'])
        MRT = float(experiment.samples['Individual'].descriptors['MRT'])
        Vss = float(experiment.samples['Individual'].descriptors['Vss'])
        thalf = float(experiment.samples['Individual'].descriptors['thalf'])
        self.assertTrue(AUC_0inf>260 and AUC_0inf<265) # Gabrielsson p.552 AUC_0inf=250
        self.assertTrue(CL_0inf>0.37 and CL_0inf<0.39) # Gabrielsson p.552 Cl=0.4
        self.assertTrue(MRT>315 and MRT<317) # Gabrielsson p.552 MRT=280
        self.assertTrue(Vss>118 and Vss<122) # Gabrielsson p.552 Vss=110
        self.assertTrue(thalf>233 and thalf<237)

        # Non-compartmental analysis IV with exponentials
        print "Performing Non-compartmental analysis of IV with exponentials ..."
        protNCAIVExp = self.newProtocol(ProtPKPDNCAIVExp,
                                        objLabel='pkpd - nca iv exponentials')
        protNCAIVExp.protExponential.set(protExponentials)
        protNCAIVExp.protElimination.set(protEliminationRate)
        self.launchProtocol(protNCAIVExp)
        self.assertIsNotNone(protNCAIVExp.outputExperiment.fnPKPD, "There was a problem with the Non-compartmental analysis ")
        self.assertIsNotNone(protNCAIVExp.outputAnalysis.fnAnalysis, "There was a problem with the Non-compartmental analysis ")
        self.validateFiles('protNCAIVExp', protNCAIVExp)

        experiment = PKPDExperiment()
        experiment.load(protNCAIVExp.outputExperiment.fnPKPD)
        AUC_0inf = float(experiment.samples['Individual'].descriptors['AUC_0inf'])
        CL = float(experiment.samples['Individual'].descriptors['CL'])
        MRT = float(experiment.samples['Individual'].descriptors['MRT'])
        Vss = float(experiment.samples['Individual'].descriptors['Vss'])
        thalf = float(experiment.samples['Individual'].descriptors['thalf'])
        self.assertTrue(AUC_0inf>270 and AUC_0inf<275) # Gabrielsson p.552 AUC_0inf=250
        self.assertTrue(CL>0.36 and CL<0.37) # Gabrielsson p.552 Cl=0.4
        self.assertTrue(MRT>308 and MRT<311) # Gabrielsson p.552 MRT=280
        self.assertTrue(Vss>112 and Vss<115) # Gabrielsson p.552 Vss=110
        self.assertTrue(thalf>233 and thalf<237)

if __name__ == "__main__":
    unittest.main()
