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


class TestGabrielssonPK06Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK06')
        cls.exptIVFn = cls.dataset.getFile('experiment_iv')
        cls.exptEVFn = cls.dataset.getFile('experiment_ev')
    
    def testGabrielssonPK06Workflow(self):
        # Import an experiment (intravenous)

        print "Import Experiment (intravenous doses)"
        protImportIV = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment iv',
                                      inputFile=self.exptIVFn)
        self.launchProtocol(protImportIV)
        self.assertIsNotNone(protImportIV.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImportIV', protImportIV)

        # Change the time unit to minute
        print "Change Units"
        protChangeTimeUnit = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (t to min)',
                                              labelToChange='t', newUnitsCategory=0, newUnitsCategoryTime=1)
        protChangeTimeUnit.inputExperiment.set(protImportIV.outputExperiment)
        self.launchProtocol(protChangeTimeUnit)
        self.assertIsNotNone(protChangeTimeUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeTimeUnit)

        # Fit a monocompartmental model to a set of measurements obtained by intravenous doses
        print "Fitting a monocompartmental model (intravenous)..."
        protIVMonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                  objLabel='pkpd - iv monocompartment',
                                                  bounds='(0.0, 0.1); (0.0, 300.0)')
        protIVMonoCompartment.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protIVMonoCompartment)
        self.assertIsNotNone(protIVMonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.validateFiles('protIVMonoCompartment', protIVMonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protIVMonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        self.assertTrue(Cl>0.095 and Cl<0.097) # Gabrielsson, p 548, Cle=6.0257 1/h=0.1004 1/min
        self.assertTrue(V>274 and V<275) # Gabrielsson, p 548, Vd=290.34
        fitting = PKPDFitting()
        fitting.load(protIVMonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-45)

        # Fit a monocompartmental model to a set of measurements obtained by intravenous doses and urine
        print "Fitting a monocompartmental model (intravenous doses and urine) ..."
        protIVMonoCompartmentUrine = self.newProtocol(ProtPKPDMonoCompartmentUrine,
                                                      objLabel='pkpd - iv monocompartment urine',
                                                      globalSearch=False,
                                                      bounds='(0.0,0.2);(150, 400.0);(0.0,1.0)')
        protIVMonoCompartmentUrine.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protIVMonoCompartmentUrine)
        self.assertIsNotNone(protIVMonoCompartmentUrine.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protIVMonoCompartmentUrine.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protIVMonoCompartmentUrine', protIVMonoCompartmentUrine)
        experiment = PKPDExperiment()
        experiment.load(protIVMonoCompartmentUrine.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        fe = float(experiment.samples['Individual'].descriptors['fe'])
        self.assertTrue(Cl>0.095 and Cl<0.097) # Gabrielsson, p 548, Cle=6.0257 1/h=0.1004 1/min
        self.assertTrue(V>274 and V<275) # Gabrielsson, p 548, Vd=290.34
        self.assertTrue(fe>0.068 and fe<0.070)
        fitting = PKPDFitting()
        fitting.load(protIVMonoCompartmentUrine.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-55)

        # Import an experiment (extravascular)
        print "Import Experiment (extravascular)"
        protImportEV = self.newProtocol(ProtImportExperiment,
                                        objLabel='pkpd - import experiment po',
                                        inputFile=self.exptEVFn)
        self.launchProtocol(protImportEV)
        self.assertIsNotNone(protImportEV.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImportEV', protImportEV)

        # Change the time unit to minute
        print "Change Units"
        protChangeTimeUnit2 = self.newProtocol(ProtPKPDChangeUnits,
                                              objLabel='pkpd - change units (t to min)',
                                              labelToChange='t', newUnitsCategory=0, newUnitsCategoryTime=1)
        protChangeTimeUnit2.inputExperiment.set(protImportEV.outputExperiment)
        self.launchProtocol(protChangeTimeUnit2)
        self.assertIsNotNone(protChangeTimeUnit2.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeTimeUnit2', protChangeTimeUnit2)

        # Fit a monocompartmental model to a set of measurements obtained by extravascular doses and urine
        print "Fitting a monocompartmental model (extravascular and urine)..."
        protEV1MonoCompartmentUrine = self.newProtocol(ProtPKPDMonoCompartmentUrine,
                                                       objLabel='pkpd - ev1 monocompartment urine',
                                                       globalSearch=False,
                                                       bounds='(0.0, 30.0); (0.0, 0.05); (0.0, 0.2); (250.0, 350.0); (0.0, 1.0)')
        protEV1MonoCompartmentUrine.inputExperiment.set(protChangeTimeUnit2.outputExperiment)
        self.launchProtocol(protEV1MonoCompartmentUrine)
        self.assertIsNotNone(protEV1MonoCompartmentUrine.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartmentUrine.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartmentUrine', protEV1MonoCompartmentUrine)
        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartmentUrine.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        fe = float(experiment.samples['Individual'].descriptors['fe'])
        Ka = float(experiment.samples['Individual'].descriptors['Oral_Ka'])
        tlag = float(experiment.samples['Individual'].descriptors['Oral_tlag'])
        self.assertTrue(Cl>0.084 and Cl<0.087) # Gabrielsson, p 548, Cle=6.0257 1/h=0.1004 1/min
        self.assertTrue(V>250 and V<275) # Gabrielsson, p 548, Vd=290.34
        self.assertTrue(fe>0.08 and fe<0.09) # Gabrielsson, p 548, fe=0.0698
        self.assertTrue(Ka>0.004 and Ka<0.006) # Gabrielsson, p 548, Ka=0.4207 1/h=0.00701 1/min
        self.assertTrue(tlag>13 and tlag<17) # Gabrielsson, p 548, tlag=0.3129 h=18.77 min
        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartmentUrine.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.965)
        self.assertTrue(fitting.sampleFits[0].AIC<-7)

        # Elimination rate of IV
        print "Elimination rate IV ..."
        protEliminationRate = self.newProtocol(ProtPKPDEliminationRate,
                                               objLabel='pkpd - elimination rate',
                                               predictor='t', predicted='Cp')
        protEliminationRate.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protEliminationRate)
        self.assertIsNotNone(protEliminationRate.outputExperiment.fnPKPD, "There was a problem with the exponential fitting")
        self.assertIsNotNone(protEliminationRate.outputFitting.fnFitting, "There was a problem with the exponential fitting")
        self.validateFiles('protEliminationRate', protEliminationRate)
        experiment = PKPDExperiment()
        experiment.load(protEliminationRate.outputExperiment.fnPKPD)
        c1 = float(experiment.samples['Individual'].descriptors['c1'])
        lambda1 = float(experiment.samples['Individual'].descriptors['lambda1'])
        self.assertAlmostEqual(c1,0.0455,2)
        self.assertAlmostEqual(lambda1,0.000351,3)

        fitting = PKPDFitting()
        fitting.load(protEliminationRate.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-45)

        # Non-compartmental analysis IV
        print "Performing Non-compartmental analysis of IV ..."
        protNCAIVObs = self.newProtocol(ProtPKPDNCAIVObs,
                                        objLabel='pkpd - nca iv observations')
        protNCAIVObs.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        protNCAIVObs.protElimination.set(protEliminationRate)
        self.launchProtocol(protNCAIVObs)
        self.assertIsNotNone(protNCAIVObs.outputExperiment.fnPKPD, "There was a problem with the Non-compartmental analysis ")
        self.assertIsNotNone(protNCAIVObs.outputAnalysis.fnAnalysis, "There was a problem with the Non-compartmental analysis ")
        self.validateFiles('protNCAIVObs', protNCAIVObs)

        experiment = PKPDExperiment()
        experiment.load(protNCAIVObs.outputExperiment.fnPKPD)
        AUC_0inf = float(experiment.samples['Individual'].descriptors['AUC_0inf'])
        self.assertTrue(AUC_0inf>120 and AUC_0inf<125)

        # Filter time variable in PO
        print "Filter time in PO"
        protFilterTime = self.newProtocol(ProtPKPDFilterMeasurements,
                                                  objLabel='pkpd - filter measurements t>800',
                                                  filterType=1, condition='$(t)>800')
        protFilterTime.inputExperiment.set(protChangeTimeUnit2.outputExperiment)
        self.launchProtocol(protFilterTime)
        self.assertIsNotNone(protFilterTime.outputExperiment.fnPKPD, "There was a problem with the filter")
        self.validateFiles('protFilterTime', protFilterTime)

        # Elimination rate of PO
        print "Elimination rate PO ..."
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
        self.assertAlmostEqual(c1,0.0780,2)
        self.assertAlmostEqual(lambda1,0.000284,3)

        fitting = PKPDFitting()
        fitting.load(protEliminationRate.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)
        self.assertTrue(fitting.sampleFits[0].AIC<-25)

        # Estimate absorption rate
        print "Estimation of the absorption rate..."
        protAbsorptionRate = self.newProtocol(ProtPKPDAbsorptionRate,
                                              objLabel='pkpd - absorption rate',
                                              bounds='(0,0.07);(200,500);(5,30)')
        protAbsorptionRate.inputExperiment.set(protChangeTimeUnit2.outputExperiment)
        protAbsorptionRate.protElimination.set(protEliminationRate)
        self.launchProtocol(protAbsorptionRate)
        self.assertIsNotNone(protAbsorptionRate.outputExperiment.fnPKPD, "There was a problem with the absorption rate estimation ")
        self.assertIsNotNone(protAbsorptionRate.outputFitting.fnFitting, "There was a problem with the absorption rate estimation ")
        self.validateFiles('protAbsorptionRate', protAbsorptionRate)

        experiment = PKPDExperiment()
        experiment.load(protAbsorptionRate.outputExperiment.fnPKPD)
        Ka = float(experiment.samples['Individual'].descriptors['Ka'])
        Ke = float(experiment.samples['Individual'].descriptors['Ke'])
        Vd = float(experiment.samples['Individual'].descriptors['Vd'])
        tlag = float(experiment.samples['Individual'].descriptors['tlag'])
        self.assertTrue(Ka>0.007 and Ka<0.045)
        self.assertTrue(Ke>0.00025 and Ke<0.00031)
        self.assertTrue(Vd>310 and Vd<470)
        self.assertTrue(tlag>15 and tlag<25)

        fitting = PKPDFitting()
        fitting.load(protAbsorptionRate.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.975)
        self.assertTrue(fitting.sampleFits[0].AIC<-27)

        # Non-compartmental analysis EV
        print "Performing Non-compartmental analysis of PO ..."
        protNCAEVObs = self.newProtocol(ProtPKPDNCAEV,
                                        objLabel='pkpd - nca ev')
        protNCAEVObs.protAbsorption.set(protAbsorptionRate)
        self.launchProtocol(protNCAEVObs)
        self.assertIsNotNone(protNCAEVObs.outputExperiment.fnPKPD, "There was a problem with the Non-compartmental analysis ")
        self.assertIsNotNone(protNCAEVObs.outputAnalysis.fnAnalysis, "There was a problem with the Non-compartmental analysis ")
        self.validateFiles('protNCAEVObs', protNCAEVObs)

        experiment = PKPDExperiment()
        experiment.load(protNCAEVObs.outputExperiment.fnPKPD)
        AUC_0inf = float(experiment.samples['Individual'].descriptors['AUC_0inf'])
        self.assertTrue(AUC_0inf>270 and AUC_0inf<273)

        # Non-compartmental analysis EV
        print "Estimating bioavailability ..."
        protBioavail = self.newProtocol(ProtPKPDNCAEstimateBioavailability,
                                        objLabel='pkpd - bioavailability nca')
        protBioavail.protNCAEV.set(protNCAEVObs)
        protBioavail.protNCAIV.set(protNCAIVObs)
        self.launchProtocol(protBioavail)
        self.assertIsNotNone(protBioavail.outputExperiment.fnPKPD, "There was a problem with the Non-compartmental analysis ")
        self.validateFiles('protBioavail', protBioavail)

        experiment = PKPDExperiment()
        experiment.load(protBioavail.outputExperiment.fnPKPD)
        bioavailability = float(experiment.samples['Individual'].descriptors['bioavailability'])
        self.assertTrue(bioavailability>1.09 and bioavailability<1.11) # Gabrielsson, p 546: F=110%

if __name__ == "__main__":
    unittest.main()
