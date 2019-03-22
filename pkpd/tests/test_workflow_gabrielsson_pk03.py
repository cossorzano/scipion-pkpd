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


class TestGabrielssonPK03Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PK03')
        cls.exptFn = cls.dataset.getFile('experiment')

    
    def testGabrielssonPK03Workflow(self):
        print "Import Experiment"
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

        # Change the concentration unit to mg/L
        print "Change Units"
        protChangeCpUnit = self.newProtocol(ProtPKPDChangeUnits,
                                            objLabel='pkpd - change units (Cp to mg/L)',
                                            labelToChange='Cp', newUnitsCategory=4, newUnitsCategoryConc=1)
        protChangeCpUnit.inputExperiment.set(protChangeTimeUnit.outputExperiment)
        self.launchProtocol(protChangeCpUnit)
        self.assertIsNotNone(protChangeCpUnit.outputExperiment.fnPKPD, "There was a problem with changing units")
        self.validateFiles('protChangeUnits', protChangeCpUnit)

        # Fit a monocompartmental model with zero order absorption
        print "Fitting monocompartmental model with zero order..."
        protEV0MonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                  objLabel='pkpd - ev0 monocompartment',
                                                  bounds='(5.0, 20.0); (0.0, 0.2); (0.4, 1.2); (50.0, 150.0)')
        protEV0MonoCompartment.inputExperiment.set(protChangeCpUnit.outputExperiment)
        self.launchProtocol(protEV0MonoCompartment)
        self.assertIsNotNone(protEV0MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV0MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV0MonoCompartment', protEV0MonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protEV0MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Rin = float(experiment.samples['Individual'].descriptors['Oral_Rin'])
        tlag = float(experiment.samples['Individual'].descriptors['Oral_tlag'])
        self.assertTrue(Cl>0.72 and Cl<0.73)
        self.assertTrue(V>95 and V<98)  # Gabrielsson, p 522, VF=96.2
        self.assertTrue(Rin>0.080 and Rin<0.085)
        self.assertTrue(tlag>15 and tlag<20)
        fitting = PKPDFitting()
        fitting.load(protEV0MonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.95)
        self.assertTrue(fitting.sampleFits[0].AIC<-20)

        # Change via to ev1
        print "Change via to ev1"
        protChangeVia1 = self.newProtocol(ProtPKPDChangeVia,
                                          objLabel='pkpd - change via ev1',
                                          viaName='Oral', newViaType="ev1", tlag="", bioavailability=1.0)
        protChangeVia1.inputExperiment.set(protChangeCpUnit.outputExperiment)
        self.launchProtocol(protChangeVia1)
        self.assertIsNotNone(protChangeVia1.outputExperiment.fnPKPD, "There was a problem with changing via")
        self.validateFiles('protChangeVia1', protChangeVia1)

        # Fit a monocompartmental model with 1st order absorption
        print "Fitting monocompartmental model with 1st order..."
        protEV1MonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                  objLabel='pkpd - ev1 monocompartment',
                                                  bounds='(10.0, 30.0); (0.0, 0.05); (0.4, 1.2); (50.0, 150.0)')
        protEV1MonoCompartment.inputExperiment.set(protChangeVia1.outputExperiment)
        self.launchProtocol(protEV1MonoCompartment)
        self.assertIsNotNone(protEV1MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV1MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV1MonoCompartment', protEV1MonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protEV1MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Ka = float(experiment.samples['Individual'].descriptors['Oral_Ka'])
        tlag = float(experiment.samples['Individual'].descriptors['Oral_tlag'])
        self.assertTrue(Cl>0.75 and Cl<0.8)
        self.assertTrue(V>91 and V<102) # Gabrielsson p 521, Solution I: FV=98.7 L
        self.assertTrue(Ka>0.0075 and Ka<0.0085) # Gabrielsson p 521, Solution I: ka=0.418 1/h=0.007 1/min
        self.assertTrue(tlag>20 and tlag<30) # Gabrielsson p 521, Solution I: tlag=0.39h=23.4 min
        fitting = PKPDFitting()
        fitting.load(protEV1MonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.92)
        self.assertTrue(fitting.sampleFits[0].AIC<-15)

        # Change via to ev01
        print "Change via to ev01"
        protChangeVia01 = self.newProtocol(ProtPKPDChangeVia,
                                           objLabel='pkpd - change via ev01',
                                           viaName='Oral', newViaType="ev01", tlag="", bioavailability=1.0)
        protChangeVia01.inputExperiment.set(protChangeCpUnit.outputExperiment)
        self.launchProtocol(protChangeVia01)
        self.assertIsNotNone(protChangeVia01.outputExperiment.fnPKPD, "There was a problem with changing via")
        self.validateFiles('protChangeVia01', protChangeVia01)

        # Fit a monocompartmental model with 0th and 1st order absorption
        print "Fitting monocompartmental model with 0th and 1st order..."
        protEV01MonoCompartment = self.newProtocol(ProtPKPDMonoCompartment,
                                                  objLabel='pkpd - ev01 monocompartment',
                                                  bounds='(10.0, 30.0); (0.04, 0.08); (220.0, 300.0); (0.0, 0.02); (0.7, 0.9); (20.0, 60.0)')
        protEV01MonoCompartment.inputExperiment.set(protChangeVia01.outputExperiment)
        self.launchProtocol(protEV01MonoCompartment)
        self.assertIsNotNone(protEV01MonoCompartment.outputExperiment.fnPKPD, "There was a problem with the monocompartmental model ")
        self.assertIsNotNone(protEV01MonoCompartment.outputFitting.fnFitting, "There was a problem with the monocompartmental model ")
        self.validateFiles('protEV01MonoCompartment', protEV01MonoCompartment)
        experiment = PKPDExperiment()
        experiment.load(protEV01MonoCompartment.outputExperiment.fnPKPD)
        Cl = float(experiment.samples['Individual'].descriptors['Cl'])
        V = float(experiment.samples['Individual'].descriptors['V'])
        Ka = float(experiment.samples['Individual'].descriptors['Oral_Ka'])
        tlag = float(experiment.samples['Individual'].descriptors['Oral_tlag'])
        Rin = float(experiment.samples['Individual'].descriptors['Oral_Rin'])
        t0 = float(experiment.samples['Individual'].descriptors['Oral_t0'])
        self.assertTrue(Cl>0.73 and Cl<0.76)
        self.assertTrue(V>32 and V<36)
        self.assertTrue(Ka>0.00045 and Ka<0.0055)
        self.assertTrue(tlag>20 and tlag<30)
        self.assertTrue(Rin>0.05 and Rin<0.06)
        self.assertTrue(t0>290 and t0<310)
        fitting = PKPDFitting()
        fitting.load(protEV01MonoCompartment.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.975)
        self.assertTrue(fitting.sampleFits[0].AIC<-25)


if __name__ == "__main__":
    unittest.main()
