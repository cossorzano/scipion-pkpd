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

class TestGabrielssonPD11Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Gabrielsson_PD11')
        cls.exptFn = cls.dataset.getFile('experiment')

    def testGabrielssonPD11Workflow(self):
        # Import an experiment (intravenous)

        print "Import Experiment"
        protImport = self.newProtocol(ProtImportExperiment,
                                      objLabel='pkpd - import experiment',
                                      inputFile=self.exptFn)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputExperiment.fnPKPD, "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Fit a PD model
        print "Fitting a Gompertz E0 PD model ..."
        protPDfitting = self.newProtocol(ProtPKPDGenericFit,
                                         objLabel='pkpd - pd Gompertz e0',
                                         predictor='distance',predicted='waterContent',
                                         modelType=4, fitType=0,
                                         bounds='(0.0, 3.0); (10.0, 50.0); (0.0, 3.0); (0.0, 0.5)')
        protPDfitting.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPDfitting)
        self.assertIsNotNone(protPDfitting.outputExperiment.fnPKPD, "There was a problem with the mono-compartmental model ")
        self.assertIsNotNone(protPDfitting.outputFitting.fnFitting, "There was a problem with the mono-compartmental model ")
        self.validateFiles('protPDfitting', protPDfitting)
        experiment = PKPDExperiment()
        experiment.load(protPDfitting.outputExperiment.fnPKPD)
        e0 = float(experiment.samples['Population measures'].descriptors['e0'])
        a = float(experiment.samples['Population measures'].descriptors['a'])
        b = float(experiment.samples['Population measures'].descriptors['b'])
        g = float(experiment.samples['Population measures'].descriptors['g'])
        self.assertTrue(e0>1.5 and e0<1.8) # Gompertz model in Gabrielsson does not have e0, which explains
        self.assertTrue(a>19 and a<21) # Gabrielsson p 995: 22.5
        self.assertTrue(b>2.7 and b<2.9) # Gabrielsson p 995: 2.1
        self.assertTrue(g>0.4 and g<0.6) # Gabrielsson p 995: 0.388
        fitting = PKPDFitting()
        fitting.load(protPDfitting.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.97)

        # Fit a PD model
        print "Fitting a Gompertz PD model ..."
        protPDfitting = self.newProtocol(ProtPKPDGenericFit,
                                         objLabel='pkpd - pd Gompertz',
                                         predictor='distance',predicted='waterContent',
                                         modelType=4, fitType=0,
                                         bounds='(0.0, 0.0); (10.0, 50.0); (0.0, 3.0); (0.0, 0.5)')
        protPDfitting.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPDfitting)
        self.assertIsNotNone(protPDfitting.outputExperiment.fnPKPD, "There was a problem with the Gompertz model ")
        self.assertIsNotNone(protPDfitting.outputFitting.fnFitting, "There was a problem with the Gompertz model ")
        self.validateFiles('protPDfitting', protPDfitting)
        experiment = PKPDExperiment()
        experiment.load(protPDfitting.outputExperiment.fnPKPD)
        e0 = float(experiment.samples['Population measures'].descriptors['e0'])
        a = float(experiment.samples['Population measures'].descriptors['a'])
        b = float(experiment.samples['Population measures'].descriptors['b'])
        g = float(experiment.samples['Population measures'].descriptors['g'])
        self.assertTrue(e0==0)
        self.assertTrue(a>22 and a<23) # Gabrielsson p 995: 22.5
        self.assertTrue(b>2 and b<2.2) # Gabrielsson p 995: 2.1
        self.assertTrue(g>0.38 and g<0.4) # Gabrielsson p 995: 0.388
        fitting = PKPDFitting()
        fitting.load(protPDfitting.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.97)

        # Fit a PD model
        print "Fitting a Logistic1 PD model ..."
        protPDfitting = self.newProtocol(ProtPKPDGenericFit,
                                         objLabel='pkpd - pd Logistic1',
                                         predictor='distance',predicted='waterContent',
                                         modelType=5, fitType=0,
                                         bounds='(0.0, 0.0); (10.0, 30.0); (2.0, 6.0); (0.2, 2.0)')
        protPDfitting.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPDfitting)
        self.assertIsNotNone(protPDfitting.outputExperiment.fnPKPD, "There was a problem with the Logistic1 model ")
        self.assertIsNotNone(protPDfitting.outputFitting.fnFitting, "There was a problem with the Logistic1 model ")
        self.validateFiles('protPDfitting', protPDfitting)
        experiment = PKPDExperiment()
        experiment.load(protPDfitting.outputExperiment.fnPKPD)
        e0 = float(experiment.samples['Population measures'].descriptors['e0'])
        a = float(experiment.samples['Population measures'].descriptors['a'])
        b = float(experiment.samples['Population measures'].descriptors['b'])
        g = float(experiment.samples['Population measures'].descriptors['g'])
        self.assertTrue(e0==0)
        self.assertTrue(a>21.4 and a<21.6) # Gabrielsson p 995: 21.5
        self.assertTrue(b>3.9 and b<4) # Gabrielsson p 995: 3.95
        self.assertTrue(g>0.6 and g<0.65) # Gabrielsson p 995: 0.62
        fitting = PKPDFitting()
        fitting.load(protPDfitting.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Fit a PD model
        print "Fitting a Weibull PD model ..."
        protPDfitting = self.newProtocol(ProtPKPDGenericFit,
                                         objLabel='pkpd - pd Weibull',
                                         predictor='distance',predicted='waterContent',
                                         modelType=11, fitType=0,
                                         bounds='(0.0, 30.0); (10.0, 30.0); (0.0, 0.004); (0.0, 5.0)')
        protPDfitting.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPDfitting)
        self.assertIsNotNone(protPDfitting.outputExperiment.fnPKPD, "There was a problem with the Weibull model ")
        self.assertIsNotNone(protPDfitting.outputFitting.fnFitting, "There was a problem with the Weibull model ")
        self.validateFiles('protPDfitting', protPDfitting)
        experiment = PKPDExperiment()
        experiment.load(protPDfitting.outputExperiment.fnPKPD)
        a = float(experiment.samples['Population measures'].descriptors['a'])
        b = float(experiment.samples['Population measures'].descriptors['b'])
        g = float(experiment.samples['Population measures'].descriptors['g'])
        d = float(experiment.samples['Population measures'].descriptors['d'])
        self.assertTrue(a>21 and a<21.2) # Gabrielsson p 996: 21.1
        self.assertTrue(b>19.7 and b<19.9) # Gabrielsson p 996: 19.81
        self.assertTrue(g>0.001 and g<0.002) # Gabrielsson p 996: 0.0018
        self.assertTrue(d>3.16 and d<3.19) # Gabrielsson p 996: 3.169
        fitting = PKPDFitting()
        fitting.load(protPDfitting.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Fit a PD model
        print "Fitting a Richards PD model ..."
        protPDfitting = self.newProtocol(ProtPKPDGenericFit,
                                         objLabel='pkpd - pd Richards',
                                         predictor='distance',predicted='waterContent',
                                         modelType=9, fitType=0,
                                         bounds='(0.0, 0.0); (10.0, 30.0); (2.0, 7.0); (0.0, 1.0); (0.0, 3.0)')
        protPDfitting.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPDfitting)
        self.assertIsNotNone(protPDfitting.outputExperiment.fnPKPD, "There was a problem with the Richards model ")
        self.assertIsNotNone(protPDfitting.outputFitting.fnFitting, "There was a problem with the Richards model ")
        self.validateFiles('protPDfitting', protPDfitting)
        experiment = PKPDExperiment()
        experiment.load(protPDfitting.outputExperiment.fnPKPD)
        e0 = float(experiment.samples['Population measures'].descriptors['e0'])
        a = float(experiment.samples['Population measures'].descriptors['a'])
        b = float(experiment.samples['Population measures'].descriptors['b'])
        g = float(experiment.samples['Population measures'].descriptors['g'])
        d = float(experiment.samples['Population measures'].descriptors['d'])
        self.assertTrue(e0==0)
        self.assertTrue(a>21.1 and a<21.3) # Gabrielsson p 997: 21.2
        self.assertTrue(b>5.5 and b<5.9) # Gabrielsson p 997: 5.7
        self.assertTrue(g>0.75 and g<0.85) # Gabrielsson p 997: 0.778
        self.assertTrue(d>1.55 and d<1.7) # Gabrielsson p 997: 1.62
        fitting = PKPDFitting()
        fitting.load(protPDfitting.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Fit a PD model
        print "Fitting a Morgan PD model ..."
        protPDfitting = self.newProtocol(ProtPKPDGenericFit,
                                         objLabel='pkpd - pd Morgan',
                                         predictor='distance',predicted='waterContent',
                                         modelType=10, fitType=0,
                                         bounds='(0.0, 0.0); (1.0, 2.0); (3000.0, 7000.0); (10.0, 30.0); (2.0, 7.0)')
        protPDfitting.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPDfitting)
        self.assertIsNotNone(protPDfitting.outputExperiment.fnPKPD, "There was a problem with the Morgan model ")
        self.assertIsNotNone(protPDfitting.outputFitting.fnFitting, "There was a problem with the Morgan model ")
        self.validateFiles('protPDfitting', protPDfitting)
        experiment = PKPDExperiment()
        experiment.load(protPDfitting.outputExperiment.fnPKPD)
        e0 = float(experiment.samples['Population measures'].descriptors['e0'])
        a = float(experiment.samples['Population measures'].descriptors['a'])
        b = float(experiment.samples['Population measures'].descriptors['b'])
        g = float(experiment.samples['Population measures'].descriptors['g'])
        d = float(experiment.samples['Population measures'].descriptors['d'])
        self.assertTrue(e0==0)
        self.assertTrue(a>21.8 and a<22.3) # Gabrielsson p 998: 22.1
        self.assertTrue(b>1.6 and b<1.8) # Gabrielsson p 998: 1.64
        self.assertTrue(g>4000 and g<7000) # Gabrielsson p 998: 5303
        self.assertTrue(d>4.3 and d<4.7) # Gabrielsson p 998: 4.5
        fitting = PKPDFitting()
        fitting.load(protPDfitting.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

        # Fit a PD model
        print "Fitting a Hill PD model ..."
        protPDfitting = self.newProtocol(ProtPKPDGenericFit,
                                         objLabel='pkpd - pd Hill',
                                         predictor='distance',predicted='waterContent',
                                         modelType=12, fitType=0,
                                         bounds='(0.0, 3.0); (0.0, 30.0); (0.0, 10.0); (0.0, 10.0)')
        protPDfitting.inputExperiment.set(protImport.outputExperiment)
        self.launchProtocol(protPDfitting)
        self.assertIsNotNone(protPDfitting.outputExperiment.fnPKPD, "There was a problem with the Hill model ")
        self.assertIsNotNone(protPDfitting.outputFitting.fnFitting, "There was a problem with the Hill model ")
        self.validateFiles('protPDfitting', protPDfitting)
        experiment = PKPDExperiment()
        experiment.load(protPDfitting.outputExperiment.fnPKPD)
        e0 = float(experiment.samples['Population measures'].descriptors['e0'])
        b = float(experiment.samples['Population measures'].descriptors['b'])
        g = float(experiment.samples['Population measures'].descriptors['g'])
        d = float(experiment.samples['Population measures'].descriptors['d'])
        self.assertTrue(e0>1.6 and e0<1.7) # Gabrielsson p 999: alpha=1.65
        self.assertTrue(b>20.3 and b<20.5) # Gabrielsson p 999: 20.4
        self.assertTrue(g>6.5 and g<6.7) # Gabrielsson p 999: 6.63
        self.assertTrue(d>4.5 and d<4.7) # Gabrielsson p 999: 4.56
        fitting = PKPDFitting()
        fitting.load(protPDfitting.outputFitting.fnFitting)
        self.assertTrue(fitting.sampleFits[0].R2>0.99)

if __name__ == "__main__":
    unittest.main()
