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
# *  e-mail address 'info@kinestat.com'
# *
# **************************************************************************

import numpy as np

import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from .protocol_pkpd_fit_base import ProtPKPDFitBase
from pkpd.objects import PKPDFitting, PKPDSampleFitBootstrap, PKPDLSOptimizer, PKPDModelBase2
from .protocol_pkpd_pdgeneric_fit import ProtPKPDGenericFit
from .protocol_pkpd_dissolution_fit import ProtPKPDDissolutionFit

# Tested by test_workflow_dissolution

class ProtPKPDFitBootstrap(ProtPKPDFitBase):
    """ Bootstrap estimate of generic fit models.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'fit bootstrap'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputFit', params.PointerParam, label="Input Fit model",
                      pointerClass='ProtPKPDGenericFit, ProtPKPDDissolutionFit',
                      help='Select a run of a fitted model')
        form.addParam('Nbootstrap', params.IntParam, label="Bootstrap samples", default=200, expertLevel=LEVEL_ADVANCED,
                      help='Number of bootstrap realizations for each sample')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runFit',self.inputFit.get().getObjId(), self.Nbootstrap.get(), self.confidenceInterval.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runFit(self, objId, Nbootstrap, confidenceInterval):
        self.protFit = self.inputFit.get()
        self.experiment = self.readExperiment(self.protFit.outputExperiment.fnPKPD)
        self.fitting = self.readFitting(self.protFit.outputFitting.fnFitting)

        # Get the X and Y variable names
        self.varNameX = self.fitting.predictor.varName
        self.varNameY = self.fitting.predicted.varName
        self.protFit.experiment = self.experiment
        self.protFit.varNameX = self.varNameX
        self.protFit.varNameY = self.varNameY

        # Setup model
        self.printSection("Model setup")
        self.protFit.model = self.protFit.createModel()
        self.protFit.model.setExperiment(self.experiment)
        self.protFit.model.setXVar(self.varNameX)
        self.protFit.model.setYVar(self.varNameY)
        self.protFit.setupFromFormParameters()
        self.protFit.setupModel()
        self.protFit.model.setBounds(self.protFit.bounds.get())
        self.protFit.model.printSetup()
        self.model = self.protFit.model

        # Setup self as model
        self.boundsList = self.model.bounds

        # Create output object
        self.fitting = PKPDFitting("PKPDSampleFitBootstrap")
        self.fitting.fnExperiment.set(self.experiment.fnPKPD.get())
        self.fitting.predictor=self.experiment.variables[self.varNameX]
        self.fitting.predicted=self.experiment.variables[self.varNameY]
        self.fitting.modelParameterUnits = None

        # Actual fitting
        if self.protFit.fitType.get()==0:
            fitType = "linear"
        elif self.protFit.fitType.get()==1:
            fitType = "log"
        elif self.protFit.fitType.get()==2:
            fitType = "relative"

        parameterNames = self.model.getParameterNames()
        for sampleName, sample in self.experiment.samples.iteritems():
            self.printSection("Fitting "+sampleName)

            x, y = sample.getXYValues(self.varNameX,self.varNameY)
            print("X= "+str(x))
            print("Y= "+str(y))
            print(" ")

            self.model.setXYValues(x, y)
            self.prepareForSampleAnalysis(sampleName)
            self.model.calculateParameterUnits(sample)
            if self.fitting.modelParameterUnits==None:
                self.fitting.modelParameterUnits = self.model.parameterUnits
            self.model.prepare()
            if self.model.bounds == None:
                continue
            print(" ")

            # Get the initial parameters
            parameters0 = []
            for parameterName in parameterNames:
                parameters0.append(float(sample.descriptors[parameterName]))
            print("Initial solution: %s"%str(parameters0))
            print(" ")

            # Output object
            sampleFit = PKPDSampleFitBootstrap()
            sampleFit.sampleName = sample.sampleName
            sampleFit.parameters = np.zeros((self.Nbootstrap.get(),len(parameters0)),np.double)
            sampleFit.xB = []
            sampleFit.yB = []

            # Bootstrap samples
            firstX = x[0]  # From [array(...)] to array(...)
            firstY = y[0]  # From [array(...)] to array(...)
            idx = [k for k in range(0,len(firstX))]
            for n in range(0,self.Nbootstrap.get()):
                ok = False
                while not ok:
                    lenToUse = len(idx)
                    idxB = sorted(np.random.choice(idx,lenToUse))
                    xB = [np.asarray([firstX[i] for i in idxB])]
                    yB = [np.asarray([firstY[i] for i in idxB])]

                    print("Bootstrap sample %d"%n)
                    print("X= "+str(xB))
                    print("Y= "+str(yB))
                    self.model.setXYValues(xB, yB)
                    self.model.parameters = parameters0

                    optimizer2 = PKPDLSOptimizer(self.model,fitType)
                    optimizer2.verbose = 0
                    try:
                        optimizer2.optimize()
                        ok=True
                    except Exception as e:
                        print(e)
                        raise(e)
                        ok=False

                # Evaluate the quality on the whole data set
                self.model.setXYValues(x, y)
                optimizer2.evaluateQuality()
                print(optimizer2.optimum)
                print("   R2 = %f R2Adj=%f AIC=%f AICc=%f BIC=%f"%(optimizer2.R2,optimizer2.R2adj,optimizer2.AIC,\
                                                                   optimizer2.AICc,optimizer2.BIC))

                # Keep this result
                sampleFit.parameters[n,:] = optimizer2.optimum
                sampleFit.xB.append(str(xB[0]))
                sampleFit.yB.append(str(yB[0]))
                sampleFit.copyFromOptimizer(optimizer2)

            self.fitting.sampleFits.append(sampleFit)

        self.fitting.modelParameters = self.getParameterNames()
        self.fitting.modelDescription = self.model.getDescription()
        self.fitting.write(self._getPath("bootstrapPopulation.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputPopulation=self.fitting)
        self._defineSourceRelation(self.inputFit.get(), self.fitting)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        msg.append("Number of bootstrap realizations: %d"%self.Nbootstrap.get())
        msg.append("Confidence interval: %f"%self.confidenceInterval.get())
        return msg

    def _validate(self):
        return []
