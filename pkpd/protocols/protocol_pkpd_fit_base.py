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

import math
from collections import OrderedDict
from itertools import izip
import numpy as np

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import (PKPDDEOptimizer, PKPDLSOptimizer, PKPDFitting,
                          PKPDSampleFit)
from pkpd.utils import parseRange

class ProtPKPDFitBase(ProtPKPD):
    """ Base fit protocol"""

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams1(self, form, defaultPredictor, defaultPredicted):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('predictor', params.StringParam, label="Predictor variable (X)", default=defaultPredictor,
                      help='Y is predicted as an exponential function of X, Y=f(X)')
        form.addParam('predicted', params.StringParam, label="Predicted variable (Y)", default=defaultPredicted,
                      help='Y is predicted as an exponential function of X, Y=f(X)')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runFit',self.getInputExperiment().getObjId(),self.getListOfFormDependencies())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def getInputExperiment(self):
        if hasattr(self,"inputExperiment"):
            return self.inputExperiment.get()
        else:
            return None

    def getXYvars(self):
        if hasattr(self,"predictor"):
            self.varNameX=self.predictor.get()
        else:
            self.varNameX=None

        if hasattr(self,"predicted"):
            self.varNameY=self.predicted.get()
        else:
            self.varNameY=None

    def setExperiment(self, experiment):
        self.experiment = experiment
        if experiment!=None:
            self.fnExperiment = experiment.fnPKPD

    def setSample(self, sample):
        self.sample = sample
        self.model.setSample(sample)

    def createModel(self):
        pass

    def parseBounds(self, boundsString):
        self.boundsList = []

        if boundsString!="" and boundsString!=None:
            tokens = boundsString.split(';')
            if len(tokens)!=self.getNumberOfParameters():
                raise Exception("The number of bound intervals does not match the number of parameters")

            for token in tokens:
                values = token.strip().split(',')
                self.boundsList.append((float(values[0][1:]),float(values[1][:-1])))

    def getBounds(self):
        return self.boundsList

    def getParameterBounds(self):
        """ Return a dictionary where the parameter name is the key
        and the bounds are its values. """
        boundsDict = OrderedDict()
        self.parseBounds(self.bounds.get()) # after this we have boundsList
        parameterNames = self.model.getParameterNames()

        for paramName, bound in izip(parameterNames, self.getBounds()):
            boundsDict[paramName] = bound

        # Set None as bound for parameters not matched
        for paramName in parameterNames:
            if paramName not in boundsDict:
                boundsDict[paramName] = None

        return boundsDict

    def setupModel(self):
        # Setup model
        self.model = self.createModel()
        self.model.setExperiment(self.experiment)
        self.setupFromFormParameters()
        self.getXYvars()
        self.model.setXVar(self.varNameX)
        self.model.setYVar(self.varNameY)

    def calculateParameterUnits(self,sample):
        self.parameterUnits = self.model.calculateParameterUnits(sample)

    def getParameterNames(self):
        return self.model.getParameterNames()

    def getNumberOfParameters(self):
        return len(self.getParameterNames())

    def setParameters(self,prm):
        self.model.setParameters(prm)

    def forwardModel(self,prm,xValues):
        return self.model.forwardModel(prm,xValues)

    def setupFromFormParameters(self):
        pass

    def prepareForAnalysis(self):
        pass

    def prepareForSampleAnalysis(self, sampleName):
        pass

    def postSampleAnalysis(self, sampleName):
        pass

    def postAnalysis(self):
        pass

    def runFit(self, objId, otherDependencies):
        self.getXYvars()
        if hasattr(self,"reportX"):
            reportX = parseRange(self.reportX.get())
        else:
            reportX = None
        self.experiment = self.readExperiment(self.getInputExperiment().fnPKPD)

        # Setup model
        self.printSection("Model setup")
        self.model = self.createModel()
        self.model.setExperiment(self.experiment)
        self.model.setXVar(self.varNameX)
        self.model.setYVar(self.varNameY)
        self.setupFromFormParameters()
        self.model.printSetup()

        # Create output object
        self.fitting = PKPDFitting()
        self.fitting.fnExperiment.set(self.getInputExperiment().fnPKPD.get())
        self.fitting.predictor=self.experiment.variables[self.varNameX]
        self.fitting.predicted=self.experiment.variables[self.varNameY]
        self.fitting.modelDescription=self.model.getDescription()
        self.fitting.modelParameters = self.model.getParameterNames()
        self.fitting.modelParameterUnits = None

        # Actual fitting
        if self.fitType.get()==0:
            fitType = "linear"
        elif self.fitType.get()==1:
            fitType = "log"
        elif self.fitType.get()==2:
            fitType = "relative"

        R2List=[]
        R2adjList=[]
        AICList=[]
        AICcList=[]
        BICList=[]
        self.prepareForAnalysis()
        for sampleName, sample in self.experiment.samples.iteritems():
            self.printSection("Fitting "+sampleName)
            x, y = sample.getXYValues(self.varNameX,self.varNameY)
            print("X= "+str(x))
            print("Y= "+str(y))
            print(" ")
            self.model.setBounds(self.bounds.get())
            self.model.setXYValues(x, y)
            self.prepareForSampleAnalysis(sampleName)
            self.model.calculateParameterUnits(sample)
            if self.fitting.modelParameterUnits==None:
                self.fitting.modelParameterUnits = self.model.parameterUnits
            self.model.prepare()
            if self.model.bounds == None:
                continue
            print(" ")

            optimizer1 = PKPDDEOptimizer(self.model,fitType)
            optimizer1.optimize()
            optimizer2 = PKPDLSOptimizer(self.model,fitType)
            optimizer2.optimize()
            optimizer2.setConfidenceInterval(self.confidenceInterval.get())
            self.setParameters(optimizer2.optimum)
            optimizer2.evaluateQuality()

            # Keep this result
            sampleFit = PKPDSampleFit()
            sampleFit.sampleName = sample.sampleName
            sampleFit.x = self.model.x
            sampleFit.y = self.model.y
            sampleFit.yp = self.model.yPredicted
            sampleFit.yl = self.model.yPredictedLower
            sampleFit.yu = self.model.yPredictedUpper
            sampleFit.parameters = self.model.parameters
            sampleFit.modelEquation = self.model.getEquation()
            sampleFit.copyFromOptimizer(optimizer2)
            self.fitting.sampleFits.append(sampleFit)

            R2List.append(sampleFit.R2)
            R2adjList.append(sampleFit.R2adj)
            AICList.append(sampleFit.AIC)
            AICcList.append(sampleFit.AICc)
            BICList.append(sampleFit.BIC)

            # Add the parameters to the sample and experiment
            for varName, varUnits, description, varValue in izip(self.model.getParameterNames(), self.model.parameterUnits, self.model.getParameterDescriptions(), self.model.parameters):
                self.experiment.addParameterToSample(sampleName, varName, varUnits, description, varValue)

            self.postSampleAnalysis(sampleName)

            if reportX!=None:
                print("Evaluation of the model at specified time points")
                yreportX = self.model.forwardModel(self.model.parameters, reportX)
                print("==========================================")
                print("X     Ypredicted     log10(Ypredicted)")
                print("==========================================")
                for n in range(0,reportX.shape[0]):
                    print("%f %f %f"%(reportX[n],yreportX[n],math.log10(yreportX[n])))
                print(' ')

        self.fitting.write(self._getPath("fitting.pkpd"))
        self.experiment.write(self._getPath("experiment.pkpd"))

        fnSummary = self._getPath("summary.txt")
        fh=open(fnSummary,"w")
        fh.write("R2    (Mean+-Std): (%f)+-(%f)\n"%(np.mean(R2List),np.std(R2List)))
        fh.write("R2adj (Mean+-Std): (%f)+-(%f)\n"%(np.mean(R2adjList),np.std(R2adjList)))
        fh.write("AIC   (Mean+-Std): (%f)+-(%f)\n"%(np.mean(AICList),np.std(AICList)))
        fh.write("AICc  (Mean+-Std): (%f)+-(%f) Recommended\n"%(np.mean(AICcList),np.std(AICcList)))
        fh.write("BIC   (Mean+-Std): (%f)+-(%f)\n"%(np.mean(BICList),np.std(BICList)))
        fh.close()

        self.postAnalysis()

    def createOutputStep(self):
        self._defineOutputs(outputFitting=self.fitting)
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.getInputExperiment(), self.fitting)
        self._defineSourceRelation(self.getInputExperiment(), self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        self.getXYvars()
        msg=['Predicting %s from %s'%(self.varNameX,self.varNameY)]
        self.addFileContentToMessage(msg, self._getPath("summary.txt"))
        return msg

    def _validate(self):
        self.getXYvars()
        errors=[]
        experiment = self.readExperiment(self.getInputExperiment().fnPKPD, False)
        if not self.varNameX in experiment.variables:
            errors.append("Cannot find %s as variable"%self.varNameX)
        if not self.varNameY in experiment.variables:
            errors.append("Cannot find %s as variable"%self.varNameY)
        return errors

    def _citations(self):
        return ['Spiess2010']

    def filterVarForWizard(self, v):
        """ Define the type of variables required (used in wizard). """
        return v.isNumeric() and (v.isMeasurement() or v.isTime())