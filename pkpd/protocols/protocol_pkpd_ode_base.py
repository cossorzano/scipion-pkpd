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

import copy
import math
from itertools import izip
from collections import OrderedDict
import numpy as np

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import (PKPDDEOptimizer, PKPDLSOptimizer, PKPDFitting,
                          PKPDSampleFit, PKPDModelBase, PKPDModelBase2)
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pkpd.utils import parseRange
from pkpd.biopharmaceutics import DrugSource


class ProtPKPDODEBase(ProtPKPD,PKPDModelBase2):
    """ Base ODE protocol"""

    def __init__(self,**kwargs):
        ProtPKPD.__init__(self,**kwargs)
        self.boundsList = None

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams1(self, form, addXY=False, defaultPredictor="", defaultPredicted=""):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        if addXY:
            form.addParam('predictor', params.StringParam, label="Predictor variable (X)", default=defaultPredictor,
                          help='Y is predicted as an exponential function of X, Y=f(X)')
            form.addParam('predicted', params.StringParam, label="Predicted variable (Y)", default=defaultPredicted,
                          help='Y is predicted as an exponential function of X, Y=f(X)')

        fromTo = form.addLine('Simulation length', expertLevel = LEVEL_ADVANCED,
                           help='Minimum and maximum time (in hours) and step size (in minutes). '
                                'If minimum and maximum are not given (set to -1), they are estimated from the sample')
        fromTo.addParam('t0', params.StringParam, default="", label='Min (h)')
        fromTo.addParam('tF', params.StringParam, default="", label='Max (h)')
        fromTo.addParam('deltaT', params.FloatParam, default=0.5, label='Step (min)')

        form.addParam('fitType', params.EnumParam, choices=["Linear","Logarithmic","Relative"], label="Fit mode", default=1,
                      expertLevel=LEVEL_ADVANCED,
                      help='Linear: sum (Cobserved-Cpredicted)^2\nLogarithmic: sum(log10(Cobserved)-log10(Cpredicted))^2\n'\
                           "Relative: sum ((Cobserved-Cpredicted)/Cobserved)^2")
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')
        form.addParam('reportX', params.StringParam, label="Evaluate at X", default="", expertLevel=LEVEL_ADVANCED,
                      help='Evaluate the model at these X values\nExample 1: [0,5,10,20,40,100]\nExample 2: 0:0.55:10, from 0 to 10 in steps of 0.5')
        form.addParam('globalSearch', params.BooleanParam, label="Global search", default=True, expertLevel=LEVEL_ADVANCED,
                      help='Global search looks for the best parameters within bounds. If it is not performed, the '
                           'middle of the bounding box is used as initial parameter for a local optimization')

    #--------------------------- INSERT steps functions --------------------------------------------
    def getListOfFormDependencies(self):
        retval = [self.fitType.get(), self.confidenceInterval.get(), self.reportX.get()]
        if hasattr(self,"predictor"):
            retval.append(self.predictor.get())
            retval.append(self.predicted.get())
        if hasattr(self,"bounds"):
            retval.append(self.bounds.get())
        return retval

    def _insertAllSteps(self):
        self._insertFunctionStep('runFit',self.getInputExperiment().getObjId(),self.getListOfFormDependencies())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def getInputExperiment(self):
        if hasattr(self,"inputExperiment"):
            return self.inputExperiment.get()
        else:
            return None

    def getConfidenceInterval(self):
        return self.confidenceInterval.get()

    def getXYvars(self):
        if hasattr(self,"predictor"):
            self.varNameX=self.predictor.get()
        else:
            self.varNameX=None

        if hasattr(self,"predicted"):
            self.varNameY=self.predicted.get()
        else:
            self.varNameY=None

    def configureSource(self, drugSource):
        pass

    def createModel(self):
        pass

    def setupModel(self):
        # Setup model
        self.model = self.createModel()
        self.model.setExperiment(self.experiment)
        self.model.setXVar(self.varNameX)
        self.model.setYVar(self.varNameY)
        self.modelList.append(self.model)

    def getResponseDimension(self):
        return self.model.getResponseDimension()

    def getStateDimension(self):
        return self.model.getStateDimension()

    def setBounds(self, sample):
        if hasattr(self,"bounds"):
            self.model.setBounds(self.bounds.get())

    def prepareForSampleAnalysis(self, sampleName):
        pass

    def postSampleAnalysis(self, sampleName):
        pass

    def setTimeRange(self, sample):
        if self.t0.get()=="" or self.tF.get()=="":
            tmin, tmax = sample.getRange(self.varNameX)
        else:
            tmin = self.t0.get()
            tmax = self.tF.get()
        if self.tmin==None:
            self.tmin=tmin
        else:
            self.tmin=min(tmin,self.tmin)
        if self.tmax==None:
            self.tmax=tmax
        else:
            self.tmax=max(tmax,self.tmax)

        if self.t0.get()=="":
            self.model.t0 = min(-10,self.tmin-10) # 10 minutes before
        else:
            self.model.t0 = min(-10,float(self.t0.get())*60)

        if self.tF.get()=="":
            self.model.tF = self.tmax+10 # 10 minutes later
        else:
            self.model.tF = float(self.tF.get())*60

        if hasattr(self,"deltaT"):
            self.model.deltaT = self.deltaT.get()

    def setVarNames(self,varNameX,varNameY):
        self.varNameX = varNameX
        self.varNameY = varNameY

    def setModel(self,model):
        self.model = model

    # As model --------------------------------------------
    def clearGroupParameters(self):
        self.tmin = None
        self.tmax = None
        self.sampleList = []
        self.modelList = []
        self.drugSourceList = []
        self.clearXYLists()

    def clearXYLists(self):
        self.XList = []
        self.YList = []

    def parseBounds(self, boundsString):
        self.boundsList = []

        if boundsString!="" and boundsString!=None:
            tokens = boundsString.split(';')
            if len(tokens)!=self.getNumberOfParameters():
                raise Exception("The number of bound intervals (%d) does not match the number of parameters (%d)"%\
                                (len(tokens),self.getNumberOfParameters()))

            for token in tokens:
                values = token.strip().split(',')
                self.boundsList.append((float(values[0][1:]),float(values[1][:-1])))

    def setBounds(self, sample):
        self.parseBounds(self.bounds.get())
        self.setBoundsFromBoundsList()

    def setBoundsFromBoundsList(self):
        Nbounds = len(self.boundsList)
        Nsource = self.drugSource.getNumberOfParameters()
        Nmodel = self.model.getNumberOfParameters()
        if Nbounds!=Nsource+Nmodel:
            raise Exception("The number of parameters (%d) and bounds (%d) are different"%(Nsource+Nmodel,Nbounds))
        self.boundsSource = self.boundsList[0:Nsource]
        self.boundsPK = self.boundsList[Nsource:]
        self.model.bounds = self.boundsPK

    def getBounds(self):
        return self.boundsList

    def getParameterBounds(self):
        """ Return a dictionary where the parameter name is the key
        and the bounds are its values. """
        boundsDict = OrderedDict()
        self.parseBounds(self.bounds.get()) # after this we have boundsList
        parameterNames = self.getParameterNames()

        for paramName, bound in izip(parameterNames, self.getBounds()):
            boundsDict[paramName] = bound

        # Set None as bound for parameters not matched
        for paramName in parameterNames:
            if paramName not in boundsDict:
                boundsDict[paramName] = None

        return boundsDict

    def setParameters(self, parameters):
        self.parameters = parameters
        self.parametersPK = self.parameters[-self.NparametersModel:]

        for n in range(len(self.modelList)):
            if self.NparametersSource>0:
                self.drugSourceList[n].setParameters(self.parameters[0:self.NparametersSource])
            self.modelList[n].setParameters(self.parametersPK)

    def setXYValues(self, x, y):
        PKPDModelBase.setXYValues(self,x,y)
        self.model.setXYValues(x, y)
        self.XList.append(x)
        self.YList.append(y)

    def mergeLists(self, iny):
        if len(self.XList)>1:
            outy=[]
            for m in range(len(iny[0])):
                outy.append(np.empty(shape=[0]))
            for m in range(len(outy)):
                for n in range(len(iny)):
                    outy[m]=np.concatenate((outy[m],iny[n][m]))
            return outy
        else:
            return iny[0]

    def separateLists(self,iny):
        outy=[]
        Nsamples=len(self.YList)
        if Nsamples==0:
            return
        Nmeasurements = self.model.getResponseDimension()
        idx=[0]*Nmeasurements
        for n in range(Nsamples):
            yn=self.YList[n]
            perSampleIn = []
            for j in range(Nmeasurements):
                ynDim = yn[j].size
                ysample = iny[j][idx[j]:(idx[j]+ynDim)]
                perSampleIn.append(ysample)
                idx[j]+=ynDim
            outy.append(perSampleIn)
        return outy

    def addSample(self, sample):
        self.sampleList.append(sample)
        self.model.setSample(sample)

    def forwardModel(self, parameters, x=None):
        self.setParameters(parameters)

        yPredictedList = []
        for n in range(len(self.modelList)):
            self.modelList[n].forwardModel(self.parametersPK,x)
            yPredictedList.append(self.modelList[n].yPredicted)
        self.yPredicted = self.mergeLists(yPredictedList)
        return copy.copy(self.yPredicted)

    def forwardModelByConvolution(self, parameters, x=None):
        self.setParameters(parameters)
        tFImpulse = None
        if hasattr(self,"tFImpulse"):
            if self.tFImpulse.get()!="":
                tFImpulse = float(self.tFImpulse.get())

        yPredictedList = []
        for n in range(len(self.modelList)):
            self.modelList[n].tFImpulse = tFImpulse
            self.modelList[n].forwardModelByConvolution(self.parametersPK,x)
            yPredictedList.append(self.modelList[n].yPredicted)
        self.yPredicted = self.mergeLists(yPredictedList)
        return copy.copy(self.yPredicted)

    def imposeConstraints(self,yt):
        self.model.imposeConstraints(yt)

    def getEquation(self):
        return self.drugSource.getEquation()+" and "+self.model.getEquation()

    def getModelEquation(self):
        return self.drugSource.getModelEquation()+" and "+self.model.getModelEquation()

    def getDescription(self):
        return self.drugSource.getDescription()+"; "+self.model.getDescription()

    def getParameterNames(self):
        retval = []
        parametersSource = self.drugSource.getParameterNames()
        self.NparametersSource = len(parametersSource)
        retval += parametersSource
        parametersModel = self.model.getParameterNames()
        self.NparametersModel = len(parametersModel)
        retval += parametersModel
        return retval

    def calculateParameterUnits(self,sample):
        retval = []
        retval += self.drugSource.calculateParameterUnits(sample)
        retval += self.model.calculateParameterUnits(sample)
        self.parameterUnits = retval

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        idx=0
        if len(self.boundsSource)>0:
            retval+=self.drugSource.areParametersSignificant(lowerBound[idx:len(self.boundsSource)],
                                                             upperBound[idx:len(self.boundsSource)])
        retval+=self.model.areParametersSignificant(lowerBound[len(self.boundsSource):],
                                                    upperBound[len(self.boundsSource):])
        return retval

    def areParametersValid(self, p):
        return self.drugSource.areParametersValid(p[0:len(self.boundsSource)]) and \
               self.model.areParametersValid(p[len(self.boundsSource):])

    def createDrugSource(self):
        self.drugSource = DrugSource()
        self.drugSourceList.append(self.drugSource)

        return self.drugSource

    # Really fit ---------------------------------------------------------
    def runFit(self, objId, otherDependencies):
        reportX = parseRange(self.reportX.get())
        self.setInputExperiment()

        # Setup model
        self.getXYvars()

        # Create output object
        self.fitting = PKPDFitting()
        self.fitting.fnExperiment.set(self._getPath("experiment.pkpd"))
        self.fitting.predictor=self.experiment.variables[self.varNameX]
        if type(self.varNameY)==list:
            self.fitting.predicted=[]
            for y in self.varNameY:
                self.fitting.predicted.append(self.experiment.variables[y])
        else:
            self.fitting.predicted=self.experiment.variables[self.varNameY]
        self.fitting.modelParameterUnits = None

        # Actual fitting
        if self.fitType.get()==0:
            fitType = "linear"
        elif self.fitType.get()==1:
            fitType = "log"
        elif self.fitType.get()==2:
            fitType = "relative"

        for groupName, group in self.experiment.groups.iteritems():
            self.printSection("Fitting "+groupName)
            self.clearGroupParameters()

            for sampleName in group.sampleList:
                print("   Sample "+sampleName)
                sample = self.experiment.samples[sampleName]

                self.createDrugSource()
                self.setupModel()

                # Get the values to fit
                x, y = sample.getXYValues(self.varNameX,self.varNameY)
                print("X= "+str(x))
                print("Y= "+str(y))
                print(" ")

                # Interpret the dose
                self.setTimeRange(sample)
                sample.interpretDose()

                self.drugSource.setDoses(sample.parsedDoseList, self.model.t0, self.model.tF)
                self.configureSource(self.drugSource)
                self.model.drugSource = self.drugSource

                # Prepare the model
                self.setBounds(sample)
                self.setXYValues(x, y)
                self.addSample(sample)
                self.prepareForSampleAnalysis(sampleName)
                self.calculateParameterUnits(sample)
                if self.fitting.modelParameterUnits==None:
                    self.fitting.modelParameterUnits = self.parameterUnits

            self.printSetup()
            self.x = self.mergeLists(self.XList)
            self.y = self.mergeLists(self.YList)

            if self.globalSearch:
                optimizer1 = PKPDDEOptimizer(self,fitType)
                optimizer1.optimize()
            else:
                self.parameters = np.zeros(len(self.boundsList),np.double)
                n = 0
                for bound in self.boundsList:
                    self.parameters[n] = 0.5*(bound[0]+bound[1])
                    n += 1
            try:
                optimizer2 = PKPDLSOptimizer(self,fitType)
                optimizer2.optimize()
            except Exception as e:
                msg="Error: "+str(e)
                msg+="\nErrors in the local optimizer may be caused by starting from a bad initial guess\n"
                msg+="Try performing a global search first or changing the bounding box"
                raise Exception("Error in the local optimizer\n"+msg)
            optimizer2.setConfidenceInterval(self.getConfidenceInterval())
            self.setParameters(optimizer2.optimum)
            optimizer2.evaluateQuality()
            self.model.printOtherParameterization()

            self.yPredictedList=self.separateLists(self.yPredicted)
            self.yPredictedLowerList=self.separateLists(self.yPredictedLower)
            self.yPredictedUpperList=self.separateLists(self.yPredictedUpper)

            n=0
            for sampleName in group.sampleList:
                sample = self.experiment.samples[sampleName]

                # Keep this result
                sampleFit = PKPDSampleFit()
                sampleFit.sampleName = sample.sampleName
                sampleFit.x = self.XList[n]
                sampleFit.y = self.YList[n]
                sampleFit.yp = self.yPredictedList[n]
                sampleFit.yl = self.yPredictedLowerList[n]
                sampleFit.yu = self.yPredictedUpperList[n]
                sampleFit.parameters = self.parameters
                sampleFit.modelEquation = self.getEquation()
                sampleFit.copyFromOptimizer(optimizer2)
                self.fitting.sampleFits.append(sampleFit)

                # Add the parameters to the sample and experiment
                for varName, varUnits, description, varValue in izip(self.getParameterNames(), self.parameterUnits, self.getParameterDescriptions(), self.parameters):
                    self.experiment.addParameterToSample(sampleName, varName, varUnits, description, varValue, rewrite=True)

                self.postSampleAnalysis(sampleName)

                if reportX!=None:
                    print("Evaluation of the model at specified time points")
                    self.model.tF = np.max(reportX)
                    yreportX = self.model.forwardModel(self.model.parameters, reportX)
                    print("==========================================")
                    print("X     Ypredicted     log10(Ypredicted)")
                    print("==========================================")
                    for n in range(0,reportX.shape[0]):
                        aux = 0
                        if yreportX[n]>0:
                            aux = math.log10(yreportX[n])
                        print("%f %f %f"%(reportX[n],yreportX[n],aux))
                    print(' ')

                n+=1

        self.fitting.modelParameters = self.getParameterNames()
        self.fitting.modelDescription=self.getDescription()
        self.fitting.write(self._getPath("fitting.pkpd"))
        self.experiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputFitting=self.fitting)
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.getInputExperiment(), self.fitting)
        self._defineSourceRelation(self.getInputExperiment(), self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        self.getXYvars()
        if self.varNameX!=None:
            msg.append('Predicting %s from %s'%(self.varNameX,self.varNameY))
        return msg

    def _validate(self):
        self.getXYvars()
        errors=[]
        if self.varNameX!=None:
            experiment = self.readExperiment(self.getInputExperiment().fnPKPD, False)
            if not self.varNameX in experiment.variables:
                errors.append("Cannot find %s as variable"%self.varNameX)
            if type(self.varNameY)==list:
                for y in self.varNameY:
                    if not y in experiment.variables:
                        errors.append("Cannot find %s as variable"%y)
            else:
                if not self.varNameY in experiment.variables:
                    errors.append("Cannot find %s as variable"%self.varNameY)
        if self.bounds.get()=="":
            errors.append("Bounds are required")
        return errors

    def _citations(self):
        return ['Spiess2010']

    def filterVarForWizard(self, v):
        """ Define the type of variables required (used in wizard). """
        return v.isNumeric() and (v.isTime() or v.isMeasurement())
