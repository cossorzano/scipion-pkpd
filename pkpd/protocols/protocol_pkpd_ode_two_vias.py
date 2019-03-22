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
from .protocol_pkpd import ProtPKPD
from pkpd.objects import (PKPDModelBase2, PKPDExperiment, PKPDFitting,
                          PKPDDEOptimizer, PKPDLSOptimizer, flattenArray,
                          PKPDSampleFit)
from pyworkflow.protocol.constants import LEVEL_ADVANCED

# TESTED in test_workflow_gabrielsson_pk10.py

class ProtPKPDODETwoVias(ProtPKPD,PKPDModelBase2):
    """ Simultaneous fit of data obtained by different vias, e.g. IV and PO, but it can be any two vias and any two
        dosing regimes, dissolution profiles, etc. It is supposed that the PK model in both cases is the same
        (e.g. two monocompartments, two two-compartments, ... """
    _label = 'ode two vias'

    def __init__(self,**kwargs):
        ProtPKPD.__init__(self,**kwargs)
        self.boundsList = None
        self.parameterNames = None

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('prot1ptr', params.PointerParam, label="Model Via 1 (e.g. IV)",
                      pointerClass='ProtPKPDMonoCompartment, ProtPKPDTwoCompartments, ProtPKPDTwoCompartmentsBoth, ProtPKPDMonoCompartmentPD, ProtPKPDTwoCompartmentsBothPD',
                      help='Select the model for the intravenous route')
        form.addParam('prot2ptr', params.PointerParam, label="Model Via 2 (e.g. PO)",
                      pointerClass='ProtPKPDMonoCompartment, ProtPKPDTwoCompartments, ProtPKPDTwoCompartmentsBoth, ProtPKPDMonoCompartmentPD, ProtPKPDTwoCompartmentsBothPD',
                      help='Select the model for the oral route. It must be of the same type as the intravenous,'
                           'and the sample names must be the same in both experiments')
        form.addParam('fitType', params.EnumParam, choices=["Linear","Logarithmic","Relative"], label="Fit mode", default=1,
                      expertLevel=LEVEL_ADVANCED,
                      help='Linear: sum (Cobserved-Cpredicted)^2\nLogarithmic: sum(log10(Cobserved)-log10(Cpredicted))^2\n'\
                           "Relative: sum ((Cobserved-Cpredicted)/Cobserved)^2")
        form.addParam('globalSearch', params.BooleanParam, label="Global search", default=False, expertLevel=LEVEL_ADVANCED,
                      help='Global search looks for the best parameters within bounds. If it is not performed, the '
                           'middle of the bounding box is used as initial parameter for a local optimization')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runFit',self.prot1ptr.get().outputExperiment.fnPKPD,
                                 self.prot2ptr.get().outputExperiment.fnPKPD)
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    # As model --------------------------------------------
    def setBounds(self, sample1, sample2):
        self.prot1.parseBounds(self.prot1.bounds.get())
        self.prot1.setBoundsFromBoundsList()

        self.prot2.parseBounds(self.prot2.bounds.get())
        self.prot2.setBoundsFromBoundsList()

    def setXYValues(self, x1, y1, x2, y2):
        self.prot1.setXYValues(x1,y1)
        self.prot2.setXYValues(x2,y2)
        self.x = x1+x2
        self.y = y1+y2

    def addSample(self, sample1, sample2):
        self.prot1.addSample(sample1)
        self.prot2.addSample(sample2)

    def prepareForSampleAnalysis(self, sampleName):
        self.prot1.prepareForSampleAnalysis(sampleName)
        self.prot2.prepareForSampleAnalysis(sampleName)

    def mergeModelParameters(self):
        self.parameterNames = []
        self.parameterProcedence = []
        self.parameterNames1=self.prot1.getParameterNames()
        self.parameterNames2=self.prot2.getParameterNames()
        i=0
        for name in self.parameterNames1:
            self.parameterNames.append(name)
            self.parameterProcedence.append([(self.prot1,i)])
            i+=1
        self.N1=i
        i=0
        for name in self.parameterNames2:
            if not name in self.parameterNames:
                self.parameterNames.append(name)
                self.parameterProcedence.append([(self.prot2,i)])
            else:
                idx = self.parameterNames.index(name)
                self.parameterProcedence[idx].append((self.prot2,i))
            i+=1
        self.N2=i

    def calculateParameterUnits(self, sample1, sample2):
        self.prot1.calculateParameterUnits(sample1)
        self.prot2.calculateParameterUnits(sample2)

        self.parameterUnits = []
        for protiList in self.parameterProcedence:
            prot, i = protiList[0]
            self.parameterUnits.append(prot.parameterUnits[i])

    def splitParameterUnits(self):
        parameterNames1=self.prot1.getParameterNames()
        parameterNames2=self.prot2.getParameterNames()
        parameterNames=self.getParameterNames()
        self.fitting1.modelParameterUnits=[]
        for name in parameterNames1:
            idx = parameterNames.index(name)
            self.fitting1.modelParameterUnits.append(self.parameterUnits[idx])
        self.fitting1.modelParameters = parameterNames1
        self.fitting1.modelDescription=self.prot1.getDescription()

        self.fitting2.modelParameterUnits=[]
        for name in parameterNames2:
            idx = parameterNames.index(name)
            self.fitting2.modelParameterUnits.append(self.parameterUnits[idx])
        self.fitting2.modelParameters = parameterNames2
        self.fitting2.modelDescription=self.prot2.getDescription()

    def getModelEquation(self):
        return "Via1: "+self.prot1.getModelEquation()+" and Via2: "+self.prot2.getModelEquation()

    def getEquation(self):
        return "Via1: "+self.prot1.getEquation()+" and Via2: "+self.prot2.getEquation()

    def getParameterNames(self):
        return self.parameterNames

    def getBounds(self):
        retval = []
        for protiList in self.parameterProcedence:
            if len(protiList)==1:
                prot, i = protiList[0]
                retval.append(prot.boundsList[i])
            else:
                prot1, i1 = protiList[0]
                prot2, i2 = protiList[1]
                m1,M1 = prot1.boundsList[i1]
                m2,M2 = prot2.boundsList[i2]
                retval.append((min(m1,m2),max(M1,M2)))
        return retval

    def setInitialSolution(self,sampleName):
        self.parameters=np.zeros(len(self.parameterNames),np.double)
        j=0
        for protiList in self.parameterProcedence:
            if len(protiList)==1:
                prot, i = protiList[0]
                if prot is self.prot1:
                    self.parameters[j]=float(self.experiment1.samples[sampleName].descriptors[self.parameterNames1[i]])
                else:
                    self.parameters[j]=float(self.experiment2.samples[sampleName].descriptors[self.parameterNames2[i]])
            else:
                prot1, i1 = protiList[0]
                prot2, i2 = protiList[1]
                self.parameters[j]=0.5*(
                                     float(self.experiment1.samples[sampleName].descriptors[self.parameterNames1[i1]])+
                                     float(self.experiment2.samples[sampleName].descriptors[self.parameterNames2[i2]]))
            j+=1

    def setParameters(self, parameters):
        self.parameters = parameters
        self.parameters1 = [None]*self.N1
        self.parameters2 = [None]*self.N2
        j=0
        for protiList in self.parameterProcedence:
            for prot, i in protiList:
                if prot is self.prot1:
                    self.parameters1[i]=parameters[j]
                else:
                    self.parameters2[i]=parameters[j]
            j+=1
        self.prot1.setParameters(self.parameters1)
        self.prot2.setParameters(self.parameters2)

    def separateBounds(self, lowerBound, upperBound, parameterNames1, parameterNames2, parameterNames):
        lower1=[]
        upper1=[]
        for name in parameterNames1:
            idx = parameterNames.index(name)
            lower1.append(lowerBound[idx])
            upper1.append(upperBound[idx])

        lower2=[]
        upper2=[]
        for name in parameterNames2:
            idx = parameterNames.index(name)
            lower2.append(lowerBound[idx])
            upper2.append(upperBound[idx])
        return lower1, upper1, lower2, upper2

    def areParametersSignificant(self, lowerBound, upperBound):
        parameterNames1=self.prot1.getParameterNames()
        parameterNames2=self.prot2.getParameterNames()
        parameterNames=self.getParameterNames()
        lower1, upper1, lower2, upper2 = self.separateBounds(lowerBound, upperBound, parameterNames1, parameterNames2,
                                                             parameterNames)

        significant1 = self.prot1.areParametersSignificant(lower1,upper1)
        significant2 = self.prot2.areParametersSignificant(lower2,upper2)

        retval=[True]*len(parameterNames)
        j=0
        for protiList in self.parameterProcedence:
            prot, i = protiList[0]
            if prot is self.prot1:
                idx=parameterNames1.index(parameterNames[j])
                retval[j]=significant1[idx]
            else:
                idx=parameterNames2.index(parameterNames[j])
                retval[j]=significant2[idx]
            j+=1
        return retval

    def setConfidenceInterval(self,lowerBound,upperBound):
        parameterNames1=self.prot1.getParameterNames()
        parameterNames2=self.prot2.getParameterNames()
        parameterNames=self.getParameterNames()
        self.lower1, self.upper1, self.lower2, self.upper2 = self.separateBounds(lowerBound, upperBound,
                                                                                 parameterNames1, parameterNames2,
                                                                                 parameterNames)
        self.prot1.setConfidenceInterval(self.lower1,self.upper1)
        self.prot2.setConfidenceInterval(self.lower2,self.upper2)

    def forwardModel(self, parameters, x=None):
        self.setParameters(parameters)
        self.yp1 = self.prot1.forwardModel(self.parameters1,x)
        self.yp2 = self.prot2.forwardModel(self.parameters2,x)
        self.yPredicted = self.yp1+self.yp2 # Merge two lists
        return self.yPredicted

    def keepSampleFit(self,fitting,prot,sampleName,x,y,yp,yplower,ypupper,prmLowerBound,prmUpperBound,optimizer2):
        sampleFit = PKPDSampleFit()
        sampleFit.sampleName = sampleName
        sampleFit.x = x
        sampleFit.y = y
        sampleFit.yp = yp
        sampleFit.yl = yplower
        sampleFit.yu = ypupper
        sampleFit.parameters = prot.parameters
        sampleFit.modelEquation = prot.getEquation()

        sampleFit.copyFromOptimizer(optimizer2)
        sampleFit.lowerBound = prmLowerBound
        sampleFit.upperBound = prmUpperBound

        fitting.sampleFits.append(sampleFit)

    # Really fit ---------------------------------------------------------
    def setupUnderlyingProtocol(self,prot):
        prot.getXYvars()
        prot.experiment=PKPDExperiment()
        prot.experiment.load(prot.outputExperiment.fnPKPD.get())
        prot.clearGroupParameters()
        prot.createDrugSource()
        prot.setupModel()
        return prot.experiment

    def setupSample(self,prot,sample,prefix):
        sample.sampleName = prefix+sample.sampleName
        prot.setTimeRange(sample)
        sample.interpretDose()
        prot.drugSource.setDoses(sample.parsedDoseList, prot.model.t0, prot.model.tF)
        prot.configureSource(prot.drugSource)

    def prepareDoseForSample(self,prot,sample):
        prot.setTimeRange(sample)
        sample.interpretDose()
        prot.drugSource.setDoses(sample.parsedDoseList, prot.model.t0, prot.model.tF)
        prot.configureSource(prot.drugSource)
        prot.model.drugSource = prot.drugSource

    def updateUnderlyingExperiments(self,sampleName):
        j=0
        for protiList in self.parameterProcedence:
            for prot, i in protiList:
                if prot is self.prot1:
                    self.experiment1.samples[sampleName].descriptors[self.parameterNames1[i]]=str(self.parameters[j])
                else:
                    self.experiment2.samples[sampleName].descriptors[self.parameterNames2[i]]=str(self.parameters[j])
            j+=1

    def createFitting(self, prot, experiment, suffix):
        # Create output object
        fitting = PKPDFitting()
        fitting.fnExperiment.set(self._getPath("experiment%d.pkpd"%suffix))
        fitting.predictor=experiment.variables[prot.varNameX]
        if type(prot.varNameY)==list:
            fitting.predicted=[]
            for y in prot.varNameY:
                fitting.predicted.append(experiment.variables[y])
        else:
            fitting.predicted=experiment.variables[prot.varNameY]
        fitting.modelParameterUnits = None
        return fitting

    def runFit(self, fn1, fn2):
        prot1 = self.prot1ptr.get()
        prot2 = self.prot2ptr.get()
        self.prot1 = prot1
        self.prot2 = prot2
        self.experiment1 = self.setupUnderlyingProtocol(prot1)
        self.experiment2 = self.setupUnderlyingProtocol(prot2)
        self.fitting1 = self.createFitting(self.prot1, self.experiment1, 1)
        self.fitting2 = self.createFitting(self.prot2, self.experiment2, 2)

        if self.fitType.get()==0:
            fitType = "linear"
        elif self.fitType.get()==1:
            fitType = "log"
        elif self.fitType.get()==2:
            fitType = "relative"

        # The fitting is performed by sampleName and not by groupName
        self.someInCommon = False
        for sample2name, sample2 in self.experiment2.samples.iteritems():
            if sample2name in self.experiment1.samples:
                self.someInCommon = True
                self.setupSample(prot2,sample2,"")

                sample1=self.experiment1.samples[sample2name]
                self.setupSample(prot1,sample1,"")

                # Get the values to fit
                x1, y1 = sample1.getXYValues(prot1.varNameX,prot1.varNameY)
                x2, y2 = sample2.getXYValues(prot2.varNameX,prot2.varNameY)
                print("Sample: "+sample2name)
                print("X1= "+str(x1))
                print("Y1= "+str(y1))
                print("X2= "+str(x2))
                print("Y2= "+str(y2))
                print(" ")

                # Interpret the dose
                self.prepareDoseForSample(prot1,sample1)
                self.prepareDoseForSample(prot2,sample2)

                # Prepare the model
                self.setBounds(sample1,sample2)
                self.setXYValues(x1, y1, x2, y2)
                self.addSample(sample1,sample2)
                self.prepareForSampleAnalysis(sample2name)
                self.mergeModelParameters()
                self.calculateParameterUnits(sample1, sample2)
                if self.fitting1.modelParameterUnits==None:
                    self.splitParameterUnits()
                self.printSetup()

                print(" ")

                # Optimize
                if self.globalSearch:
                    optimizer1 = PKPDDEOptimizer(self,fitType)
                    optimizer1.optimize()
                else:
                    self.setInitialSolution(sample2name)
                optimizer2 = PKPDLSOptimizer(self,fitType)
                optimizer2.optimize()
                optimizer2.setConfidenceInterval(self.prot1.confidenceInterval.get())
                self.setParameters(optimizer2.optimum)
                optimizer2.evaluateQuality()

                # Set the variables back to the experiments
                self.updateUnderlyingExperiments(sample2name)

                # Update fittings
                self.keepSampleFit(self.fitting1, self.prot1, sample2name,
                                   x1,y1,self.prot1.yPredicted,
                                   self.prot1.yPredictedLower,self.prot1.yPredictedUpper,
                                   self.lower1,self.upper1,optimizer2)
                self.keepSampleFit(self.fitting2, self.prot2, sample2name,
                                   x2,y2,self.prot2.yPredicted,
                                   self.prot2.yPredictedLower,self.prot2.yPredictedUpper,
                                   self.lower2,self.upper2,optimizer2)

        if self.someInCommon:
            self.experiment1.write(self._getPath("experiment1.pkpd"))
            self.experiment2.write(self._getPath("experiment2.pkpd"))
            self.fitting1.write(self._getPath("fitting1.pkpd"))
            self.fitting2.write(self._getPath("fitting2.pkpd"))

    def createOutputStep(self):
        if self.someInCommon:
            self._defineOutputs(outputExperiment1=self.experiment1)
            self._defineOutputs(outputExperiment2=self.experiment2)
            self._defineOutputs(outputFitting1=self.fitting1)
            self._defineOutputs(outputFitting2=self.fitting2)
            self._defineSourceRelation(self.prot1ptr.get().outputExperiment, self.experiment1)
            self._defineSourceRelation(self.prot2ptr.get().outputExperiment, self.experiment1)
            self._defineSourceRelation(self.prot1ptr.get().outputExperiment, self.experiment2)
            self._defineSourceRelation(self.prot2ptr.get().outputExperiment, self.experiment2)
            self._defineSourceRelation(self.prot1ptr.get().outputExperiment, self.fitting1)
            self._defineSourceRelation(self.prot2ptr.get().outputExperiment, self.fitting1)
            self._defineSourceRelation(self.prot1ptr.get().outputExperiment, self.fitting2)
            self._defineSourceRelation(self.prot2ptr.get().outputExperiment, self.fitting2)
        else:
            print("Error: There are no samples in common, nothing is produced")

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        return msg

# Falta construir salida experiment y fitting
