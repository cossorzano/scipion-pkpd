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

from itertools import izip

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDFitting, PKPDSampleFit, PKPDLSOptimizer
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from .protocol_pkpd_ode_base import ProtPKPDODEBase

class ProtPKPDODERefine(ProtPKPDODEBase):
    """ Refinement of an ODE protocol. The parameters are reestimated with a finer sampling rate.
    """

    _label = 'ODE refinement'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputODE', params.PointerParam, label="Input ODE model",
                      pointerClass='ProtPKPDMonoCompartment, ProtPKPDMonoCompartmentUrine, ProtPKPDTwoCompartments, '\
                                   'ProtPKPDTwoCompartmentsAutoinduction, ProtPKPDTwoCompartmentsClint, '\
                                   'ProtPKPDTwoCompartmentsClintMetabolite, ProtPKPDTwoCompartmentsUrine, '\
                                   'ProtPKPDODERefine',
                      help='Select a run of an ODE model')
        form.addParam('deltaT', params.FloatParam, default=0.5, label='Step (min)', expertLevel=LEVEL_ADVANCED)
        form.addParam('fitType', params.EnumParam, choices=["Linear","Logarithmic","Relative","Same as previous protocol"], label="Fit mode", default=3,
                      help='Linear: sum (Cobserved-Cpredicted)^2\nLogarithmic: sum(log10(Cobserved)-log10(Cpredicted))^2\n'\
                           "Relative: sum ((Cobserved-Cpredicted)/Cobserved)^2")
        form.addParam('bounds', params.StringParam, label="Parameter bounds", default="",
                      help="If empty, same bounds as for the previous protocol")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runFit',self.inputODE.get().getObjId(), self.deltaT.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def parseBounds(self, boundsString):
        if boundsString!="" and boundsString!=None:
            tokens=boundsString.split(';')
            if len(tokens)!=self.getNumberOfParameters():
                raise Exception("The number of bound intervals does not match the number of parameters")
            self.boundsList=[]
            for token in tokens:
                values = token.strip().split(',')
                self.boundsList.append((float(values[0][1:]),float(values[1][:-1])))

    def setBounds(self,sample):
        self.parseBounds(self.protODE.bounds.get() if self.bounds.get()=="" else self.bounds.get())
        self.setBoundsFromBoundsList()

    def setBoundsFromBoundsList(self):
        Nbounds = len(self.boundsList)
        Nsource = self.drugSource.getNumberOfParameters()
        Nmodel = self.model.getNumberOfParameters()
        if Nbounds!=Nsource+Nmodel:
            raise "The number of parameters (%d) and bounds (%d) are different"%(Nsource+Nmodel,Nbounds)
        self.boundsSource = self.boundsList[0:Nsource]
        self.boundsPK = self.boundsList[Nsource:]
        self.model.bounds = self.boundsPK

    def getBounds(self):
        return self.boundsList

    def createModel(self):
        return self.inputODE.get().createModel()

    def setTimeRange(self, sample):
        self.inputODE.get().setTimeRange(sample)

    def setVarNames(self,varNameX,varNameY):
        ProtPKPDODEBase.setVarNames(self,varNameX,varNameY)
        self.inputODE.get().setVarNames(varNameX,varNameY)

    def setModel(self,model):
        ProtPKPDODEBase.setModel(self,model)
        self.inputODE.get().setModel(model)

    def setupModel(self):
        ProtPKPDODEBase.setupModel(self)
        self.inputODE.get().setModel(self.model)

    def clearGroupParameters(self):
        ProtPKPDODEBase.clearGroupParameters(self)
        self.inputODE.get().clearGroupParameters()

    def getConfidenceInterval(self):
        return self.inputODE.get().getConfidenceInterval()

    def runFit(self, objId, deltaT):
        self.protODE = self.inputODE.get()
        self.experiment = self.readExperiment(self.protODE.outputExperiment.fnPKPD)
        self.fitting = self.readFitting(self.protODE.outputFitting.fnFitting)

        # Get the X and Y variable names
        self.varNameX = self.fitting.predictor.varName
        if type(self.fitting.predicted)==list:
            self.varNameY = [v.varName for v in self.fitting.predicted]
        else:
            self.varNameY = self.fitting.predicted.varName
        self.protODE.experiment = self.experiment
        self.protODE.setVarNames(self.varNameX, self.varNameY)

        # Create output object
        self.fitting = PKPDFitting()
        self.fitting.fnExperiment.set(self.experiment.fnPKPD.get())
        self.fitting.predictor=self.experiment.variables[self.varNameX]
        self.fitting.predicted=self.experiment.variables[self.varNameY]
        self.fitting.modelParameterUnits = None

        # Actual fitting
        fitTypeN = self.protODE.fitType.get() if self.fitType.get()==3 else self.fitType.get()
        if fitTypeN==0:
            fitType = "linear"
        elif fitTypeN==1:
            fitType = "log"
        elif fitTypeN==2:
            fitType = "relative"

        parameterNames = None
        for groupName, group in self.experiment.groups.iteritems():
            self.printSection("Fitting "+groupName)
            self.protODE.clearGroupParameters()
            self.clearGroupParameters()

            for sampleName in group.sampleList:
                print("   Sample "+sampleName)
                sample = self.experiment.samples[sampleName]

                self.protODE.createDrugSource()
                self.protODE.setupModel()

                # Setup self model
                self.drugSource = self.protODE.drugSource
                self.drugSourceList = self.protODE.drugSourceList
                self.model = self.protODE.model
                self.modelList = self.protODE.modelList
                self.model.deltaT = self.deltaT.get()
                self.model.setXVar(self.varNameX)
                self.model.setYVar(self.varNameY)

                # Get the values to fit
                x, y = sample.getXYValues(self.varNameX,self.varNameY)
                print("X= "+str(x))
                print("Y= "+str(y))

                # Interpret the dose
                self.protODE.varNameX = self.varNameX
                self.protODE.varNameY = self.varNameY
                self.protODE.model = self.model
                self.protODE.setTimeRange(sample)
                sample.interpretDose()

                self.drugSource.setDoses(sample.parsedDoseList, self.model.t0, self.model.tF)
                self.protODE.configureSource(self.drugSource)
                self.model.drugSource = self.drugSource

                # Prepare the model
                self.model.setSample(sample)
                self.calculateParameterUnits(sample)
                if self.fitting.modelParameterUnits==None:
                    self.fitting.modelParameterUnits = self.parameterUnits

                # Get the initial parameters
                if parameterNames==None:
                    parameterNames = self.getParameterNames()
                parameters0 = []
                for parameterName in parameterNames:
                    parameters0.append(float(sample.descriptors[parameterName]))
                print("Initial solution: %s"%str(parameters0))
                print(" ")

                # Set bounds
                self.setBounds(sample)
                self.setXYValues(x, y)
                self.parameters = parameters0

            self.printSetup()
            self.x = self.mergeLists(self.XList)
            self.y = self.mergeLists(self.YList)

            optimizer2 = PKPDLSOptimizer(self,fitType)
            optimizer2.optimize()
            optimizer2.setConfidenceInterval(self.protODE.getConfidenceInterval())
            self.setParameters(optimizer2.optimum)

            n=0
            for sampleName in group.sampleList:
                sample = self.experiment.samples[sampleName]

                # Keep this result
                sampleFit = PKPDSampleFit()
                sampleFit.sampleName = sample.sampleName
                sampleFit.x = x
                sampleFit.y = y
                sampleFit.yp = self.yPredicted
                sampleFit.yl = self.yPredictedLower
                sampleFit.yu = self.yPredictedUpper
                sampleFit.parameters = self.parameters
                sampleFit.modelEquation = self.getEquation()
                sampleFit.copyFromOptimizer(optimizer2)
                self.fitting.sampleFits.append(sampleFit)

                # Add the parameters to the sample and experiment
                for varName, varUnits, description, varValue in izip(self.getParameterNames(), self.parameterUnits, self.getParameterDescriptions(), self.parameters):
                    self.experiment.addParameterToSample(sampleName, varName, varUnits, description, varValue, rewrite=True)

                n+=1

        self.fitting.modelParameters = self.getParameterNames()
        self.fitting.modelDescription = self.getDescription()
        self.fitting.write(self._getPath("fitting.pkpd"))
        self.experiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputFitting=self.fitting)
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputODE.get(), self.fitting)
        self._defineSourceRelation(self.inputODE.get(), self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        msg.append("New sampling rate: %f"%self.deltaT.get())
        return msg

    def _validate(self):
        return []
