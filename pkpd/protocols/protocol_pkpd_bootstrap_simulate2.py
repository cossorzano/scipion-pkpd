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
import numpy as np
import math
import random

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDDose, PKPDSample, PKPDVariable
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pkpd.pkpd_units import createUnit, multiplyUnits, strUnit, PKPDUnit
from .protocol_pkpd_ode_base import ProtPKPDODEBase

class ProtPKPDODESimulate2(ProtPKPDODEBase):
    """ Simulate a complex sequence of doses. Each dose has its own PK profile, and the overall profile is the
        addition of all the individual profiles. Only the parameters from the first fitting of each model are
        considered.

        The first dose goes with the first model, the second dose with the second model, ...

        A single profile is simulated that is the addition of the different simulations for each one of the input
        models."""

    _label = 'PK simulate (complex)'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputODEs', params.MultiPointerParam, label="Input ODE models",
                      pointerClass='ProtPKPDMonoCompartment, ProtPKPDTwoCompartments, ProtPKPDMonoCompartmentPD, ProtPKPDTwoCompartmentsBothPD, '\
                                   'ProtPKPDODERefine', help='Select a run of an ODE model')
        form.addParam('doses', params.TextParam, label="Doses", height=5, width=50,
                      default="",
                      help="Structure: [Dose Name] ; [Via] ; [Dose type] ; [time] ; [dose] \n"\
                           "The dose name should have no space or special character\n"\
                           "The via should be one present in the input experiment to the ODE model.\n"\
                           "Valid units are: h, mg, ug, ...\n"\
                           "The description is either a bolus or an infusion as shown in the examples\n"\
                           "\nIt is important that there are two semicolons.\n"\
                           "Examples:\n"\
                           "Infusion0; via=Intravenous; infusion; t=0.500000:0.750000 h; d=60 mg/h\n"\
                           "Bolus0; via=Oral; bolus; t=0.000000 min; d=60 mg\n"\
                           "RepeatedBolus; via=Oral; repeated_bolus; t=0:24:120 h; d=60 mg")

        form.addParam('t0', params.FloatParam, label="Initial time (see help)", default=0,
                      help="Same units as input experiment")
        form.addParam('tF', params.FloatParam, label="Final time (see help)", default=24*7,
                      help="Same units as input experiment")
        form.addParam('deltaT', params.FloatParam, label="Time step (see help)", default=0.5, expertLevel=LEVEL_ADVANCED,
                      help="Same units as input experiment")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulate')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addSample(self, sampleName, doseList, simulationsX, y):
        newSample = PKPDSample()
        newSample.sampleName = sampleName
        newSample.variableDictPtr = self.outputExperiment.variables
        newSample.doseDictPtr = self.outputExperiment.doses
        newSample.descriptors = {}
        newSample.doseList = doseList
        if type(self.varNameY)!=list:
            newSample.addMeasurementPattern([self.varNameY])
            newSample.addMeasurementColumn("t", simulationsX)
            newSample.addMeasurementColumn(self.varNameY,y)
        else:
            for j in range(len(self.varNameY)):
                newSample.addMeasurementPattern([self.varNameY[j]])
            newSample.addMeasurementColumn("t", simulationsX)
            for j in range(len(self.varNameY)):
                newSample.addMeasurementColumn(self.varNameY[j], y[j])
        self.outputExperiment.samples[sampleName] = newSample

    def runSimulate(self):
        # Take first experiment
        self.protODE = self.inputODEs[0].get()

        if hasattr(self.protODE, "outputExperiment"):
            self.experiment = self.readExperiment(self.protODE.outputExperiment.fnPKPD, show=False)
        elif hasattr(self.protODE, "outputExperiment1"):
            self.experiment = self.readExperiment(self.protODE.outputExperiment1.fnPKPD, show=False)
        else:
            raise Exception("Cannot find an outputExperiment in the input ODE")

        if hasattr(self.protODE, "outputFitting"):
            self.fitting = self.readFitting(self.protODE.outputFitting.fnFitting, show=False)
        elif hasattr(self.protODE, "outputFitting1"):
            self.fitting = self.readFitting(self.protODE.outputFitting1.fnFitting, show=False)

        self.varNameX = self.fitting.predictor.varName
        if type(self.fitting.predicted)!=list:
            self.varNameY = self.fitting.predicted.varName
        else:
            self.varNameY = [var.varName for var in self.fitting.predicted]

        # Create output object
        self.outputExperiment = PKPDExperiment()
        tvar = PKPDVariable()
        tvar.varName = "t"
        tvar.varType = PKPDVariable.TYPE_NUMERIC
        tvar.role = PKPDVariable.ROLE_TIME
        tvar.units = createUnit(self.experiment.getTimeUnits().unit)
        Nsamples = int(math.ceil((self.tF.get() - self.t0.get()) / self.deltaT.get())) + 1
        if tvar.units == PKPDUnit.UNIT_TIME_MIN:
            Nsamples *= 60

        self.outputExperiment.variables[self.varNameX] = tvar
        if type(self.fitting.predicted) != list:
            self.outputExperiment.variables[self.varNameY] = self.experiment.variables[self.varNameY]
        else:
            for varName in self.varNameY:
                self.outputExperiment.variables[varName] = self.experiment.variables[varName]
        self.outputExperiment.general["title"] = "Simulated ODE response"
        self.outputExperiment.general["comment"] = "Simulated ODE response"
        self.outputExperiment.vias = self.experiment.vias

        # Read the doses
        doseLines = []
        for line in self.doses.get().replace('\n', ';;').split(';;'):
            doseLines.append(line)

        doseIdx = 0
        simulationsY = None
        doseList = []
        for protODEPtr in self.inputODEs:
            self.protODE = protODEPtr.get()

            if hasattr(self.protODE, "outputExperiment"):
                self.experiment = self.readExperiment(self.protODE.outputExperiment.fnPKPD, show=False)
            elif hasattr(self.protODE, "outputExperiment1"):
                self.experiment = self.readExperiment(self.protODE.outputExperiment1.fnPKPD, show=False)
            else:
                raise Exception("Cannot find an outputExperiment in the input ODE")

            if hasattr(self.protODE, "outputFitting"):
                self.fitting = self.readFitting(self.protODE.outputFitting.fnFitting, show=False)
            elif hasattr(self.protODE, "outputFitting1"):
                self.fitting = self.readFitting(self.protODE.outputFitting1.fnFitting, show=False)

            for viaName in self.experiment.vias:
                if not viaName in self.outputExperiment.vias:
                    self.outputExperiment.vias[viaName]=copy.copy(self.experiment.vias[viaName])

            # Create drug source
            self.clearGroupParameters()
            self.createDrugSource()

            # Setup model
            self.model = self.protODE.createModel()
            self.model.setExperiment(self.outputExperiment)
            self.model.deltaT = self.deltaT.get()
            self.model.setXVar(self.varNameX)
            self.model.setYVar(self.varNameY)

            self.model.x = [self.t0.get()+i*self.model.deltaT for i in range(0,Nsamples)]
            self.modelList.append(self.model)

            tokens = doseLines[doseIdx].split(';')
            if len(tokens) < 5:
                print("Skipping dose: ", line)
                continue
            dosename = tokens[0].strip()
            self.outputExperiment.doses[dosename] = PKPDDose()
            self.outputExperiment.doses[dosename].parseTokens(tokens,self.outputExperiment.vias)
            doseList.append(dosename)

            auxSample = PKPDSample()
            auxSample.descriptors = {}
            auxSample.doseDictPtr = self.outputExperiment.doses
            auxSample.variableDictPtr = self.outputExperiment.variables
            auxSample.doseList = [dosename]
            auxSample.interpretDose()
            self.drugSource.setDoses(auxSample.parsedDoseList, self.t0.get()-10, self.tF.get()+10)
            self.model.drugSource = self.drugSource

            # Simulate the different responses
            simulationsX = self.model.x
            self.setTimeRange(None)

            parameters = self.fitting.sampleFits[0].parameters
            print("Input model: %s"%self.protODE.getObjLabel())
            print("Sample name: %s"%self.fitting.sampleFits[0].sampleName)
            print("Parameters: ",parameters)
            print("Dose: %s"%doseLines[doseIdx])
            print(" ")

            # Prepare source and this object
            self.protODE.configureSource(self.drugSource)
            parameterNames = self.getParameterNames() # Necessary to count the number of source and PK parameters

            # Prepare the model
            y = self.forwardModel(parameters, [simulationsX]*self.getResponseDimension())

            # Keep results
            if simulationsY is None:
                simulationsY=copy.copy(y)
            else:
                for j in range(self.getResponseDimension()):
                    simulationsY[j] += y[j]

            doseIdx += 1

        self.addSample("Simulation", doseList, simulationsX, simulationsY[0])

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        for protODEPtr in self.inputODEs:
            self._defineSourceRelation(protODEPtr.get(), self.outputExperiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        return msg

    def _validate(self):
        msg=[]
        return msg
