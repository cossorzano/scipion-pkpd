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
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.pkpd_units import createUnit
from .protocol_pkpd_ode_base import ProtPKPDODEBase
from pkpd.biopharmaceutics import DrugSource
from pkpd.utils import twoWayUniqueFloatValues

# Tested in test_workflow_deconvolution
# Tested by test_workflow_levyplot
# Tested in test_workflow_deconvolution2


class ProtPKPDDeconvolve(ProtPKPDODEBase):
    """ Deconvolve the drug dissolution from a compartmental model."""

    _label = 'dissol deconv'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputODE', params.PointerParam, label="Input ODE model",
                      pointerClass='ProtPKPDMonoCompartment, ProtPKPDTwoCompartments,ProtPKPDODERefine',
                      help='Select a run of an ODE model')
        form.addParam('normalize', params.BooleanParam, label="Normalize by dose", default=True,
                      help='Normalize the output by the input dose, so that a total absorption is represented by 100.')
        form.addParam('removeTlag', params.BooleanParam, label="Remove tlag effect", default=True,
                      help='If set to True, then the deconvolution is performed ignoring the the tlag in the absorption.'
                           'This homogeneizes the different responses.')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('deconvolve',self.inputODE.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addSample(self, sampleName, t, y):
        newSample = PKPDSample()
        newSample.sampleName = sampleName
        newSample.variableDictPtr = self.outputExperiment.variables
        newSample.descriptors = {}
        newSample.addMeasurementPattern(["A"])
        tUnique, yUnique = twoWayUniqueFloatValues(t,y)
        newSample.addMeasurementColumn("t", tUnique)
        newSample.addMeasurementColumn("A",yUnique)
        self.outputExperiment.samples[sampleName] = newSample

    def deconvolve(self, objId):
        self.protODE = self.inputODE.get()
        self.experiment = self.readExperiment(self.protODE.outputExperiment.fnPKPD)
        self.fitting = self.readFitting(self.protODE.outputFitting.fnFitting)
        self.varNameX = self.fitting.predictor.varName
        self.varNameY = self.fitting.predicted.varName

        # Create drug source
        self.clearGroupParameters()
        self.createDrugSource()

        # Create output object
        self.outputExperiment = PKPDExperiment()
        tvar = PKPDVariable()
        tvar.varName = "t"
        tvar.varType = PKPDVariable.TYPE_NUMERIC
        tvar.role = PKPDVariable.ROLE_TIME
        tvar.units = createUnit("min")

        Avar = PKPDVariable()
        Avar.varName = "A"
        Avar.varType = PKPDVariable.TYPE_NUMERIC
        Avar.role = PKPDVariable.ROLE_MEASUREMENT
        if self.normalize.get():
            Avar.units = createUnit("none")
        else:
            Avar.units = createUnit(self.experiment.getDoseUnits())

        self.outputExperiment.variables[tvar.varName] = tvar
        self.outputExperiment.variables[Avar.varName] = Avar
        self.outputExperiment.general["title"]="Deconvolution of the amount released"
        self.outputExperiment.general["comment"]="Amount released at any time t"

        # Simulate the different responses
        timeRange = self.experiment.getRange(self.varNameX)
        deltaT = 0.5
        t = np.arange(0.0,timeRange[1],deltaT)
        for sampleName, sample in self.experiment.samples.iteritems():
            self.printSection("Deconvolving "+sampleName)
            sample.interpretDose()
            drugSource = DrugSource()
            drugSource.setDoses(sample.parsedDoseList, 0.0, timeRange[1])

            p=[]
            tlag=0
            for paramName in drugSource.getParameterNames():
                p.append(float(sample.getDescriptorValue(paramName)))
                if paramName.endswith('_tlag') and self.removeTlag.get():
                    tlag=float(sample.getDescriptorValue(paramName))
            drugSource.setParameters(p)

            cumulatedDose=0.0
            A=t*0.0 # Allocate memory
            totalReleased = drugSource.getAmountReleasedUpTo(10*t[-1])/100 # Divided by 100 to have a number between 0 and 100
            print("t(min) A(%s)"%Avar.units._toString())
            for i in range(t.size):
                cumulatedDose+=drugSource.getAmountReleasedAt(t[i],deltaT)
                A[i]=cumulatedDose
                if self.normalize.get():
                    A[i] /= totalReleased
                print("%f %f"%(t[i],A[i]))
                # print("%f %f %f %f"%(t[i], A[i], drugSource.getAmountReleasedAt(t[i], 0.5), drugSource.getAmountReleasedUpTo(t[i] + 0.5)))
            self.addSample(sampleName,t-tlag,A)

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        self._defineSourceRelation(self.inputODE.get(), self.outputExperiment)

    def _validate(self):
        return []

    def _summary(self):
        return []