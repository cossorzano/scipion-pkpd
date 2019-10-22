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
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.fftpack import fft, ifft

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.pkpd_units import createUnit
from pkpd.utils import uniqueFloatValues
from .protocol_pkpd_ode_base import ProtPKPDODEBase
from pkpd.biopharmaceutics import DrugSource, createDeltaDose, createVia

# Tested in test_workflow_deconvolution.py

class ProtPKPDDeconvolveFourier(ProtPKPDODEBase):
    """ Deconvolve the drug dissolution from a compartmental model. It does the deconvolution in
        Fourier so that it only uses the impulse response of the compartmental model. This impulse
        response only depends on the distribution, metabolism and excretion (DME) part of the ADME
        properties, meaning that it overcomes the limitations of a poor modelling of the raise
        of the concentration. On the other side, it has the disadvantage of considering the noise
        as true fluctuations due to the absorption."""

    _label = 'dissol deconv Fourier'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputODE', params.PointerParam, label="Input ODE model",
                      pointerClass='ProtPKPDMonoCompartment, ProtPKPDTwoCompartments, ProtPKPDODERefine',
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
        newSample.addMeasurementColumn("t", t)
        newSample.addMeasurementColumn("A",y)
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

        # Create PK model
        timeRange = self.experiment.getRange(self.varNameX)
        deltaT = 0.5
        t = np.arange(0.0,timeRange[1]*5+deltaT,deltaT)

        model = self.protODE.createModel()
        model.setExperiment(self.experiment)
        model.deltaT = deltaT
        model.setXVar(self.varNameX)
        model.setYVar(self.varNameY)
        model.x = t
        model.t0=0
        model.tF=np.max(t)
        prmNames = model.getParameterNames()
        print(prmNames)

        if len(self.experiment.samples)==0:
            print("Cannot find any sample in the experiment")
            return

        # Create unit dose
        drugSource = DrugSource()
        dose=None
        for sampleName, sample in self.experiment.samples.iteritems():
            sample.interpretDose()
            if len(sample.parsedDoseList)>0:
                dose = createDeltaDose(1.0, via=createVia("Intravenous; iv", self.experiment),
                                       dunits=sample.parsedDoseList[0].dunits)
        if dose is None:
            print("Cannot find any dose among the samples")
            return
        drugSource.setDoses([dose], 0.0, timeRange[1])
        model.drugSource = drugSource

        # Simulate the different responses
        for sampleName, sample in self.experiment.samples.iteritems():
            self.printSection("Deconvolving "+sampleName)
            drugSourceSample = DrugSource()
            drugSourceSample.setDoses(sample.parsedDoseList, 0.0, timeRange[1])

            p=[]
            tlag=0
            for paramName in drugSourceSample.getParameterNames():
                p.append(float(sample.getDescriptorValue(paramName)))
                if paramName.endswith('_tlag') and self.removeTlag.get():
                    tlag=float(sample.getDescriptorValue(paramName))
            drugSourceSample.setParameters(p)

            parameters=[]
            for prmName in prmNames:
                parameters.append(float(sample.descriptors[prmName]))
            model.setParameters(parameters)
            h = model.forwardModel(parameters, [t]*model.getResponseDimension())[0]

            ts,Cs = sample.getXYValues(self.varNameX,self.varNameY)
            ts=np.insert(ts,0,0.0)
            Cs=np.insert(Cs,0,0.0)
            ts, Cs = uniqueFloatValues(ts,Cs)
            B=InterpolatedUnivariateSpline(ts, Cs, k=1)
            C=np.clip(B(t),0.0,None)

            Fabs=np.clip(np.real(ifft(np.divide(fft(C),fft(h)))),0.0,None)

            cumulatedDose=0.0
            A=t*0.0 # Allocate memory
            totalReleased = drugSourceSample.getAmountReleasedUpTo(10*t[-1])/100 # Divided by 100 to have a number between 0 and 100
            print("t(min) A(%s)"%Avar.units._toString())
            for i in range(t.size):
                cumulatedDose+=Fabs[i]
                A[i]=cumulatedDose
                if self.normalize.get():
                    A[i] /= totalReleased
                print("%f %f"%(t[i],A[i]))
                # print("%f %f %f %f"%(t[i], A[i], drugSource.getAmountReleasedAt(t[i], 0.5), drugSource.getAmountReleasedUpTo(t[i] + 0.5)))
            As, ts = uniqueFloatValues(np.clip(A,0,100),t-tlag)
            self.addSample(sampleName,ts,As)

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        self._defineSourceRelation(self.inputODE.get(), self.outputExperiment)

    def _validate(self):
        return []

    def _summary(self):
        return []
