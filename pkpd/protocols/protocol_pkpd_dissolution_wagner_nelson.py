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

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.pkpd_units import createUnit
from pkpd.utils import uniqueFloatValues, calculateAUC0t, smoothPchip

# tested in test_workflow_levyplot.py
# tested in test_workflow_ivivc.py

class ProtPKPDDeconvolutionWagnerNelson(ProtPKPD):
    """ Calculate the absorption profile of an in vivo concentration profile using
        the Wagner-Nelson approach. This is only valid for profiles that have been
        modelled with a monocompartment PK model.

        The formula is Fabs(t)=(Cp(t)+Ke*AUC0t(t))/(Ke*AUC0inf)
        where Ke=Cl/V

        In this implementation it is assumed that AUC0inf is the last AUC0t observed,
        meaning that Cp(t) has almost vanished in the last samples"""

    _label = 'deconvolution Wagner Nelson'

    SAME_INPUT = 0
    ANOTHER_INPUT = 1

    BIOAVAIL_NONE = 0
    BIOAVAIL_MULT = 1
    BIOAVAIL_DIV = 2

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="In-vivo profiles",
                      pointerClass='PKPDExperiment', help='Make sure that it has a clearance parameter (Cl) and central volume (V)')
        form.addParam('externalIV', params.EnumParam, choices=['Get impulse response from the same input fit',
                                                               'Get impulse response from another fit'],
                      label='Impulse response source',
                      default=self.SAME_INPUT, help="The impulse response is an estimate of the intravenous response")
        form.addParam('externalIVODE', params.PointerParam, label="External impulse response ODE model",
                      condition='externalIV==1',
                      pointerClass='ProtPKPDMonoCompartment, ProtPKPDTwoCompartments,ProtPKPDODERefine,' \
                                   ' ProtPKPDTwoCompartmentsClint, ProtPKPDTwoCompartmentsClintCl',
                      help='Select a run of an ODE model. It should be ideally the intravenous response.')
        form.addParam('timeVar', params.StringParam, label="Time variable", default="t",
                      help='Which variable contains the time stamps.')
        form.addParam('concVar', params.StringParam, label="Concentration variable", default="Cp",
                      help='Which variable contains the plasma concentration.')
        form.addParam('normalize', params.BooleanParam, label="Normalize by dose", default=True,
                      help='Normalize the output by AUC0inf, so that a total absorption is represented by 100.')
        form.addParam('saturate', params.BooleanParam, label="Saturate at 100%", default=True, condition='normalize',
                      help='Saturate the absorption so that there cannot be values beyond 100')
        form.addParam('considerBioaval', params.EnumParam, label="Consider bioavailability", default=self.BIOAVAIL_NONE,
                      choices=['Do not correct','Multiply deconvolution by bioavailability','Divide deconvolution by bioavailability'],
                      help='Take into account the bioavailability')
        form.addParam('resampleT', params.FloatParam, label="Resample profiles (time step)", default=-1,
                      help='Resample the input profiles at this time step (make sure it is in the same units as the input). '
                           'Leave it to -1 for no resampling')
        form.addParam('smooth', params.BooleanParam, label="Monotonic smooth", default=True,
                      help='Apply a Pchip interpolation to make sure that the Adissolved is monotonically increasing')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('deconvolve',self.inputExperiment.get().getObjId())
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

    def estimateAUCrightTail(self, t, Cp, Cl, V):
        idx = Cp>0
        x = t[idx]
        y = np.log10(Cp[idx])
        ydiff = np.diff(y)
        idx = -1
        while ydiff[idx]<0:
            idx-=1
        if idx<-3:
            p = np.polyfit(x[idx:], y[idx:], 1)
            Ke = Cl/V
            lastCt=np.power(10,p[0]*x[-1]+p[1])
            return lastCt/np.abs(Ke)
        else:
            return 0

    def deconvolve(self, objId1):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        # Create output object
        self.outputExperiment = PKPDExperiment()
        tvar = PKPDVariable()
        tvar.varName = "t"
        tvar.varType = PKPDVariable.TYPE_NUMERIC
        tvar.role = PKPDVariable.ROLE_TIME
        tvar.units = createUnit(self.experiment.getTimeUnits().unit)

        Avar = PKPDVariable()
        Avar.varName = "A"
        Avar.varType = PKPDVariable.TYPE_NUMERIC
        Avar.role = PKPDVariable.ROLE_MEASUREMENT
        Avar.units = createUnit("none")

        self.outputExperiment.variables[tvar.varName] = tvar
        self.outputExperiment.variables[Avar.varName] = Avar
        self.outputExperiment.general["title"]="Deconvolution of the amount released"
        self.outputExperiment.general["comment"]="Amount released at any time t"

        # Get the input sample from another experiment if necessary
        sampleFrom = None
        if self.externalIV.get() == self.ANOTHER_INPUT:
            anotherExperiment = self.readExperiment(self.externalIVODE.get().outputExperiment.fnPKPD)
            for _, sampleFrom in anotherExperiment.samples.items(): # Take the first sample from the reference
                break

        timeRange = self.experiment.getRange(self.timeVar.get())
        for sampleName, sample in self.experiment.samples.items():
            # Get t, Cp
            t=np.asarray(sample.getValues(self.timeVar.get()),dtype=np.float64)
            Cp=np.asarray(sample.getValues(self.concVar.get()),dtype=np.float64)
            Cp=np.clip(Cp,0.0,None)
            t=np.insert(t,0,0) # Add (0,0) to the profile
            Cp=np.insert(Cp,0,0)
            t, Cp = uniqueFloatValues(t, Cp)
            if self.resampleT.get()>0:
                B = InterpolatedUnivariateSpline(t, Cp, k=1)
                t = np.arange(np.min(t),np.max(t)+self.resampleT.get(),self.resampleT.get())
                Cp = B(t)

            # Calculate AUC0t
            AUC0t=calculateAUC0t(t,Cp)

            # Deconvolve
            if self.externalIV.get()==self.SAME_INPUT:
                sampleFrom = sample
            Cl=float(sampleFrom.descriptors['Cl'])
            V=float(sampleFrom.descriptors['V'])
            Ke=Cl/V
            AUC0inf = float(AUC0t[-1])+self.estimateAUCrightTail(t,Cp, Cl, V)
            A = (Cp + Ke * AUC0t) / Ke
            if self.normalize.get():
                A *= 100/AUC0inf
                if self.saturate.get():
                    A = np.clip(A,None,100.0)
            A = np.clip(A,0,None)
            if self.smooth:
                if t[0]>0:
                    t = np.insert(t, 0, 0)
                    A = np.insert(A, 0, 0)
                if self.saturate.get() and self.normalize.get():
                    A = np.clip(A, None, 100.0)
                A = np.clip(smoothPchip(t, A),0,None)
                if self.saturate.get() and self.normalize.get():
                    A = np.clip(A, None, 100.0)

            if self.considerBioaval.get()==self.BIOAVAIL_DIV:
                A /= sample.getBioavailability()
            elif self.considerBioaval.get()==self.BIOAVAIL_MULT:
                A *= sample.getBioavailability()
            self.addSample(sampleName,t,A)

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        self._defineSourceRelation(self.inputExperiment.get(), self.outputExperiment)

    def _validate(self):
        retval = []
        if self.externalIV.get() == self.ANOTHER_INPUT:
            sample = self.readExperiment(self.externalIVODE.get().outputExperiment.fnPKPD).getFirstSample()
        else:
            sample = self.readExperiment(self.inputExperiment.get().fnPKPD).getFirstSample()
        if sample is None:
            reval.append('Cannot find a sample in the input experiment')
        else:
            if sample.getDescriptorValue('Cl') is None:
                retval.append('Cannot find Cl in the input experiment')
            if sample.getDescriptorValue('V') is None:
                retval.append('Cannot find V in the input experiment')
        return retval

    def _summary(self):
        retval = []
        return retval
