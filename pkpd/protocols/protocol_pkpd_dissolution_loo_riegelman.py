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
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.pkpd_units import createUnit
from pkpd.utils import uniqueFloatValues, calculateAUC0t, smoothPchip

# tested in test_workflow_levyplot.py

class ProtPKPDDeconvolutionLooRiegelman(ProtPKPD):
    """ Calculate the absorption profile of an in vivo concentration profile using
        the Loo-Riegelman approach. This is only valid for profiles that have been
        modelled with a two-compartments PK model.

        The formula is Fabs(t)=(Cp(t)+Cperipheral(t)+K10*AUC0t(t))/(K10*AUC0inf)
        where K10=Cl/V and
        Cperipheral(t_n)=k12*Delta Cp*Delta t/2+k12/k21 * Cp(t_n-1)(1-exp(-k21*Delta t)+Cperipheral(t_n-1)*exp(-k21*Delta t)

        In this implementation it is assumed that AUC0inf is the last AUC0t observed,
        meaning that Cp(t) has almost vanished in the last samples.

        This protocol is much more accurate when the input Cp(t) is reinterpolated to a small time step like 0.5 minutes.

        Reference: Leon Shargel, Susanna Wu-Pong, Andrew B.C. Yu. Applied Biopharmaceutics & Pharmacokinetics, 6e. McGraw Hill, 1999. Chap. 7"""

    _label = 'deconvolution Loo-Riegelman'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="In-vivo profiles",
                      pointerClass='PKPDExperiment', help='Make sure that it has a clearance parameter (Cl) and central volume (V)')
        form.addParam('timeVar', params.StringParam, label="Time variable", default="t",
                      help='Which variable contains the time stamps.')
        form.addParam('concVar', params.StringParam, label="Concentration variable", default="Cp",
                      help='Which variable contains the plasma concentration.')
        form.addParam('resampleT', params.FloatParam, label="Resample profiles (time step)", default=0.5,
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

    def calculateCperipheral(self,t,Cp,k12,k21):
        Cperipheral=np.zeros(t.shape)
        for n in range(1,Cperipheral.size):
            DeltaCp=np.abs(Cp[n]-Cp[n-1])
            DeltaT=t[n]-t[n-1]
            Cperipheral[n] = k12*DeltaCp*DeltaT*0.5 + \
                             k12 / k21 * Cp[n-1]*(1 - np.exp(-k21 * DeltaT)) + \
                             Cperipheral[n-1] * np.exp(-k21 * DeltaT)
            if Cperipheral[n]<0.0:
                Cperipheral[n]=0.0
        return Cperipheral

    def deconvolve(self, objId1):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

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
        Avar.units = createUnit("none")

        self.outputExperiment.variables[tvar.varName] = tvar
        self.outputExperiment.variables[Avar.varName] = Avar
        self.outputExperiment.general["title"]="Deconvolution of the amount released"
        self.outputExperiment.general["comment"]="Amount released at any time t"

        timeRange = self.experiment.getRange(self.timeVar.get())
        for sampleName, sample in self.experiment.samples.iteritems():
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
            AUC0inf = float(AUC0t[-1])

            # Calculate peripheral
            Cl=float(sample.descriptors['Cl'])
            V=float(sample.descriptors['V'])
            Clp=float(sample.descriptors['Clp'])
            Vp=float(sample.descriptors['Vp'])
            k12 = Clp / V
            k21 = Clp / Vp
            Cperipheral = self.calculateCperipheral(t,Cp,k12,k21)

            # Deconvolve
            k10=Cl/V
            A = (Cp + Cperipheral + k10 * AUC0t) / (k10 * AUC0inf) * 100
            A = np.clip(A,0,100)
            if self.smooth:
                if t[0]>0:
                    t = np.insert(t, 0, 0)
                    A = np.insert(A, 0, 0)
                A = np.clip(smoothPchip(t, A),0,100)

            self.addSample(sampleName,t,A)

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        self._defineSourceRelation(self.inputExperiment.get(), self.outputExperiment)

    def _validate(self):
        return []

    def _summary(self):
        retval = []
        return retval
