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
from pkpd.pkpd_units import createUnit, divideUnits
from pkpd.utils import uniqueFloatValues, calculateAUC0t, smoothPchip

class ProtPKPDDeconvolutionLooRiegelmanInverse(ProtPKPD):
    """ Given a profile of amount absorbed, find the central and peripheral concentrations that gave raise to it.
        This is an inverse Loo-Riegelman problem.

        Reference: Humbert, H., Cabiac, M.D., Bosshardt, H. In vitro-in vivo correlation of a modified-release oral
                   form of ketotigen: in vitro dissolution rate specification. J. Pharm Sci 1994, 83, 131-136
    """

    _label = 'inverse Loo-Riegelman'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="In-vivo profiles",
                      pointerClass='PKPDExperiment', help='It should have the amount absorbed')
        form.addParam('timeVar', params.StringParam, label="Time variable", default="t",
                      help='Which variable contains the time stamps.')
        form.addParam('amountVar', params.StringParam, label="Absorbed amount variable", default="A",
                      help='Which variable contains the amount absorbed.')
        form.addParam('resampleT', params.FloatParam, label="Resample profiles (time step)", default=0.5,
                      help='Resample the input profiles at this time step (make sure it is in the same units as the input). '
                           'Leave it to -1 for no resampling')
        form.addParam('k10', params.FloatParam, label="Elimination rate (k10)",
                      help='Units t^-1')
        form.addParam('k12', params.FloatParam, label="Rate from central to peripheral (k12)",
                      help='Units t^-1')
        form.addParam('k21', params.FloatParam, label="Rate from peripheral to central (k21)",
                      help='Units t^-1')
        form.addParam('V', params.FloatParam, label="Central volume (V) [L]")
        form.addParam('Vp', params.FloatParam, label="Peripheral volume (Vp) [L]")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('solve',self.inputExperiment.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addSample(self, sampleName, t, C, Cp):
        newSample = PKPDSample()
        newSample.sampleName = sampleName
        newSample.variableDictPtr = self.outputExperiment.variables
        newSample.descriptors = {}
        newSample.addMeasurementPattern(["C","Cp"])
        newSample.addMeasurementColumn("t", t)
        newSample.addMeasurementColumn("C", C)
        newSample.addMeasurementColumn("Cp", Cp)
        self.outputExperiment.samples[sampleName] = newSample

    def calculateConcentrations(self,t,A):
        k10 = float(self.k10.get())
        k12 = float(self.k12.get())
        k21 = float(self.k21.get())
        V = float(self.V.get())
        Vp = float(self.Vp.get())

        C = np.zeros(t.shape)
        Cp = np.zeros(t.shape)
        AUC0t = np.zeros(t.shape)
        for n in range(1,C.size):
            DeltaA=np.abs(A[n]-A[n-1])
            DeltaT=t[n]-t[n-1]
            DeltaT2=DeltaT/2
            K1 = np.exp(-k21 * DeltaT)
            K2 = k10/(1+k12*DeltaT2+k10*DeltaT2)

            C[n] = DeltaA/V - Vp*Cp[n-1]/V*K1 + C[n-1]*k12*DeltaT2-C[n-1]*k12/k21*(1-K1)\
                   -C[n-1]*k10*DeltaT2-AUC0t[n-1]*K2

            DeltaC = C[n]-C[n-1]

            # Update AUC0t
            if DeltaC>0:  # Trapezoidal in the raise
                AUC0t[n] = AUC0t[n-1] + DeltaT/2 * (C[n - 1] + C[n])
            else:  # Log-trapezoidal in the decay
                if C[n-1] > 0 and C[n] > 0:
                    decrement = C[n-1] / C[n]
                    K = math.log(decrement)
                    AUC0t[n] = AUC0t[n - 1] + DeltaT * (C[n-1] - C[n]) / K
                else:
                    AUC0t[n] = AUC0t[n - 1]

            Cp[n] = Cp[n-1]*K1+k12/k21*C[n-1]*(1-K1)+k12*DeltaC*DeltaT2

            if C[n]<0.0:
                C[n]=0.0
            if Cp[n]<0.0:
                Cp[n]=0.0
        return C, Cp

    def solve(self, objId1):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        # Create output object
        self.outputExperimet = None

        timeRange = self.experiment.getRange(self.timeVar.get())
        for sampleName, sample in self.experiment.samples.iteritems():
            # Get t, A
            t=np.asarray(sample.getValues(self.timeVar.get()),dtype=np.float64)
            A=np.asarray(sample.getValues(self.amountVar.get()),dtype=np.float64)
            if t[0]>0:
                t=np.insert(t,0,0) # Add (0,0) to the profile
                A=np.insert(A,0,0)
            t, A = uniqueFloatValues(t, A)
            if self.resampleT.get()>0:
                B = InterpolatedUnivariateSpline(t, A, k=1)
                t = np.arange(np.min(t),np.max(t)+self.resampleT.get(),self.resampleT.get())
                A = B(t)

            # Find C and Cp
            C, Cp = self.calculateConcentrations(t,A)
            print(C, Cp)

            if self.outputExperimet is None:
                self.outputExperiment = PKPDExperiment()
                tvar = PKPDVariable()
                tvar.varName = "t"
                tvar.varType = PKPDVariable.TYPE_NUMERIC
                tvar.role = PKPDVariable.ROLE_TIME
                tvar.units = createUnit(self.experiment.getTimeUnits().unit)

                Cvar = PKPDVariable()
                Cvar.varName = "C"
                Cvar.varType = PKPDVariable.TYPE_NUMERIC
                Cvar.role = PKPDVariable.ROLE_MEASUREMENT
                Lunits = createUnit("L")
                Cvar.units = createUnit(divideUnits(self.experiment.getVarUnits(self.amountVar.get()),Lunits.unit))
                Cvar.comment = "Concentration central compartment"

                Cpvar = PKPDVariable()
                Cpvar.varName = "Cp"
                Cpvar.varType = PKPDVariable.TYPE_NUMERIC
                Cpvar.role = PKPDVariable.ROLE_MEASUREMENT
                Cpvar.units = createUnit(divideUnits(self.experiment.getVarUnits(self.amountVar.get()),Lunits.unit))
                Cpvar.comment = "Concentration peripheral compartment"

                self.outputExperiment.variables[tvar.varName] = tvar
                self.outputExperiment.variables[Cvar.varName] = Cvar
                self.outputExperiment.variables[Cpvar.varName] = Cpvar
                self.outputExperiment.general["title"]="Inverse Loo-Riegelman"
                self.outputExperiment.general["comment"]=""

            self.addSample(sampleName,t,C,Cp)

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        self._defineSourceRelation(self.inputExperiment.get(), self.outputExperiment)

    def _validate(self):
        return []

    def _summary(self):
        retval = []
        return retval
