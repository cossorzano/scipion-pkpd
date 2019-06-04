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
from pkpd.objects import PKPDVariable
from pkpd.pkpd_units import createUnit, PKPDUnit, multiplyUnits, strUnit
import math


class ProtPKPDNCANumeric(ProtPKPD):
    """ Non-compartmental analysis just based on the samples. The results are valid only up to T.\n
        It is valid for any kind of via (intravenous, extravascular, ...)\n.
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'nca numeric'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select the experiment with the measurements you want to analyze')
        form.addParam('xVar', params.StringParam, label="Variable to analyze", default="Cp")

    #--------------------------- STEPS functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runAnalysis',self.inputExperiment.getObjId(),self.xVar.get())
        self._insertFunctionStep('createOutputStep')

    def analyzeSample(self, sample, t, C):
        t=np.insert(t,0,0.0)
        C=np.insert(C,0,0.0)
        if sample.descriptors is None:
            sample.descriptors = {}

        self.AUC0t = 0
        self.AUMC0t = 0
        t0 = t[0]
        tperiod0=0 # Time at which the dose was given
        for idx in range(0,t.shape[0]-1):
            dt = (t[idx+1]-t[idx])
            if C[idx+1]>=C[idx]: # Trapezoidal in the raise
                self.AUC0t  += 0.5*dt*(C[idx]+C[idx+1])
                self.AUMC0t += 0.5*dt*(C[idx]*t[idx]+C[idx+1]*t[idx+1])
            else: # Log-trapezoidal in the decay
                decrement = C[idx]/C[idx+1]
                K = math.log(decrement)
                B = K/dt
                self.AUC0t  += dt*(C[idx]-C[idx+1])/K
                self.AUMC0t += (C[idx]*(t[idx]-tperiod0)-C[idx+1]*(t[idx+1]-tperiod0))/B-(C[idx+1]-C[idx])/(B*B)

            if idx==0:
                self.Cmax=C[idx]
                self.Tmax=t[idx]-t0
            else:
                if C[idx]>self.Cmax:
                    self.Cmax=C[idx]
                    self.Tmax=t[idx]-t0

        self.MRT = self.AUMC0t/self.AUC0t

        print("   Cmax=%f [%s]"%(self.Cmax,strUnit(self.Cunits.unit)))
        print("   Tmax=%f [min]"%self.Tmax)
        print("   AUC0t=%f [%s]"%(self.AUC0t,strUnit(self.AUCunits)))
        print("   AUMC0t=%f [%s]"%(self.AUMC0t,strUnit(self.AUMCunits)))
        print("   MRT=%f [min]"%self.MRT)


        sample.descriptors["Tmax"]=self.Tmax
        sample.descriptors["Cmax"]=self.Cmax
        sample.descriptors["MRT"] = self.MRT
        sample.descriptors["AUC0t"] = self.AUC0t
        sample.descriptors["AUMC0t"] = self.AUMC0t

    def runAnalysis(self,objId,xvarName):
        self.outputExperiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        tvarName = None
        for varName in self.outputExperiment.variables:
            if self.outputExperiment.variables[varName].role == PKPDVariable.ROLE_TIME:
                tvarName = varName
                break
        if tvarName is None:
            raise Exception("Cannot find a time variable in this experiment")

        tvar = PKPDVariable()
        tvar.varName = "t"
        tvar.varType = PKPDVariable.TYPE_NUMERIC
        tvar.role = PKPDVariable.ROLE_TIME
        tvar.units = createUnit("min")

        self.Cunits = self.outputExperiment.variables[xvarName].units
        self.AUCunits = multiplyUnits(tvar.units.unit, self.Cunits.unit)
        self.AUMCunits = multiplyUnits(tvar.units.unit, self.AUCunits)

        AUCvar = PKPDVariable()
        AUCvar.varName = "AUC0t"
        AUCvar.varType = PKPDVariable.TYPE_NUMERIC
        AUCvar.role = PKPDVariable.ROLE_LABEL
        AUCvar.units = createUnit(strUnit(self.AUCunits))

        AUMCvar = PKPDVariable()
        AUMCvar.varName = "AUMC0t"
        AUMCvar.varType = PKPDVariable.TYPE_NUMERIC
        AUMCvar.role = PKPDVariable.ROLE_LABEL
        AUMCvar.units = createUnit(strUnit(self.AUMCunits))

        MRTvar = PKPDVariable()
        MRTvar.varName = "MRT"
        MRTvar.varType = PKPDVariable.TYPE_NUMERIC
        MRTvar.role = PKPDVariable.ROLE_LABEL
        MRTvar.units = createUnit("min")

        Cmaxvar = PKPDVariable()
        Cmaxvar.varName = "Cmax"
        Cmaxvar.varType = PKPDVariable.TYPE_NUMERIC
        Cmaxvar.role = PKPDVariable.ROLE_LABEL
        Cmaxvar.units = createUnit(strUnit(self.Cunits.unit))

        Tmaxvar = PKPDVariable()
        Tmaxvar.varName = "Tmax"
        Tmaxvar.varType = PKPDVariable.TYPE_NUMERIC
        Tmaxvar.role = PKPDVariable.ROLE_LABEL
        Tmaxvar.units = createUnit("min")

        self.outputExperiment.variables["AUC0t"] = AUCvar
        self.outputExperiment.variables["AUMC0t"] = AUMCvar
        self.outputExperiment.variables["MRT"] = MRTvar
        self.outputExperiment.variables["Cmax"] = Cmaxvar
        self.outputExperiment.variables["Tmax"] = Tmaxvar

        for sampleName, sample in self.outputExperiment.samples.iteritems():
            [t,Cp] = sample.getXYValues(tvarName,xvarName)
            self.analyzeSample(sample, t[0], Cp[0])

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        self._defineSourceRelation(self.inputExperiment.get(), self.outputExperiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        msg.append("Non-compartmental analysis for the observations of the variable %s"%self.xVar.get())
        return msg