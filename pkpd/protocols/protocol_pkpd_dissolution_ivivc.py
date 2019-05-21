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

from math import sqrt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.utils import uniqueFloatValues
from pkpd.pkpd_units import createUnit, PKPDUnit


# tested in test_workflow_levyplot

from .protocol_pkpd_dissolution_levyplot import ProtPKPDDissolutionLevyPlot

class ProtPKPDDissolutionIVIVC(ProtPKPDDissolutionLevyPlot):
    """ Calculate the in vitro-in vivo correlation between two experiments. Each experiment may have
        several profiles and all vs all profiles are calculated. You may scale the time between the two
        sets of experiments"""

    _label = 'dissol ivivc'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputInVitro', params.PointerParam, label="Dissolution profiles in vitro",
                      pointerClass='ProtPKPDDissolutionFit', help='Select an experiment with dissolution profiles')
        form.addParam('inputInVivo', params.PointerParam, label="Dissolution profiles in vivo",
                      pointerClass='ProtPKPDDeconvolve', help='Select an experiment with dissolution profiles')
        form.addParam('timeScale', params.EnumParam, label="Time scaling",
                      choices=["None (Fabs(t)=Adissol(t))",
                               "tlag (Fabs(t-tlag)=Adissol(t)",
                               "Linear scale (Fabs(k*t)=Adissol(t))",
                               "Affine transformation (Fabs(k*(t-tlag))=Adissol(t)"], default=0,
                      help='Fabs is the fraction absorbed in vivo, while Adissol is the amount dissolved in vitro.')
        form.addParam('responseScale', params.EnumParam, label="Response scaling",
                      choices=["Linear scale (Fabs(t)=A*Adissol(t))",
                               "Affine transformation (Fabs(t)=A*Adissol(t)+B"], default=0,
                      help='Fabs is the fraction absorbed in vivo, while Adissol is the amount dissolved in vitro.')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllIvIvC',self.inputInVitro.get().getObjId(),self.inputInVivo.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addSample(self, sampleName, Fabs, Adissol,tlag,k,A,B,R):
        newSample = PKPDSample()
        newSample.sampleName = sampleName
        newSample.variableDictPtr = self.outputExperiment.variables
        newSample.descriptors = {}
        newSample.addMeasurementPattern(["Adissol"])
        newSample.addMeasurementColumn("Adissol", Adissol)
        newSample.addMeasurementColumn("Fabs",Fabs)
        self.outputExperiment.samples[sampleName] = newSample
        if self.timeScale.get()==1:
            timeScaleMsg="tvivo=tvitro-tlag"
        elif self.timeScale.get()==2:
            timeScaleMsg = "tvivo=k*tvitro"
        elif self.timeScale.get()==3:
            timeScaleMsg = "tvivo=k*(tvitro-tlag)"
        if self.responseScale.get()==0:
            responseMsg="Fabs(t)=A*Adissol(t)"
        elif self.responseScale.get()==1:
            responseMsg="Fabs(t)=A*Adissol(t)+B"
        if not tlag is None:
            self.outputExperiment.addParameterToSample(sampleName, "tlag", PKPDUnit.UNIT_TIME_MIN, timeScaleMsg, tlag)
        if not k is None:
            self.outputExperiment.addParameterToSample(sampleName, "k", PKPDUnit.UNIT_NONE, timeScaleMsg, k)
        if not A is None:
            self.outputExperiment.addParameterToSample(sampleName, "A", PKPDUnit.UNIT_NONE, responseMsg, A)
        if not B is None:
            self.outputExperiment.addParameterToSample(sampleName, "B", PKPDUnit.UNIT_NONE, responseMsg, B)
        self.outputExperiment.addParameterToSample(sampleName, "R", PKPDUnit.UNIT_NONE, "IVIV Correlation coefficient", R)

    def summarize(self,fh,x,msg):
        if len(x)>0:
            p = np.percentile(x,[2.5, 50, 97.5],axis=0)
            self.doublePrint(fh,"%s: median=%f; 95%% Confidence interval=[%f,%f]"%(msg,p[1],p[0],p[2]))

    def calculateAllIvIvC(self, objId1, objId2):
        parametersInVitro=self.getInVitroModels()
        profilesInVivo=self.getInVivoProfiles()

        self.outputExperiment = PKPDExperiment()
        AdissolVar = PKPDVariable()
        AdissolVar.varName = "Adissol"
        AdissolVar.varType = PKPDVariable.TYPE_NUMERIC
        AdissolVar.role = PKPDVariable.ROLE_TIME
        AdissolVar.units = createUnit(self.experimentInVitro.getVarUnits(self.varNameY))
        AdissolVar.comment = "Amount disolved in vitro"

        FabsVar = PKPDVariable()
        FabsVar.varName = "Fabs"
        FabsVar.varType = PKPDVariable.TYPE_NUMERIC
        FabsVar.role = PKPDVariable.ROLE_MEASUREMENT
        FabsVar.units = createUnit(self.experimentInVivo.getVarUnits("A"))
        FabsVar.comment = "Amount absorbed in vivo"

        self.outputExperiment.variables[AdissolVar.varName] = AdissolVar
        self.outputExperiment.variables[FabsVar.varName] = FabsVar
        self.outputExperiment.general["title"] = "In-vitro In-vivo correlation"
        self.outputExperiment.general["comment"] = "Time in vivo vs time in vitro"

        i=1
        allTlag=[]
        allk=[]
        allA=[]
        allB=[]
        allR=[]
        for parameterInVitro in parametersInVitro:
            for t,Fabs in profilesInVivo:
                tvitro, tvivo, Adissol = self.produceLevyPlot(t, parameterInVitro, Fabs)

                # Time scale
                tlag=None
                k=None
                if self.timeScale.get()==1: # tvivo=tvitro-tlag
                    tlag=np.mean(tvitro)-np.mean(tvivo)
                    allTlag.append(tlag)
                    tvivo = tvitro-tlag
                elif self.timeScale.get()==2: # tvivo=k*tvitro
                    k=np.mean(tvivo)/np.mean(tvitro)
                    allk.append(k)
                    tvivo = k*tvitro
                elif self.timeScale.get()==3: # tvivo=k*(tvitro-tlag)
                    p=np.polyfit(tvitro,tvivo,deg=1)
                    k=p[0]
                    tlag=p[1]/k
                    allTlag.append(tlag)
                    allk.append(k)
                    tvivo = k * (tvitro - tlag)

                # Reinterpolate Fabs
                if self.timeScale.get()>0:
                    tvitroUnique, FabsUnique = uniqueFloatValues(tvitro, Fabs)
                    B = InterpolatedUnivariateSpline(tvitroUnique, FabsUnique, k=1)
                    Fabs = B(tvivo)

                # Response scale
                A=None
                B=None
                if self.responseScale.get()==0: # Linear
                    A=np.mean(Fabs)/np.mean(Adissol)
                    allA.append(A)
                    Fabs=A*Adissol
                elif self.responseScale.get()==1: # Affine
                    p=np.polyfit(Adissol,Fabs,deg=1)
                    A=p[0]
                    B=p[1]
                    allA.append(A)
                    allB.append(B)
                    Fabs=A*Adissol+B

                # Evaluate correlation
                residuals = Fabs-Adissol
                R2=1-np.var(residuals)/np.var(Fabs)
                R=sqrt(R2)
                allR.append(R)

                self.addSample("ivivc_%04d"%i,Adissol,Fabs,tlag,k,A,B,R)
                i+=1

        fh=open(self._getPath("summary.txt"),"w")
        self.summarize(fh,allTlag,"tlag")
        self.summarize(fh,allk,"k")
        self.summarize(fh,allA,"A")
        self.summarize(fh,allB,"B")
        self.doublePrint(fh," ")
        self.summarize(fh,allR,"Correlation coefficient (R)")
        fh.close()

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def _validate(self):
        return []

    def _summary(self):
        retval = []
        if self.timeScale.get()==0:
            retval.append("Time scaling: none")
        elif self.timeScale.get()==1:
            retval.append("Time scaling: tlag (tvivo=tvitro-tlag)")
        elif self.timeScale.get()==2:
            retval.append("Time scaling: linear transformation (tvivo=k*tvitro)")
        elif self.timeScale.get()==3:
            retval.append("Time scaling: affine transformation (tvivo=k*(tvitro-tlag))")
        if self.responseScale.get()==0:
            retval.append("Response scaling: Fabs(t)=A*Adissol(t)")
        elif self.responseScale.get()==1:
            retval.append("Response scaling: Fabs(t)=A*Adissol(t)+B")
        self.addFileContentToMessage(retval,self._getPath("summary.txt"))
        return retval
