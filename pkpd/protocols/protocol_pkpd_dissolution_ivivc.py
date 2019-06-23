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
from scipy.optimize import differential_evolution

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
                      pointerClass='ProtPKPDDeconvolve,ProtPKPDDeconvolutionWagnerNelson', help='Select an experiment with dissolution profiles')
        form.addParam('timeScale', params.EnumParam, label="Time scaling",
                      choices=["None (Fabs(t)=Adissol(t))",
                               "t0 (Fabs(t)=Adissol(t-t0)",
                               "Linear scale (Fabs(t)=Adissol(k*t))",
                               "Affine transformation (Fabs(t)=Adissol(k*(t-t0))"], default=0,
                      help='Fabs is the fraction absorbed in vivo, while Adissol is the amount dissolved in vitro.')
        form.addParam('t0Bounds',params.StringParam,label='Bounds t0',default='[-100,100]',
                      condition='timeScale==1 or timeScale==3',
                      help='Make sure it is in the same time units as the inputs')
        form.addParam('kBounds',params.StringParam,label='Bounds k',default='[0.1,10]',
                      condition='timeScale==2 or timeScale==3')
        form.addParam('responseScale', params.EnumParam, label="Response scaling",
                      choices=["Linear scale (Fabs(t)=A*Adissol(t))",
                               "Affine transformation (Fabs(t)=A*Adissol(t)+B"], default=0,
                      help='Fabs is the fraction absorbed in vivo, while Adissol is the amount dissolved in vitro.')
        form.addParam('ABounds',params.StringParam,label='Bounds A',default='[0.1,10]')
        form.addParam('BBounds',params.StringParam,label='Bounds B',default='[-50,50]',
                      condition='responseScale==1')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllIvIvC',self.inputInVitro.get().getObjId(),self.inputInVivo.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addSample(self, sampleName, individualFrom, vesselFrom, Fabs, Adissol,optimum,R):
        newSample = PKPDSample()
        newSample.sampleName = sampleName
        newSample.variableDictPtr = self.outputExperiment.variables
        newSample.descriptors = {}
        newSample.addMeasurementPattern(["Adissol"])
        newSample.addMeasurementColumn("Adissol", Adissol)
        newSample.addMeasurementColumn("Fabs",Fabs)
        self.outputExperiment.samples[sampleName] = newSample
        if self.timeScale.get()==1:
            timeScaleMsg="tvitro=tvivo-t0"
        elif self.timeScale.get()==2:
            timeScaleMsg = "tvitro=k*tvivo"
        elif self.timeScale.get()==3:
            timeScaleMsg = "tvitro=k*(tvivo-t0)"
        if self.responseScale.get()==0:
            responseMsg="Fabs(t)=A*Adissol(t)"
        elif self.responseScale.get()==1:
            responseMsg="Fabs(t)=A*Adissol(t)+B"

        t0=None
        k=None
        A=None
        B=None
        i=0
        for prm in self.parameters:
            exec ("%s=%f" % (prm, optimum[i]))
            i+=1

        self.outputExperiment.addLabelToSample(sampleName, "from", "individual---vesel", "%s---%s"%(individualFrom,vesselFrom))
        if not t0 is None:
            self.outputExperiment.addParameterToSample(sampleName, "t0", PKPDUnit.UNIT_TIME_MIN, timeScaleMsg, t0)
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

    def goalFunction(self,x):
        t0=0.0
        k=1.0
        A=1.0
        B=0.0
        i=0
        try:
            for prm in self.parameters:
                exec ("%s=%f" % (prm, x[i]))
                i+=1
            tvitroUnique=np.clip(k*(self.tvivoUnique-t0),self.tvitroMin,self.tvitroMax)
            self.AdissolReinterpolatedUnique = np.clip(A*self.B(tvitroUnique)+B,0.0,None)
            self.residuals = self.AdissolReinterpolatedUnique-self.FabsUnique
            error=np.sqrt(np.mean(self.residuals**2))
            if error<self.bestError:
                print("New minimum error=%f"%error,"x=%s"%np.array2string(x,max_line_width=1000))
                self.bestError=error
            if self.verbose:
                for i in range(len(self.AdissolReinterpolatedUnique)):
                    print("i=",i,"tvitro[i]=",tvitroUnique[i],"tvivo[i]=",self.tvivoUnique[i],"AdissolReinterpolated[i]=",self.AdissolReinterpolatedUnique[i],"Fabs[i]",self.FabsUnique[i])
        except:
            return 1e38
        return error

    def parseBounds(self,bounds):
        tokens=bounds.strip()[1:-1].split(',')
        return np.asarray(tokens,dtype=np.float64)

    def produceAdissol(self,parameterInVitro,tmax):
        tvitro = np.arange(0,tmax+1,1)
        self.protFit.model.x = tvitro
        Avitro = self.protFit.model.forwardModel(parameterInVitro)[0]
        return (tvitro, Avitro)

    def calculateAllIvIvC(self, objId1, objId2):
        parametersInVitro, vesselNames=self.getInVitroModels()
        profilesInVivo, sampleNames=self.getInVivoProfiles()

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
        allt0=[]
        allk=[]
        allA=[]
        allB=[]
        allR=[]
        self.parameters=[]
        self.bounds=[]
        if self.timeScale.get() == 1:  # tvitro=tvivo-t0
            self.parameters.append('t0')
            self.bounds.append(self.parseBounds(self.t0Bounds.get()))
        elif self.timeScale.get() == 2:  # tvitro=k*tvivo
            self.parameters.append('k')
            self.bounds.append(self.parseBounds(self.kBounds.get()))
        elif self.timeScale.get() == 3:  # tvitro=k*(tvivo-t0)
            self.parameters.append('k')
            self.parameters.append('t0')
            self.bounds.append(self.parseBounds(self.kBounds.get()))
            self.bounds.append(self.parseBounds(self.t0Bounds.get()))
        self.parameters.append('A')
        self.bounds.append(self.parseBounds(self.ABounds.get()))
        if self.responseScale.get() == 1:  # Affine
            self.parameters.append('B')
            self.bounds.append(self.parseBounds(self.BBounds.get()))

        invitroIdx=0
        self.verbose=False
        for parameterInVitro in parametersInVitro:
            invivoIdx=0
            for self.tvivo,self.Fabs in profilesInVivo:
                print("New combination %d"%i)
                self.FabsUnique, self.tvivoUnique = uniqueFloatValues(self.Fabs, self.tvivo)
                self.tvitro, self.Adissol=self.produceAdissol(parameterInVitro,np.max(self.tvivoUnique*10))
                self.AdissolUnique, self.tvitroUnique = uniqueFloatValues(self.Adissol, self.tvitro)
                self.tvitroMin=np.min(self.tvitroUnique)
                self.tvitroMax=np.max(self.tvitroUnique)
                # for i in range(len(self.tvivoUnique)):
                #    print("i=",i,"tvivo[i]=",self.tvivoUnique[i],"Fabs[i]",self.FabsUnique[i])
                # for i in range(len(self.tvitroUnique)):
                #    print("i=",i,"tvitro[i]=",self.tvitroUnique[i],"Adissol[i]",self.AdissolUnique[i])
                self.B = InterpolatedUnivariateSpline(self.tvitroUnique, self.AdissolUnique, k=1)

                self.bestError = 1e38
                optimum = differential_evolution(self.goalFunction,self.bounds,popsize=50)
                # self.verbose=True
                self.goalFunction(optimum.x)

                j = 0
                t0 = 0.0
                k = 1.0
                A = 1.0
                B = 0.0
                for prm in self.parameters:
                    exec ("%s=%f" % (prm, optimum.x[j]))
                    exec ("all%s.append(%f)" % (prm, optimum.x[j]))
                    j += 1

                # Evaluate correlation
                R2=np.clip(1-np.var(self.residuals)/np.var(self.Fabs),0.0,1.0)
                R=sqrt(R2)
                allR.append(R)

                self.addSample("ivivc_%04d"%i,sampleNames[invivoIdx],vesselNames[invitroIdx],self.AdissolReinterpolatedUnique,self.FabsUnique,optimum.x,R)
                i+=1
                invivoIdx+=1
            invitroIdx+=1

        fh=open(self._getPath("summary.txt"),"w")
        self.summarize(fh,allt0,"t0")
        self.summarize(fh,allk,"k")
        self.summarize(fh,allA,"A")
        self.summarize(fh,allB,"B")
        self.doublePrint(fh," ")
        self.summarize(fh,allR,"Correlation coefficient (R)")
        self.doublePrint(fh," ")

        if self.timeScale.get() == 0:
            timeStr = "t"
        elif self.timeScale.get() == 1:
            timeStr = "t-t0"
        elif self.timeScale.get() == 2:
            timeStr = "k*t"
        elif self.timeScale.get() == 3:
            timeStr = "k*(t-t0)"
        if self.responseScale.get() == 0:
            eqStr = "Fabs(%s)=A*Adissol(t)" % timeStr
        elif self.responseScale.get() == 1:
            eqStr = "Fabs(%s)=A*Adissol(t)+B" % timeStr
        self.doublePrint(fh,"IVIVC equation: %s"%eqStr)
        fh.close()

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def _validate(self):
        return []

    def _summary(self):
        retval = []
        if self.timeScale.get()==0:
            retval.append("Time scaling: none")
        elif self.timeScale.get()==1:
            retval.append("Time scaling: t0 (tvivo=tvitro-t0)")
        elif self.timeScale.get()==2:
            retval.append("Time scaling: linear transformation (tvivo=k*tvitro)")
        elif self.timeScale.get()==3:
            retval.append("Time scaling: affine transformation (tvivo=k*(tvitro-t0))")
        if self.responseScale.get()==0:
            retval.append("Response scaling: Fabs(t)=A*Adissol(t)")
        elif self.responseScale.get()==1:
            retval.append("Response scaling: Fabs(t)=A*Adissol(t)+B")
        self.addFileContentToMessage(retval,self._getPath("summary.txt"))
        return retval
