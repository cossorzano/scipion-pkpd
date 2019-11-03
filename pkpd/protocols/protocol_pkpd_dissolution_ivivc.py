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
from math import sqrt
import numpy as np
import sys
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import differential_evolution

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.utils import uniqueFloatValues, computeXYmean, smoothPchip
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
                      pointerClass='ProtPKPDDeconvolve,ProtPKPDDeconvolutionWagnerNelson,ProtPKPDDeconvolutionLooRiegelman, PKPDExperiment',
                      help='Select an experiment with dissolution profiles')
        form.addParam('removeInVitroTlag',params.BooleanParam, label="Remove in vitro tlag", default=True,
                      help="If there is an in vitro tlag, set it to 0")
        form.addParam('timeScale', params.EnumParam, label="Time scaling",
                      choices=["None (Fabs(t)=Adissol(t))",
                               "t0 (Fabs(t)=Adissol(t-t0)",
                               "Linear scale (Fabs(t)=Adissol(k*t))",
                               "Affine transformation (Fabs(t)=Adissol(k*(t-t0))",
                               "Power scale (Fabs(t)=Adissol(k*t^alpha)",
                               "Delayed power scale (Fabs(t)=Adissol(k*(t-t0)^alpha)"
                               ], default=0,
                      help='Fabs is the fraction absorbed in vivo, while Adissol is the amount dissolved in vitro.')
        form.addParam('t0Bounds',params.StringParam,label='Bounds t0',default='[-100,100]',
                      condition='timeScale==1 or timeScale==3 or timeScale==5',
                      help='Make sure it is in the same time units as the inputs')
        form.addParam('kBounds',params.StringParam,label='Bounds k',default='[0.1,10]',
                      condition='timeScale==2 or timeScale==3 or timeScale==4 or timeScale==5')
        form.addParam('alphaBounds',params.StringParam,label='Bounds alpha',default='[0.1,10]',
                      condition='timeScale==4 or timeScale==5')
        form.addParam('responseScale', params.EnumParam, label="Response scaling",
                      choices=["None",
                               "Linear scale (Fabs(t)=A*Adissol(t))",
                               "Affine transformation (Fabs(t)=A*Adissol(t)+B"], default=0,
                      help='Fabs is the fraction absorbed in vivo, while Adissol is the amount dissolved in vitro.')
        form.addParam('ABounds',params.StringParam,label='Bounds A',default='[0.8,1.2]', condition='responseScale>=1')
        form.addParam('BBounds',params.StringParam,label='Bounds B',default='[-50,50]',
                      condition='responseScale==2')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllIvIvC',self.inputInVitro.get().getObjId(),self.inputInVivo.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addParametersToExperiment(self, outputExperiment, sampleName, individualFrom, vesselFrom, optimum, R):
        if self.timeScale.get()==1:
            timeScaleMsg="tvitro=tvivo-t0"
        elif self.timeScale.get()==2:
            timeScaleMsg = "tvitro=k*tvivo"
        elif self.timeScale.get()==3:
            timeScaleMsg = "tvitro=k*(tvivo-t0)"
        elif self.timeScale.get()==4:
            timeScaleMsg = "tvitro=k*tvivo^alpha"
        elif self.timeScale.get()==5:
            timeScaleMsg = "tvitro=k*(tvivo-t0)^alpha"
        if self.responseScale.get()==1:
            responseMsg="Fabs(t)=A*Adissol(t)"
        elif self.responseScale.get()==2:
            responseMsg="Fabs(t)=A*Adissol(t)+B"

        t0=None
        k=None
        alpha=None
        A=None
        B=None
        i=0
        for prm in self.parameters:
            exec ("%s=%f" % (prm, optimum[i]))
            i+=1

        outputExperiment.addLabelToSample(sampleName, "from", "individual---vesel", "%s---%s"%(individualFrom,vesselFrom))
        if not t0 is None:
            outputExperiment.addParameterToSample(sampleName, "t0", PKPDUnit.UNIT_TIME_MIN, timeScaleMsg, t0)
        if not k is None:
            outputExperiment.addParameterToSample(sampleName, "k", PKPDUnit.UNIT_NONE, timeScaleMsg, k)
        if not alpha is None:
            outputExperiment.addParameterToSample(sampleName, "alpha", PKPDUnit.UNIT_NONE, timeScaleMsg, alpha)
        if not A is None:
            outputExperiment.addParameterToSample(sampleName, "A", PKPDUnit.UNIT_NONE, responseMsg, A)
        if not B is None:
            outputExperiment.addParameterToSample(sampleName, "B", PKPDUnit.UNIT_NONE, responseMsg, B)
        outputExperiment.addParameterToSample(sampleName, "R", PKPDUnit.UNIT_NONE, "IVIV Correlation coefficient", R)
        outputExperiment.addLabelToSample(sampleName, "from", "individual---vesel", "%s---%s"%(individualFrom,vesselFrom))

    def addSample(self, sampleName, individualFrom, vesselFrom, optimum, R, set=1):
        newSampleFabs = PKPDSample()
        newSampleFabs.sampleName = sampleName
        newSampleFabs.variableDictPtr = self.outputExperimentFabs.variables
        newSampleFabs.descriptors = {}
        newSampleFabs.addMeasurementColumn("tvitroReinterpolated", self.tvitroReinterpolated)
        newSampleFabs.addMeasurementColumn("AdissolReinterpolated", self.AdissolReinterpolated)
        newSampleFabs.addMeasurementColumn("tvivo", self.tvivoUnique)
        newSampleFabs.addMeasurementColumn("FabsPredicted", self.FabsPredicted)
        newSampleFabs.addMeasurementColumn("Fabs",self.FabsUnique)
        if set==1:
            self.outputExperimentFabs.samples[sampleName] = newSampleFabs
            self.addParametersToExperiment(self.outputExperimentFabs, sampleName, individualFrom, vesselFrom, optimum, R)
        else:
            self.outputExperimentFabsSingle.samples[sampleName] = newSampleFabs
            self.addParametersToExperiment(self.outputExperimentFabsSingle, sampleName, individualFrom, vesselFrom, optimum, R)

        if set==1:
            newSampleAdissol = PKPDSample()
            newSampleAdissol.sampleName = sampleName
            newSampleAdissol.variableDictPtr = self.outputExperimentAdissol.variables
            newSampleAdissol.descriptors = {}
            newSampleAdissol.addMeasurementColumn("tvivoReinterpolated", self.tvivoReinterpolated)
            newSampleAdissol.addMeasurementColumn("FabsReinterpolated", self.FabsReinterpolated)
            newSampleAdissol.addMeasurementColumn("tvitro", self.tvitroUnique)
            newSampleAdissol.addMeasurementColumn("AdissolPredicted", self.AdissolPredicted)
            newSampleAdissol.addMeasurementColumn("Adissol",self.AdissolUnique)
            self.outputExperimentAdissol.samples[sampleName] = newSampleAdissol
            self.addParametersToExperiment(self.outputExperimentAdissol, sampleName, individualFrom, vesselFrom, optimum, R)

    def summarize(self,fh,x,msg):
        if len(x)>0:
            p = np.percentile(x,[2.5, 50, 97.5],axis=0)
            self.doublePrint(fh,"%s: median=%f; 95%% Confidence interval=[%f,%f]"%(msg,p[1],p[0],p[2]))

    def guaranteeMonotonicity(self):
            self.FabsPredicted = np.clip(smoothPchip(self.FabsUnique, self.FabsPredicted),0,100)
            self.AdissolPredicted = np.clip(smoothPchip(self.AdissolUnique, self.AdissolPredicted),0,100)

    def calculateIndividualError(self, x, tvitroUnique, tvivoUnique):
        self.guaranteeMonotonicity()
        self.residualsForward = self.FabsPredicted - self.FabsUnique
        self.residualsBackward = self.AdissolUnique - self.AdissolPredicted

        idx = np.isnan(self.residualsForward)
        self.residualsForward[idx] = self.FabsUnique[idx]
        idx = np.isnan(self.residualsBackward)
        self.residualsBackward[idx] = self.AdissolUnique[idx]

        errorForward = np.mean(self.residualsForward ** 2)
        errorBackward = np.mean(self.residualsBackward ** 2)
        if errorForward > errorBackward:
            self.residuals = self.residualsForward
        else:
            self.residuals = self.residualsBackward
        diff=errorBackward-errorForward
        error = 0.5*(errorBackward+errorForward)#+np.abs(diff)

        # error=np.sqrt(np.mean(self.residuals**2))
        # error=np.sqrt(np.sum(self.residuals**2))
        return error, errorBackward, errorForward

    def calculateError(self, x, tvitroUnique, tvivoUnique):
        error, errorBackward, errorForward=self.calculateIndividualError(x, tvitroUnique, tvivoUnique)
        if error < self.bestError:
            print("New minimum error=%f (back=%f, forw=%f) R=%f" % (error,errorBackward,errorForward,self.calculateR()),
                  "x=%s" % np.array2string(x, max_line_width=1000))
            self.bestError = error
            sys.stdout.flush()
        if self.verbose:
            print("Forward error")
            for i in range(len(self.FabsPredicted)):
                print("i=", i, "tvitro[i]=", tvitroUnique[i], "tvivo[i]=", self.tvivoUnique[i], "FabsPredicted[i]=",
                      self.FabsPredicted[i], "Fabs[i]", self.FabsUnique[i])
            print("Backward error")
            for i in range(len(self.AdissolPredicted)):
                print("i=", i, "tvitro[i]=", self.tvitroUnique[i], "tvivo[i]=", tvivoUnique[i], "Adissol[i]=",
                      self.AdissolUnique[i], "AdissolPredicted[i]", self.AdissolPredicted[i])
            print("Error forward", np.sqrt(np.mean(self.residualsForward ** 2)), np.sum(self.residualsForward ** 2))
            print("Error backward", np.sqrt(np.mean(self.residualsBackward ** 2)), np.sum(self.residualsBackward ** 2))
        return error

    def goalFunction(self,x):
        t0=0.0
        k=1.0
        alpha=1.0
        A=1.0
        B=0.0
        i=0
        try:
            for prm in self.parameters:
                exec ("%s=%f" % (prm, x[i]))
                i+=1

            self.tvitroReinterpolated=np.clip(k*np.power(np.clip(self.tvivoUnique-t0,0,None),alpha),self.tvitroMin,self.tvitroMax)
            self.AdissolReinterpolated = self.BAdissol(self.tvitroReinterpolated)
            self.FabsPredicted = np.clip(A*self.AdissolReinterpolated+B,0.0,100)

            self.tvivoReinterpolated=np.clip(np.power(self.tvitroUnique/k,1.0/alpha)+t0,self.tvivoMin,self.tvivoMax)
            self.FabsReinterpolated = self.BFabs(self.tvivoReinterpolated)
            self.AdissolPredicted = np.clip((self.FabsReinterpolated-B)/A,0.0,100)

            error = self.calculateError(x, self.tvitroReinterpolated, self.tvivoReinterpolated)
        except:
            return 1e38
        return error

    def parseBounds(self,bounds):
        tokens=bounds.strip()[1:-1].split(',')
        return np.asarray(tokens,dtype=np.float64)

    def produceAdissol(self,parameterInVitro,tmax):
        tvitro = np.arange(0,tmax+1,1)
        if self.removeInVitroTlag:
            i=0
            for prmName in self.protFit.model.getParameterNames():
                if "tlag" in prmName:
                    parameterInVitro[i]=0.0
                i+=1
        self.protFit.model.x = tvitro
        Avitro = self.protFit.model.forwardModel(parameterInVitro)[0]
        return (tvitro, Avitro)

    def createOutputExperiments(self, set):
        tvitroVar = PKPDVariable()
        tvitroVar.varName = "tvitro"
        tvitroVar.varType = PKPDVariable.TYPE_NUMERIC
        tvitroVar.role = PKPDVariable.ROLE_TIME
        tvitroVar.units = createUnit("min")
        tvitroVar.comment = "tvitro"

        tvivoVar = PKPDVariable()
        tvivoVar.varName = "tvivo"
        tvivoVar.varType = PKPDVariable.TYPE_NUMERIC
        tvivoVar.role = PKPDVariable.ROLE_TIME
        tvivoVar.units = createUnit("min")
        tvivoVar.comment = "tvivo"

        tvitroReinterpolatedVar = PKPDVariable()
        tvitroReinterpolatedVar.varName = "tvitroReinterpolated"
        tvitroReinterpolatedVar.varType = PKPDVariable.TYPE_NUMERIC
        tvitroReinterpolatedVar.role = PKPDVariable.ROLE_TIME
        tvitroReinterpolatedVar.units = createUnit("min")
        tvitroReinterpolatedVar.comment = "tvitro reinterpolated"

        tvivoReinterpolatedVar = PKPDVariable()
        tvivoReinterpolatedVar.varName = "tvivoReinterpolated"
        tvivoReinterpolatedVar.varType = PKPDVariable.TYPE_NUMERIC
        tvivoReinterpolatedVar.role = PKPDVariable.ROLE_TIME
        tvivoReinterpolatedVar.units = createUnit("min")
        tvivoReinterpolatedVar.comment = "tvivo reinterpolated"

        AdissolVar = PKPDVariable()
        AdissolVar.varName = "Adissol"
        AdissolVar.varType = PKPDVariable.TYPE_NUMERIC
        AdissolVar.role = PKPDVariable.ROLE_MEASUREMENT
        AdissolVar.units = createUnit(self.experimentInVitro.getVarUnits(self.varNameY))
        AdissolVar.comment = "Amount disolved in vitro"

        FabsVar = PKPDVariable()
        FabsVar.varName = "Fabs"
        FabsVar.varType = PKPDVariable.TYPE_NUMERIC
        FabsVar.role = PKPDVariable.ROLE_MEASUREMENT
        FabsVar.units = createUnit(self.experimentInVivo.getVarUnits("A"))
        FabsVar.comment = "Amount absorbed in vivo"

        AdissolReinterpolatedVar = PKPDVariable()
        AdissolReinterpolatedVar.varName = "AdissolReinterpolated"
        AdissolReinterpolatedVar.varType = PKPDVariable.TYPE_NUMERIC
        AdissolReinterpolatedVar.role = PKPDVariable.ROLE_MEASUREMENT
        AdissolReinterpolatedVar.units = createUnit(self.experimentInVitro.getVarUnits(self.varNameY))
        AdissolReinterpolatedVar.comment = "Time reinterpolated amount disolved in vitro"

        FabsReinterpolatedVar = PKPDVariable()
        FabsReinterpolatedVar.varName = "FabsReinterpolated"
        FabsReinterpolatedVar.varType = PKPDVariable.TYPE_NUMERIC
        FabsReinterpolatedVar.role = PKPDVariable.ROLE_MEASUREMENT
        FabsReinterpolatedVar.units = createUnit(self.experimentInVivo.getVarUnits("A"))
        FabsReinterpolatedVar.comment = "Time reinterpolated amount absorbed in vivo"

        AdissolPredictedVar = PKPDVariable()
        AdissolPredictedVar.varName = "AdissolPredicted"
        AdissolPredictedVar.varType = PKPDVariable.TYPE_NUMERIC
        AdissolPredictedVar.role = PKPDVariable.ROLE_MEASUREMENT
        AdissolPredictedVar.units = createUnit(self.experimentInVitro.getVarUnits(self.varNameY))
        AdissolPredictedVar.comment = "Predicted amount disolved in vitro"

        FabsPredictedVar = PKPDVariable()
        FabsPredictedVar.varName = "FabsPredicted"
        FabsPredictedVar.varType = PKPDVariable.TYPE_NUMERIC
        FabsPredictedVar.role = PKPDVariable.ROLE_MEASUREMENT
        FabsPredictedVar.units = createUnit(self.experimentInVivo.getVarUnits("A"))
        FabsPredictedVar.comment = "Predicted amount absorbed in vivo"

        if set==1:
            self.outputExperimentFabs = PKPDExperiment()
            self.outputExperimentFabs.variables[tvitroReinterpolatedVar.varName] = tvitroReinterpolatedVar
            self.outputExperimentFabs.variables[AdissolReinterpolatedVar.varName] = AdissolReinterpolatedVar
            self.outputExperimentFabs.variables[tvivoVar.varName] = tvivoVar
            self.outputExperimentFabs.variables[FabsVar.varName] = FabsVar
            self.outputExperimentFabs.variables[FabsPredictedVar.varName] = FabsPredictedVar
            self.outputExperimentFabs.general["title"] = "In-vitro In-vivo correlation"
            self.outputExperimentFabs.general["comment"] = "Fabs vs Predicted Fabs"

            self.outputExperimentAdissol = PKPDExperiment()
            self.outputExperimentAdissol.variables[tvivoReinterpolatedVar.varName] = tvivoReinterpolatedVar
            self.outputExperimentAdissol.variables[FabsReinterpolatedVar.varName] = FabsReinterpolatedVar
            self.outputExperimentAdissol.variables[tvitroVar.varName] = tvitroVar
            self.outputExperimentAdissol.variables[AdissolVar.varName] = AdissolVar
            self.outputExperimentAdissol.variables[AdissolPredictedVar.varName] = AdissolPredictedVar
            self.outputExperimentAdissol.general["title"] = "In-vitro In-vivo correlation"
            self.outputExperimentAdissol.general["comment"] = "Adissol vs Predicted Adissol"
        else:
            self.outputExperimentFabsSingle = PKPDExperiment()
            self.outputExperimentFabsSingle.variables[tvitroReinterpolatedVar.varName] = tvitroReinterpolatedVar
            self.outputExperimentFabsSingle.variables[AdissolReinterpolatedVar.varName] = AdissolReinterpolatedVar
            self.outputExperimentFabsSingle.variables[tvivoVar.varName] = tvivoVar
            self.outputExperimentFabsSingle.variables[FabsVar.varName] = FabsVar
            self.outputExperimentFabsSingle.variables[FabsPredictedVar.varName] = FabsPredictedVar
            self.outputExperimentFabsSingle.general["title"] = "In-vitro In-vivo correlation"
            self.outputExperimentFabsSingle.general["comment"] = "Fabs vs Predicted Fabs"


    def calculateR(self):
        R2forward = np.clip(1 - np.var(self.residualsForward) / np.var(self.FabsUnique), 0.0, 1.0)
        R2backward = np.clip(1 - np.var(self.residualsBackward) / np.var(self.AdissolUnique), 0.0, 1.0)
        R2 = max(R2forward, R2backward)
        R = sqrt(R2)
        return R

    def calculateAllIvIvC(self, objId1, objId2):
        self.parametersInVitro, self.vesselNames=self.getInVitroModels()
        self.profilesInVivo, self.sampleNames=self.getInVivoProfiles()

        self.createOutputExperiments(set=1)

        i=1
        allt0=[]
        allk=[]
        allalpha=[]
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
        elif self.timeScale.get() == 4:  # tvitro=k*tvivo^alpha
            self.parameters.append('k')
            self.parameters.append('alpha')
            self.bounds.append(self.parseBounds(self.kBounds.get()))
            self.bounds.append(self.parseBounds(self.alphaBounds.get()))
        elif self.timeScale.get() == 5:  # tvitro=k*(tvivo-t0)^alpha
            self.parameters.append('k')
            self.parameters.append('alpha')
            self.parameters.append('t0')
            self.bounds.append(self.parseBounds(self.kBounds.get()))
            self.bounds.append(self.parseBounds(self.alphaBounds.get()))
            self.bounds.append(self.parseBounds(self.t0Bounds.get()))
        if self.responseScale.get() >= 1: # Linear
            self.parameters.append('A')
            self.bounds.append(self.parseBounds(self.ABounds.get()))
        if self.responseScale.get() == 2: # Affine
            self.parameters.append('B')
            self.bounds.append(self.parseBounds(self.BBounds.get()))

        # Compute all pairs
        invitroIdx=0
        self.verbose=False
        vitroList = []
        vivoList = []
        for parameterInVitro, vesselName in izip(self.parametersInVitro,self.vesselNames):
            invivoIdx=0
            if "tvitroMax" in self.experimentInVitro.variables:
                tvitroMax=float(self.experimentInVitro.samples[vesselName].getDescriptorValue("tvitroMax"))
            else:
                tvitroMax=1e38
            for self.tvivo,self.Fabs in self.profilesInVivo:
                print("New combination %d"%i)
                self.FabsUnique, self.tvivoUnique = uniqueFloatValues(self.Fabs, self.tvivo)
                self.tvitro, self.Adissol=self.produceAdissol(parameterInVitro,min(np.max(self.tvivoUnique*10),tvitroMax))
                self.AdissolUnique, self.tvitroUnique = uniqueFloatValues(self.Adissol, self.tvitro)
                self.tvivoMin=np.min(self.tvivoUnique)
                self.tvivoMax=np.max(self.tvivoUnique)
                self.tvitroMin=np.min(self.tvitroUnique)
                self.tvitroMax=np.max(self.tvitroUnique)
                #for i in range(len(self.tvivoUnique)):
                #   print("i=",i,"tvivo[i]=",self.tvivoUnique[i],"Fabs[i]",self.FabsUnique[i])
                #for i in range(len(self.tvitroUnique)):
                #   print("i=",i,"tvitro[i]=",self.tvitroUnique[i],"Adissol[i]",self.AdissolUnique[i])

                # Make sure they are sorted in x
                self.tvivoUnique, self.FabsUnique = uniqueFloatValues(self.tvivoUnique, self.FabsUnique)
                self.tvitroUnique, self.AdissolUnique = uniqueFloatValues(self.tvitroUnique, self.AdissolUnique)
                vivoList.append((self.tvivoUnique, self.FabsUnique))
                vitroList.append((self.tvitroUnique, self.AdissolUnique))

                self.BAdissol = InterpolatedUnivariateSpline(self.tvitroUnique, self.AdissolUnique, k=1)
                self.BFabs = InterpolatedUnivariateSpline(self.tvivoUnique, self.FabsUnique, k=1)

                self.bestError = 1e38
                if len(self.bounds)>0:
                    optimum = differential_evolution(self.goalFunction,self.bounds,popsize=50)
                    x = optimum.x
                else:
                    x = None
                # self.verbose=True
                self.goalFunction(x)

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
                R = self.calculateR()
                allR.append(R)

                self.addSample("ivivc_%04d" % i, self.sampleNames[invivoIdx], self.vesselNames[invitroIdx], x, R, set=1)
                i+=1
                invivoIdx+=1
            invitroIdx+=1

        fh=open(self._getPath("summary.txt"),"w")
        self.summarize(fh,allt0,"t0")
        self.summarize(fh,allk,"k")
        self.summarize(fh,allalpha,"alpha")
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
        elif self.timeScale.get() == 4:
            timeStr = "k*t^alpha"
        elif self.timeScale.get() == 5:
            timeStr = "k*(t-t0)^alpha"
        if self.responseScale.get() == 0:
            eqStr = "Fabs(%s)=Adissol(t)" % timeStr
        elif self.responseScale.get() == 1:
            eqStr = "Fabs(%s)=A*Adissol(t)" % timeStr
        elif self.responseScale.get() == 2:
            eqStr = "Fabs(%s)=A*Adissol(t)+B" % timeStr
        self.doublePrint(fh,"IVIVC equation: %s"%eqStr)
        fh.close()

        self.outputExperimentFabs.write(self._getPath("experimentFabs.pkpd"))
        self.outputExperimentAdissol.write(self._getPath("experimentAdissol.pkpd"))

        # Compute single
        print("Single IVIVC")
        self.tvivoUnique, self.FabsUnique = computeXYmean(vivoList)
        self.tvitroUnique, self.AdissolUnique = computeXYmean(vitroList)
        self.tvivoUnique, self.FabsUnique = uniqueFloatValues(self.tvivoUnique, self.FabsUnique)
        self.tvitroUnique, self.AdissolUnique = uniqueFloatValues(self.tvitroUnique, self.AdissolUnique)

        self.BAdissol = InterpolatedUnivariateSpline(self.tvitroUnique, self.AdissolUnique, k=1)
        self.BFabs = InterpolatedUnivariateSpline(self.tvivoUnique, self.FabsUnique, k=1)
        self.bestError = 1e38
        if len(self.bounds) > 0:
            optimum = differential_evolution(self.goalFunction, self.bounds, popsize=50)
            x = optimum.x
        else:
            x = None
        self.goalFunction(x)
        R = self.calculateR()
        self.createOutputExperiments(set=2)
        self.addSample("ivivc_single", "AvgVivo", "AvgVitro", x, R, set=2)
        self.outputExperimentFabsSingle.write(self._getPath("experimentFabsSingle.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperimentFabs=self.outputExperimentFabs)
        self._defineSourceRelation(self.inputInVitro.get(), self.outputExperimentFabs)
        self._defineSourceRelation(self.inputInVivo.get(), self.outputExperimentFabs)

        self._defineOutputs(outputExperimentAdissol=self.outputExperimentAdissol)
        self._defineSourceRelation(self.inputInVitro.get(), self.outputExperimentAdissol)
        self._defineSourceRelation(self.inputInVivo.get(), self.outputExperimentAdissol)

        self._defineOutputs(outputExperimentFabsSingle=self.outputExperimentFabsSingle)
        self._defineSourceRelation(self.inputInVitro.get(), self.outputExperimentFabsSingle)
        self._defineSourceRelation(self.inputInVivo.get(), self.outputExperimentFabsSingle)

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
        elif self.timeScale.get()==4:
            retval.append("Time scaling: power transformation (tvivo=k*tvitro^alpha)")
        elif self.timeScale.get()==5:
            retval.append("Time scaling: delayed power transformation (tvivo=k*(tvitro-t0)^alpha)")
        if self.responseScale.get()==0:
            retval.append("Response scaling: Fabs(t)=A*Adissol(t)")
        elif self.responseScale.get()==1:
            retval.append("Response scaling: Fabs(t)=A*Adissol(t)+B")
        self.addFileContentToMessage(retval,self._getPath("summary.txt"))
        return retval
