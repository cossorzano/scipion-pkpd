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
from pkpd.utils import uniqueFloatValues, parseOperation
from pkpd.pkpd_units import createUnit, PKPDUnit


# tested in test_workflow_levyplot

from .protocol_pkpd_dissolution_levyplot import ProtPKPDDissolutionLevyPlot

class ProtPKPDDissolutionIVIVCGeneric(ProtPKPDDissolutionLevyPlot):
    """ Calculate the in vitro-in vivo correlation between two experiments. Each experiment may have
        several profiles and all vs all profiles are calculated. You may scale the time between the two
        sets of experiments"""

    _label = 'dissol ivivc generic'

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
        form.addParam('timeScale', params.StringParam, label="Time scaling", default="$(t)",
                      help='Write a time scaling function like $(t), or $[a0]+$[a1]*$(t)+$[a2]*np.power($(t),2). '
                           'The meaning is tvitro=a0+a1*tvivo+a2*tvivo^2 in the last example. '
                           'You have all numpy functions at your disposal: np.sqrt, np.log10, np.log, np.exp')
        form.addParam('timeBounds',params.StringParam,label='Time coefficients bounds',default='',
                      help='a0: [-100,100]; a1: [0.0001,100]; a2: [1e-6,1e-4]')
        form.addParam('responseScale', params.StringParam, label="Response scaling", default="$(Adissol)",
                      help='Write a response scaling function like $[k1]*$(Adissol)+$[k2]*np.pow($(Adissol),2)$. '
                           'Make sure that the time and response scaling coefficients have different names. '
                           'The meaning is Fabs=k1*Adissol+k2*Adissol^2 in the last example. '
                           'You have all numpy functions at your disposal: np.sqrt, np.log10, np.log, np.exp')
        form.addParam('responseBounds',params.StringParam,label='Response coefficients bounds',default='',
                      help='k1: [0.0001,100]; k2: [1e-6,1e-4]')

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
        self.outputExperiment.addLabelToSample(sampleName, "from", "individual---vesel", "%s---%s"%(individualFrom,vesselFrom))

        i=0
        timeScaleMsg="tvitro=%s"%self.timeScale.get().replace('$(t)','$(tvivo)')
        for prm in self.coeffTimeList:
            self.outputExperiment.addParameterToSample(sampleName, prm, PKPDUnit.UNIT_NONE, timeScaleMsg, optimum[i])
            i+=1
        responseMsg="Fabs(t)=%s"%self.responseScale.get()
        for prm in self.coeffResponseList:
            self.outputExperiment.addParameterToSample(sampleName, prm, PKPDUnit.UNIT_NONE, responseMsg, optimum[i])
            i+=1

        self.outputExperiment.addParameterToSample(sampleName, "R", PKPDUnit.UNIT_NONE, "IVIV Correlation coefficient", R)

    def summarize(self,fh,x,msg):
        if len(x)>0:
            p = np.percentile(x,[2.5, 50, 97.5],axis=0)
            self.doublePrint(fh,"%s: median=%f; 95%% Confidence interval=[%f,%f]"%(msg,p[1],p[0],p[2]))

    def goalFunction(self,x):
        i=0
        try:
            for prm in self.parameters:
                exec ("%s=%f" % (prm, x[i]))
                i+=1

            t=self.tvivoUnique
            tvitroUnique=eval(self.parsedTimeOperation)
            tvitroUnique=np.clip(tvitroUnique,self.tvitroMin,self.tvitroMax)
            self.AdissolPredicted = self.BAdissol(tvitroUnique)
            Adissol = self.AdissolPredicted
            FabsPredicted = eval(self.parsedResponseOperation)
            self.FabsPredicted = np.clip(FabsPredicted,0.0,None)
            self.residualsForward = self.FabsPredicted-self.FabsUnique

            errorForward  = np.sum(self.residualsForward**2)
            self.residuals = self.residualsForward
            error = errorForward
            if error<self.bestError:
                print("New minimum error=%f"%error,"x=%s"%np.array2string(x,max_line_width=1000))
                self.bestError=error
            if self.verbose:
                print("Forward error")
                for i in range(len(tvitroUnique)):
                    print("i=",i,"tvitro[i]=",tvitroUnique[i],"tvivo[i]=",self.tvivoUnique[i],"AdissolPredicted[i]=",self.AdissolPredicted[i],"Fabs[i]",self.FabsPredicted[i])
                print("Error forward",np.sqrt(np.mean(self.residualsForward**2)),np.sum(self.residualsForward**2))
        except:
           return 1e38
        return error

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

    def constructBounds(self, coeffList):
        boundsDict = {}
        def parseBounds(allBoundsStr):
            try:
                for prmBound in allBoundsStr.split(';'):
                    tokens = prmBound.split(':')
                    coeffName = tokens[0].strip()
                    boundStr = tokens[1].strip()
                    values = boundStr.strip().split(',')
                    boundsDict[coeffName]=(float(values[0][1:]), float(values[1][:-1]))
            except:
                return

        parseBounds(self.timeBounds.get())
        parseBounds(self.responseBounds.get())
        if len(boundsDict)!=len(coeffList):
            raise Exception("The number of bounds (%d) is different from the number of parameters (%d)"%\
                            (len(boundsDict),len(coeffList)))

        boundList=[]
        for coeffName in coeffList:
            if coeffName in boundsDict:
                boundList.append(boundsDict[coeffName])
            else:
                raise Exception("Cannot find the bound of %s"%coeffName)

        return boundList

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
        self.outputExperiment.general["title"] = "In-vitro In-vivo correlation generic"
        self.outputExperiment.general["comment"] = "Time in vivo vs time in vitro"

        self.parsedTimeOperation, self.varTimeList, self.coeffTimeList = parseOperation(self.timeScale.get())
        self.parsedResponseOperation, self.varResponseList, self.coeffResponseList = parseOperation(self.responseScale.get())

        self.parameters = self.coeffTimeList + self.coeffResponseList
        self.bounds = self.constructBounds(self.parameters)

        i=1
        allCoeffs = []
        allR=[]

        invitroIdx=0
        self.verbose=False
        for parameterInVitro in parametersInVitro:
            invivoIdx=0
            for self.tvivo,self.Fabs in profilesInVivo:
                print("New combination %d"%i)
                self.FabsUnique, self.tvivoUnique = uniqueFloatValues(self.Fabs, self.tvivo)
                self.tvitro, self.Adissol=self.produceAdissol(parameterInVitro,np.max(self.tvivoUnique*10))
                self.AdissolUnique, self.tvitroUnique = uniqueFloatValues(self.Adissol, self.tvitro)
                self.tvivoMin=np.min(self.tvivoUnique)
                self.tvivoMax=np.max(self.tvivoUnique)
                self.tvitroMin=np.min(self.tvitroUnique)
                self.tvitroMax=np.max(self.tvitroUnique)

                # Make sure they are sorted in x
                self.tvivoUnique, self.FabsUnique = uniqueFloatValues(self.tvivoUnique, self.FabsUnique)
                self.tvitroUnique, self.AdissolUnique = uniqueFloatValues(self.tvitroUnique, self.AdissolUnique)

                self.BAdissol = InterpolatedUnivariateSpline(self.tvitroUnique, self.AdissolUnique, k=1)

                self.bestError = 1e38
                optimum = differential_evolution(self.goalFunction,self.bounds,popsize=50)
                # self.verbose=True
                allCoeffs.append(optimum.x.tolist())


                # Evaluate correlation
                self.goalFunction(optimum.x)
                R2=np.clip(1-np.var(self.residuals)/np.var(self.Fabs),0.0,1.0)
                R=sqrt(R2)
                allR.append(R)

                self.addSample("ivivc_%04d" % i, sampleNames[invivoIdx], vesselNames[invitroIdx], self.FabsUnique, self.AdissolPredicted, optimum.x, R)
                i+=1
                invivoIdx+=1
            invitroIdx+=1

        fh=open(self._getPath("summary.txt"),"w")
        allCoeffs = np.asarray(allCoeffs)
        i=0
        for prm in self.parameters:
            self.summarize(fh,allCoeffs[:,i],prm)
            i+=1
        self.doublePrint(fh," ")
        self.summarize(fh,allR,"Correlation coefficient (R)")
        self.doublePrint(fh," ")

        self.doublePrint(fh,"Time scale: tvitro=%s"%self.timeScale.get().replace('$(t)','$(tvivo)'))
        self.doublePrint(fh,"Response scale: Fabs(t)=%s"%self.responseScale.get())
        fh.close()

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def _validate(self):
        return []

    def _summary(self):
        retval = []
        self.addFileContentToMessage(retval,self._getPath("summary.txt"))
        return retval
