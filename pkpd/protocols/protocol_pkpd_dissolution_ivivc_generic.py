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
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import differential_evolution

import pyworkflow.protocol.params as params
from pkpd.utils import uniqueFloatValues, parseOperation, computeXYmean
from pkpd.pkpd_units import PKPDUnit


# tested in test_workflow_levyplot

from .protocol_pkpd_dissolution_ivivc import ProtPKPDDissolutionIVIVC

class ProtPKPDDissolutionIVIVCGeneric(ProtPKPDDissolutionIVIVC):
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
                      help='Write a response scaling function like $[k1]*$(Adissol)+$[k2]*np.power($(Adissol),2)$. '
                           'Make sure that the time and response scaling coefficients have different names. '
                           'The meaning is Fabs=k1*Adissol+k2*Adissol^2 in the last example. '
                           'You have all numpy functions at your disposal: np.sqrt, np.log10, np.log, np.exp, np.power')
        form.addParam('responseBounds',params.StringParam,label='Response coefficients bounds',default='',
                      help='k1: [0.0001,100]; k2: [1e-6,1e-4]')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllIvIvC',self.inputInVitro.get().getObjId(),self.inputInVivo.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addParametersToExperiment(self, outputExperiment, sampleName, individualFrom, vesselFrom, optimum, R):
        i=0
        timeScaleMsg="tvitro=%s"%self.timeScale.get().replace('$(t)','$(tvivo)')
        for prm in self.coeffTimeList:
            self.outputExperimentFabs.addParameterToSample(sampleName, prm, PKPDUnit.UNIT_NONE, timeScaleMsg, optimum[i])
            self.outputExperimentAdissol.addParameterToSample(sampleName, prm, PKPDUnit.UNIT_NONE, timeScaleMsg, optimum[i])
            i+=1
        responseMsg="Fabs(t)=%s"%self.responseScale.get()
        for prm in self.coeffResponseList:
            self.outputExperimentFabs.addParameterToSample(sampleName, prm, PKPDUnit.UNIT_NONE, responseMsg, optimum[i])
            self.outputExperimentAdissol.addParameterToSample(sampleName, prm, PKPDUnit.UNIT_NONE, responseMsg, optimum[i])
            i+=1

        self.outputExperimentFabs.addParameterToSample(sampleName, "R", PKPDUnit.UNIT_NONE, "IVIV Correlation coefficient", R)
        self.outputExperimentAdissol.addParameterToSample(sampleName, "R", PKPDUnit.UNIT_NONE, "IVIV Correlation coefficient", R)
        self.outputExperimentFabs.addLabelToSample(sampleName, "from", "individual---vesel", "%s---%s"%(individualFrom,vesselFrom))
        self.outputExperimentAdissol.addLabelToSample(sampleName, "from", "individual---vesel", "%s---%s"%(individualFrom,vesselFrom))

    def goalFunction(self,x):
        i=0
        try:
            for prm in self.parameters:
                exec ("%s=%f" % (prm, x[i]))
                i+=1

            t=self.tvivoUnique
            tvitroUnique=eval(self.parsedTimeOperation)
            self.tvitroReinterpolated=np.clip(tvitroUnique,self.tvitroMin,self.tvitroMax)
            self.AdissolReinterpolated = self.BAdissol(self.tvitroReinterpolated)
            Adissol = self.AdissolReinterpolated
            FabsPredicted = eval(self.parsedResponseOperation)
            self.FabsPredicted = np.clip(FabsPredicted,0.0,100)

            tvitroAux, tvivoAux = uniqueFloatValues(tvitroUnique,self.tvivoUnique)
            Btinv = InterpolatedUnivariateSpline(tvitroAux, tvivoAux, k=1)
            self.tvivoReinterpolated = np.clip(Btinv(self.tvitroUnique), self.tvivoMin, self.tvivoMax)
            self.FabsReinterpolated = self.BFabs(self.tvivoReinterpolated)
            FabsPredictedAux, AdissolAux = uniqueFloatValues(self.FabsPredicted, self.AdissolReinterpolated)
            Bfinv = InterpolatedUnivariateSpline(FabsPredictedAux, AdissolAux, k=1)
            self.AdissolPredicted = np.clip(Bfinv(self.FabsReinterpolated), 0.0, 100)

            error = self.calculateError(x, self.tvitroReinterpolated, self.tvivoReinterpolated)
        except:
           return 1e38
        return error

    def constructBounds(self, coeffList, boundsStr):
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

        parseBounds(boundsStr)
        if len(boundsDict)!=len(coeffList):
            print("Coefficient list: ",coeffList)
            print("Bounds: ",boundsDict)
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
        self.parametersInVitro, self.vesselNames=self.getInVitroModels()
        self.profilesInVivo, self.sampleNames=self.getInVivoProfiles()

        self.createOutputExperiments(set=1)

        self.parsedTimeOperation, self.varTimeList, self.coeffTimeList = parseOperation(self.timeScale.get())
        self.parsedResponseOperation, self.varResponseList, self.coeffResponseList = parseOperation(self.responseScale.get())

        self.parameters = self.coeffTimeList + self.coeffResponseList
        self.bounds = self.constructBounds(self.coeffTimeList, self.timeBounds.get()) + \
                      self.constructBounds(self.coeffResponseList, self.responseBounds.get())
        self.calculateAllIvIvCLoop()

    def calculateAllIvIvCLoop(self):
        i=1
        allCoeffs = []
        allR=[]

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

                # Make sure they are sorted in x
                self.tvivoUnique, self.FabsUnique = uniqueFloatValues(self.tvivoUnique, self.FabsUnique)
                self.tvitroUnique, self.AdissolUnique = uniqueFloatValues(self.tvitroUnique, self.AdissolUnique)
                self.FabsMax = np.max(self.FabsUnique)
                self.AdissolMax = np.max(self.AdissolUnique)
                vivoList.append((self.tvivoUnique, self.FabsUnique))
                vitroList.append((self.tvitroUnique, self.AdissolUnique))

                self.BAdissol = InterpolatedUnivariateSpline(self.tvitroUnique, self.AdissolUnique, k=1)
                self.BFabs = InterpolatedUnivariateSpline(self.tvivoUnique, self.FabsUnique, k=1)

                self.bestError = 1e38
                optimum = differential_evolution(self.goalFunction,self.bounds,popsize=50)
                # self.verbose=True
                allCoeffs.append(optimum.x.tolist())


                # Evaluate correlation
                self.goalFunction(optimum.x)
                R = self.calculateR()
                allR.append(R)

                self.addSample("ivivc_%04d" % i, self.sampleNames[invivoIdx], self.vesselNames[invitroIdx], optimum.x, R, set=1)
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
        self.printFormulas(fh)
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

    def printFormulas(self, fh):
        self.doublePrint(fh,"Time scale: tvitro=%s"%self.timeScale.get().replace('$(t)','$(tvivo)'))
        self.doublePrint(fh,"Response scale: Fabs(t)=%s"%self.responseScale.get())

    def _summary(self):
        retval = []
        self.addFileContentToMessage(retval,self._getPath("summary.txt"))
        return retval
