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
from pkpd.utils import uniqueFloatValues
from pkpd.pkpd_units import PKPDUnit
from .protocol_pkpd_dissolution_ivivc_generic import ProtPKPDDissolutionIVIVCGeneric

# tested in test_workflow_levyplot2

class ProtPKPDDissolutionIVIVCSplines(ProtPKPDDissolutionIVIVCGeneric):
    """ Calculate the in vitro-in vivo correlation between two experiments. Each experiment may have
        several profiles and all vs all profiles are calculated. You may scale the time between the two
        sets of experiments"""

    _label = 'dissol ivivc splines'

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
        form.addParam('timeScale', params.EnumParam, label="Time scaling", default=2,
                      choices=["SplineXY0","SplineXY1","SplineXY2","SplineXY3","SplineXY4","SplineXY5"])
        form.addParam('responseScale', params.EnumParam, label="Response scaling", default=2,
                      choices=["SplineXY0","SplineXY1","SplineXY2","SplineXY3","SplineXY4","SplineXY5"])

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllIvIvC',self.inputInVitro.get().getObjId(),self.inputInVivo.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addParametersToExperiment(self, outputExperiment, sampleName, individualFrom, vesselFrom, optimum, R):
        i=0
        timeScaleMsg="SplineXY%d"%self.timeSplinesN
        for prm in self.coeffTimeList:
            outputExperiment.addParameterToSample(sampleName, prm, PKPDUnit.UNIT_TIME_MIN, timeScaleMsg, optimum[i])
            i+=1
        responseMsg="SplineXY%d"%self.responseSplinesN
        for prm in self.coeffResponseList:
            outputExperiment.addParameterToSample(sampleName, prm, PKPDUnit.UNIT_NONE, responseMsg, optimum[i])
            i+=1

        outputExperiment.addParameterToSample(sampleName, "R", PKPDUnit.UNIT_NONE, "IVIV Correlation coefficient", R)
        outputExperiment.addLabelToSample(sampleName, "from", "individual---vesel", "%s---%s"%(individualFrom,vesselFrom))

    def getParameters(self, x, i0, parameterList, vitroPrefix):
        # print("Unsorted x",x)
        # Get parameters from x
        vitroX = []
        vivoX = []
        i = i0
        for prm in parameterList:
            if prm.startswith(vitroPrefix):
                vitroX.append(x[i])
            else:
                vivoX.append(x[i])
            i += 1
        vitroX = np.sort(vitroX)
        vivoX = np.sort(vivoX)
        i = i0
        it = 0
        iv = 0 # Copy the sorted vector back to x
        for prm in parameterList:
            if prm.startswith(vitroPrefix):
                x[i] = vitroX[it]
                it += 1
            else:
                x[i] = vivoX[iv]
                iv += 1
            i += 1
        vitroX = np.concatenate([[0.0],vitroX,[1.0]])
        vivoX = np.concatenate([[0.0],vivoX,[1.0]])
        tvivoXUnique, tvitroXUnique = uniqueFloatValues(vivoX, vitroX)
        # print("Sorted x",x)
        return tvivoXUnique, tvitroXUnique, i

    def goalFunction(self,x):
        try:
            tvivoXUnique, tvitroXUnique, i0 = self.getParameters(x,0,self.coeffTimeList,'tvitro')
            tvivoXUnique = self.tvivoMax * tvivoXUnique
            tvitroXUnique = self.tvitroMax * tvitroXUnique

            fabsXUnique, adissolXUnique, _ = self.getParameters(x,i0,self.coeffResponseList,'adissol')
            adissolXUnique = self.AdissolMax * adissolXUnique
            fabsXUnique = self.FabsMax * fabsXUnique

            Bt = InterpolatedUnivariateSpline(tvivoXUnique,tvitroXUnique,k=1)
            tvitroUnique=Bt(self.tvivoUnique)
            self.tvitroReinterpolated=np.clip(tvitroUnique,self.tvitroMin,self.tvitroMax)
            self.AdissolReinterpolated = self.BAdissol(self.tvitroReinterpolated)
            BA = InterpolatedUnivariateSpline(adissolXUnique,fabsXUnique,k=1)
            FabsPredicted = BA(self.AdissolReinterpolated)
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
        return [(0,1)]*len(coeffList)

    def constructTimeCoeffs(self, N):
        coeffTimeList=[]
        for i in range(N):
            coeffTimeList+=['tvivo%d'%i,'tvitro%d'%i]
        return coeffTimeList

    def constructResponseCoeffs(self, N):
        coeffResponseList=[]
        for i in range(N):
            coeffResponseList+=['adissol%d'%i,'fabs%d'%i]
        return coeffResponseList

    def calculateAllIvIvC(self, objId1, objId2):
        self.parametersInVitro, self.vesselNames=self.getInVitroModels()
        self.profilesInVivo, self.sampleNames=self.getInVivoProfiles()

        self.createOutputExperiments(set=1)
        self.timeSplinesN = self.timeScale.get()
        self.responseSplinesN = self.responseScale.get()

        self.coeffTimeList=self.constructTimeCoeffs(self.timeSplinesN)
        self.coeffResponseList=self.constructResponseCoeffs(self.responseSplinesN)
        self.parameters=self.coeffTimeList+self.coeffResponseList
        self.bounds = self.constructBounds(self.parameters,"")

        self.calculateAllIvIvCLoop()

    def printFormulas(self, fh):
        self.doublePrint(fh,"Time scale: tvitro=SplineXY%d(tvivo)"%self.timeSplinesN)
        self.doublePrint(fh,"Response scale: Fabs(t)=SplineXY%d(Adissol)"%self.responseSplinesN)
