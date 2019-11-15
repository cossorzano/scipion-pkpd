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

import copy
from itertools import izip
import numpy as np
import sys
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import differential_evolution

import pyworkflow.protocol.params as params
from .protocol_pkpd_dissolution_ivivc_splines import ProtPKPDDissolutionIVIVCSplines, ProtPKPDDissolutionIVIVCGeneric
from pkpd.pkpd_units import PKPDUnit
from pkpd.objects import PKPDSample
from pkpd.utils import parseOperation, uniqueFloatValues


class ProtPKPDDissolutionIVIVCJoinRecalculate(ProtPKPDDissolutionIVIVCSplines):
    """ Join several IVIVCs into a single one. The strategy is to compute the average of all the plots involved in the
        IVIVC process: 1) tvivo -> tvitro; 2) tvitro -> Adissol; 3) Adissol->FabsPredicted. The plot tvivo-Fabs comes
        after the IVIVC process, while the plot tvivo-FabsOrig is the observed one in the input files. These two
        plots need not be exactly the same. """

    _label = 'dissol ivivc join recalculate'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputIVIVCs', params.MultiPointerParam, label="IVIVCs",
                      pointerClass='ProtPKPDDissolutionIVIVC, ProtPKPDDissolutionIVIVCGeneric, ProtPKPDDissolutionIVIVCSplines', help='Choose runs of IVIV correlations')
        form.addParam('removeInVitroTlag',params.BooleanParam, label="Remove in vitro tlag", default=True,
                      help="If there is an in vitro tlag, set it to 0")
        form.addParam('timeScale', params.EnumParam, label="Time scaling", default=0,
                      choices=["None (Fabs(t)=Adissol(t))",
                               "t0 (Fabs(t)=Adissol(t-t0)",
                               "Linear scale (Fabs(t)=Adissol(k*t))",
                               "Affine transformation (Fabs(t)=Adissol(k*(t-t0))",
                               "Power scale (Fabs(t)=Adissol(k*t^alpha)",
                               "Delayed power scale (Fabs(t)=Adissol(k*(t-t0)^alpha)",
                               "Generic","SplineXY0","SplineXY1","SplineXY2","SplineXY3","SplineXY4","SplineXY5"],
                      help='Fabs is the fraction absorbed in vivo, while Adissol is the amount dissolved in vitro.')
        form.addParam('timeScaleOperation', params.StringParam, label="Time scaling", default="$(t)", condition="timeScale==6",
                      help='Write a time scaling function like $(t), or $[a0]+$[a1]*$(t)+$[a2]*np.power($(t),2). '
                           'The meaning is tvitro=a0+a1*tvivo+a2*tvivo^2 in the last example. '
                           'You have all numpy functions at your disposal: np.sqrt, np.log10, np.log, np.exp')
        form.addParam('timeBounds',params.StringParam,label='Time coefficients bounds',default='', condition="timeScale>=1 and timeScale<=6",
                      help='a0: [-100,100]; a1: [0.0001,100]; a2: [1e-6,1e-4]')
        form.addParam('responseScale', params.EnumParam, label="Response scaling",
                      choices=["None",
                               "Linear scale (Fabs(t)=A*Adissol(t))",
                               "Affine transformation (Fabs(t)=A*Adissol(t)+B",
                               "Generic","SplineXY0","SplineXY1","SplineXY2","SplineXY3","SplineXY4","SplineXY5"], default=0,
                      help='Fabs is the fraction absorbed in vivo, while Adissol is the amount dissolved in vitro.')
        form.addParam('responseScaleOperation', params.StringParam, label="Response scaling", default="$(Adissol)", condition="responseScale==3",
                      help='Write a response scaling function like $[k1]*$(Adissol)+$[k2]*np.power($(Adissol),2)$. '
                           'Make sure that the time and response scaling coefficients have different names. '
                           'The meaning is Fabs=k1*Adissol+k2*Adissol^2 in the last example. '
                           'You have all numpy functions at your disposal: np.sqrt, np.log10, np.log, np.exp')
        form.addParam('responseBounds',params.StringParam,label='Response coefficients bounds',default='', condition="responseScale>=1 and responseScale<=3",
                      help='k1: [0.0001,100]; k2: [1e-6,1e-4]')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllIvIvC')
        self._insertFunctionStep('createOutputStep')

    def predict(self, tvivoUnique, tvitroUnique, tvivoXUnique, tvitroXUnique, adissolXUnique0, fabsXUnique0, FabsUnique, AdissolUnique, tvivoMin, tvivoMax):
        self.tvivoUnique = tvivoUnique
        self.tvitroUnique = tvitroUnique

        # Forward error
        if self.parsedTimeOperation is None:
            # Splines for time
            Bt = InterpolatedUnivariateSpline(tvivoXUnique, tvitroXUnique, k=1)
            self.tvitroUnique1 = Bt(tvivoUnique)
        else:
            t = tvivoUnique
            self.tvitroUnique1 = eval(self.parsedTimeOperation)

        self.tvitroReinterpolated = np.clip(self.tvitroUnique1, self.tvitroMin, self.tvitroMax)
        self.AdissolReinterpolated = self.BAdissol(self.tvitroReinterpolated)

        if self.parsedResponseOperation is None:
            adissolXUnique = self.AdissolMaxx * adissolXUnique0
            fabsXUnique = self.FabsMaxx * fabsXUnique0
            BA = InterpolatedUnivariateSpline(adissolXUnique, fabsXUnique, k=1)
            FabsPredicted = BA(self.AdissolReinterpolated)
        else:
            Adissol = self.AdissolReinterpolated
            FabsPredicted = eval(self.parsedResponseOperation)
        self.FabsUnique = FabsUnique
        self.AdissolUnique = AdissolUnique
        self.FabsPredicted = np.clip(FabsPredicted, 0.0, 100)

        # Backward error
        tvitroAux, tvivoAux = uniqueFloatValues(self.tvitroUnique1, tvivoUnique)
        Btinv = InterpolatedUnivariateSpline(tvitroAux, tvivoAux, k=1)
        self.tvivoReinterpolated = np.clip(Btinv(tvitroUnique), tvivoMin, tvivoMax)
        self.FabsReinterpolated = self.BFabs(self.tvivoReinterpolated)
        FabsPredictedAux, AdissolAux = uniqueFloatValues(self.FabsPredicted, self.AdissolReinterpolated)
        Bfinv = InterpolatedUnivariateSpline(FabsPredictedAux, AdissolAux, k=1)
        self.AdissolPredicted = np.clip(Bfinv(self.FabsReinterpolated), 0.0, 100)

    def getAxes(self,x):
        # Get time parameters
        if self.parsedTimeOperation is None:
            # Splines for time
            tvivoXUnique, tvitroXUnique, i0 = ProtPKPDDissolutionIVIVCSplines.getParameters(self, x, 0,
                                                                                            self.coeffTimeList,
                                                                                            'tvitro')
            tvivoXUnique = self.tvivoMaxx * tvivoXUnique
            tvitroXUnique = self.tvitroMaxx * tvitroXUnique
        else:
            i = 0
            for prm in self.coeffTimeList:
                exec ("%s=%f" % (prm, x[i]))
                i += 1
            i0 = i
            tvivoXUnique = None
            tvitroXUnique = None

        # Get response parameters
        if self.parsedResponseOperation is None:
            # Splines for response
            fabsXUnique0, adissolXUnique0, _ = ProtPKPDDissolutionIVIVCSplines.getParameters(self, x, i0,
                                                                                             self.coeffResponseList,
                                                                                             'adissol')
        else:
            i = i0
            for prm in self.coeffResponseList:
                exec ("%s=%f" % (prm, x[i]))
                i += 1
            adissolXUnique0 = None
            fabsXUnique0 = None
        return tvivoXUnique, tvitroXUnique, adissolXUnique0, fabsXUnique0

    def goalFunction(self,x):
        try:
            tvivoXUnique, tvitroXUnique, adissolXUnique0, fabsXUnique0=self.getAxes(x)

            error=0
            errorBackward=0
            errorForward=0
            N=0
            self.allR=[]
            i=0
            for tvivoUnique, tvitroUnique, FabsUnique, AdissolUnique, BAdissol, BFabs, tvitroMin, tvitroMax, tvivoMin, tvivoMax, AdissolMax, FabsMax, _, _ in self.allPairs:
                i+=1
                self.predict(tvivoUnique, tvitroUnique, tvivoXUnique, tvitroXUnique, adissolXUnique0, fabsXUnique0, FabsUnique, AdissolUnique, tvivoMin, tvivoMax)

                errorn, errorBackwardn, errorForwardn = self.calculateIndividualError(x, self.tvitroReinterpolated, self.tvivoReinterpolated)
                error+=errorn
                errorBackward+=errorBackwardn
                errorForward+=errorForwardn
                R=ProtPKPDDissolutionIVIVCSplines.calculateR(self)
                self.allR.append(R)
                N += 1
        except:
            return 1e38

        if N>0:
            error/=N
            errorForward/=N
            errorBackward/=N
            if error < self.bestError:
                print("New minimum error=%f (back=%f, forw=%f) R=%f" % (
                error, errorBackward, errorForward, self.calculateR()),
                      "x=%s" % np.array2string(x, max_line_width=1000))
                self.bestError = error
                sys.stdout.flush()
            return error
        else:
            return 1e38

    def calculateR(self):
        return np.mean(self.allR)

    #--------------------------- STEPS functions --------------------------------------------
    def getTimeMsg(self):
        if self.timeScale.get()==0:
            return "tvitro=tvivo"
        elif self.timeScale.get()==1:
            return "tvitro=tvivo-t0"
        elif self.timeScale.get() == 2:
            return "tvitro=k*tvivo"
        elif self.timeScale.get()==3:
            return "tvitro=k*(tvivo-t0)"
        elif self.timeScale.get()==4:
            return "tvitro=k*tvivo^alpha"
        elif self.timeScale.get()==5:
            return "tvitro=k*(tvivo-t0)^alpha"
        elif self.timeScale.get()==6:
            return "tvitro=%s"%self.timeScaleOperation.get().replace('$(t)', '$(tvivo)')
        else:
            return "tvitro=SplineXY%d(tvivo)"%(self.timeScale.get()-7)

    def getResponseMsg(self):
        if self.responseScale.get()==0:
            return "FabsPredicted=Adissol",
        elif self.responseScale.get()==1:
            return "Fabs=A*Adissol",
        elif self.responseScale.get()==2:
            return "Fabs=A*Adissol+B"
        elif self.responseScale.get()==3:
            return "Fabs=%s"%self.responseScaleOperation.get()
        else:
            return "Fabs=SplineXY%d(Adissol)" % (self.responseScale.get() - 4)

    def addParametersToExperiment(self, outputExperiment, sampleName, individualFrom, vesselFrom, optimum, R):
        i=0
        timeScaleMsg=self.getTimeMsg()
        for prm in self.coeffTimeList:
            outputExperiment.addParameterToSample(sampleName, prm, PKPDUnit.UNIT_NONE, timeScaleMsg, optimum[i])
            i+=1
        responseMsg=self.getResponseMsg()
        for prm in self.coeffResponseList:
            outputExperiment.addParameterToSample(sampleName, prm, PKPDUnit.UNIT_NONE, responseMsg, optimum[i])
            i+=1

        outputExperiment.addParameterToSample(sampleName, "R", PKPDUnit.UNIT_NONE, "IVIV Correlation coefficient", R)
        outputExperiment.addLabelToSample(sampleName, "from", "individual---vesel", "%s---%s"%(individualFrom,vesselFrom))

    def addSample(self, outputExperiment, sampleName, individualFrom, vesselFrom, optimum, R, addFabs=False):
        newSampleFabs = PKPDSample()
        newSampleFabs.sampleName = sampleName
        newSampleFabs.variableDictPtr = self.outputExperimentFabsSingle.variables
        newSampleFabs.descriptors = {}
        newSampleFabs.addMeasurementColumn("tvitroReinterpolated", self.tvitroReinterpolated)
        newSampleFabs.addMeasurementColumn("AdissolReinterpolated", self.AdissolReinterpolated)
        newSampleFabs.addMeasurementColumn("tvivo", self.tvivoUnique)
        newSampleFabs.addMeasurementColumn("FabsPredicted", self.FabsPredicted)
        if addFabs:
            newSampleFabs.addMeasurementColumn("Fabs", self.FabsUnique)
        outputExperiment.samples[sampleName] = newSampleFabs
        self.addParametersToExperiment(outputExperiment, sampleName, individualFrom, vesselFrom, optimum, R)

    def calculateAllIvIvC(self):
        # Get the PK and dissolution profiles from the input
        self.parametersInVitro = []
        self.vesselNames = []
        self.profilesInVivo = []
        self.sampleNames = []
        self.experimentsInVitro = []
        idx=1
        for ptrProt in self.inputIVIVCs:
            parametersInVitro, vesselNames = ptrProt.get().getInVitroModels()
            profilesInVivo, sampleNames = ptrProt.get().getInVivoProfiles()
            self.parametersInVitro.append(parametersInVitro)
            self.vesselNames.append(vesselNames)
            self.profilesInVivo.append(profilesInVivo)
            self.sampleNames.append(sampleNames)

            self.experimentInVitro = ptrProt.get().experimentInVitro
            self.experimentInVivo = ptrProt.get().experimentInVivo
            if idx==1:
                self.varNameX = ptrProt.get().fitting.predictor.varName
                self.varNameY = ptrProt.get().fitting.predicted.varName
                self.protFit = ptrProt.get().protFit
            self.experimentsInVitro.append(self.experimentInVitro)
            idx+=1

        # Prepare all data and pairs for processing
        self.allPairs = []
        self.tvitroMaxx = -1e38
        self.tvivoMaxx = -1e38
        self.FabsMaxx = -1e38
        self.AdissolMaxx = -1e38
        for block in range(len(self.sampleNames)):
            for parameterInVitro, vesselName in izip(self.parametersInVitro[block], self.vesselNames[block]):
                if "tvitroMax" in self.experimentsInVitro[block].variables:
                    tvitroMax = float(self.experimentsInVitro[block].samples[vesselName].getDescriptorValue("tvitroMax"))
                else:
                    tvitroMax = 1e38
                j=0
                for self.tvivo,self.Fabs in self.profilesInVivo[block]:
                    self.FabsUnique, self.tvivoUnique = uniqueFloatValues(self.Fabs, self.tvivo)
                    self.tvitro, self.Adissol = self.produceAdissol(parameterInVitro, min(np.max(self.tvivoUnique * 10), tvitroMax))
                    self.AdissolUnique, self.tvitroUnique = uniqueFloatValues(self.Adissol, self.tvitro)
                    self.tvivoMin = np.min(self.tvivoUnique)
                    self.tvivoMax = np.max(self.tvivoUnique)
                    self.tvitroMin = np.min(self.tvitroUnique)
                    self.tvitroMax = np.max(self.tvitroUnique)
                    self.tvitroMaxx = max(self.tvitroMaxx,self.tvitroMax)
                    self.tvivoMaxx = max(self.tvivoMaxx,self.tvivoMax)
                    self.FabsMaxx = max(self.FabsMaxx,np.max(self.FabsUnique))
                    self.AdissolMaxx = max(self.AdissolMaxx,np.max(self.AdissolUnique))

                    # Make sure they are sorted in x
                    self.tvivoUnique, self.FabsUnique = uniqueFloatValues(self.tvivoUnique, self.FabsUnique)
                    self.tvitroUnique, self.AdissolUnique = uniqueFloatValues(self.tvitroUnique, self.AdissolUnique)
                    self.FabsMax = np.max(self.FabsUnique)
                    self.AdissolMax = np.max(self.AdissolUnique)

                    self.BAdissol = InterpolatedUnivariateSpline(self.tvitroUnique, self.AdissolUnique, k=1)
                    self.BFabs = InterpolatedUnivariateSpline(self.tvivoUnique, self.FabsUnique, k=1)
                    self.allPairs.append([self.tvivoUnique, self.tvitroUnique, self.FabsUnique, self.AdissolUnique,
                                          self.BAdissol, self.BFabs,
                                          self.tvitroMin, self.tvitroMax, self.tvivoMin, self.tvivoMax,
                                          self.AdissolMax, self.FabsMax, vesselName, self.sampleNames[j]])
                    j+=1

        # Prepare the parameter names and bounds
        if self.timeScale.get()<=6:
            if self.timeScale.get()==0:
                timeScaleOperation="$(t)"
            elif self.timeScale.get()==1:
                timeScaleOperation = "$(t)-$[t0]"
            elif self.timeScale.get() == 2:
                timeScaleOperation = "$[k]*$(t)"
            elif self.timeScale.get() == 3:
                timeScaleOperation = "$[k]*($(t)-$[t0])"
            elif self.timeScale.get() == 4:
                timeScaleOperation = "$[k]*np.power($(t),$[alpha])"
            elif self.timeScale.get() == 5:
                timeScaleOperation = "$[k]*np.power($(t)-$[t0],$[alpha])"
            else:
                timeScaleOperation = self.timeScaleOperation.get()
            self.parsedTimeOperation, self.varTimeList, self.coeffTimeList = parseOperation(timeScaleOperation)
            self.timeBoundsList = ProtPKPDDissolutionIVIVCGeneric.constructBounds(self, self.coeffTimeList, self.timeBounds.get())
        else:
            self.parsedTimeOperation = None # It is a spline
            self.timeSplinesN = self.timeScale.get()-6
            self.coeffTimeList = self.constructTimeCoeffs(self.timeSplinesN)
            self.timeBoundsList = ProtPKPDDissolutionIVIVCSplines.constructBounds(self, self.coeffTimeList, "")

        if self.responseScale.get()<=3:
            if self.responseScale.get()==0:
                responseScaleOperation="$(Adissol)"
            elif self.responseScale.get()==1:
                responseScaleOperation="$[A]*$(Adissol)"
            elif self.responseScale.get()==2:
                responseScaleOperation="$[A]*$(Adissol)+$[B]"
            elif self.responseScale.get()==3:
                responseScaleOperation=self.responseScaleOperation.get()
            self.parsedResponseOperation, self.varResponseList, self.coeffResponseList = parseOperation(responseScaleOperation)
            self.responseBoundsList = ProtPKPDDissolutionIVIVCGeneric.constructBounds(self, self.coeffResponseList, self.responseBounds.get())
        else:
            self.parsedResponseOperation = None # It is a spline
            self.responseSplinesN = self.responseScale.get()-3
            self.coeffResponseList = self.constructResponseCoeffs(self.responseSplinesN)
            self.responseBoundsList = ProtPKPDDissolutionIVIVCSplines.constructBounds(self, self.coeffResponseList, "")

        self.parameters = self.coeffTimeList + self.coeffResponseList
        self.bounds = self.timeBoundsList + self.responseBoundsList

        self.createOutputExperiments(set=2)

        self.verbose=False
        self.bestError = 1e38
        optimum = differential_evolution(self.goalFunction, self.bounds, popsize=50)

        # Evaluate correlation
        self.goalFunction(optimum.x)
        R = self.calculateR()

        self.addSample(self.outputExperimentFabsSingle, "ivivc_all", "allSamples", "allVessels", optimum.x, R)
        self.outputExperimentFabsSingle.write(self._getPath("experimentFabsSingle.pkpd"))

        fh = open(self._getPath("summary.txt"), "w")
        self.doublePrint(fh, "Correlation coefficient (R) %f"%R)
        self.printFormulas(fh)
        fh.close()

        # Generate Fabs and FabsPredicted for all pairs
        self.createOutputExperiments(set=1)
        tvivoXUnique, tvitroXUnique, adissolXUnique0, fabsXUnique0 = self.getAxes(optimum.x)
        i=1
        for tvivoUnique, tvitroUnique, FabsUnique, AdissolUnique, BAdissol, BFabs, tvitroMin, tvitroMax, tvivoMin, tvivoMax, AdissolMax, FabsMax, vesselName, sampleName in self.allPairs:
            self.predict(tvivoUnique, tvitroUnique, tvivoXUnique, tvitroXUnique, adissolXUnique0, fabsXUnique0, FabsUnique, AdissolUnique, tvivoMin, tvivoMax)
            R = ProtPKPDDissolutionIVIVCSplines.calculateR(self)
            self.addSample(self.outputExperimentFabs, "ivivc_%d"%i, sampleName, vesselName, optimum.x, R, True)
            i+=1
        self.outputExperimentFabs.write(self._getPath("experimentFabs.pkpd"))

    def printFormulas(self, fh):
        self.doublePrint(fh, "Time scale: %s" % self.getTimeMsg())
        self.doublePrint(fh, "Response scale: %s" % self.getResponseMsg())

    def createOutputStep(self):
        self._defineOutputs(outputExperimentFabsSingle=self.outputExperimentFabsSingle)
        self._defineOutputs(outputExperimentFabs=self.outputExperimentFabs)
        for ptrProt in self.inputIVIVCs:
            self._defineSourceRelation(ptrProt.get().inputInVitro.get(), self.outputExperimentFabsSingle)
            self._defineSourceRelation(ptrProt.get().inputInVivo.get(), self.outputExperimentFabsSingle)
            self._defineSourceRelation(ptrProt.get().inputInVitro.get(), self.outputExperimentFabs)
            self._defineSourceRelation(ptrProt.get().inputInVivo.get(), self.outputExperimentFabs)

