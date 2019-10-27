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
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.protocol.params as params
from .protocol_pkpd_dissolution_ivivc_splines import ProtPKPDDissolutionIVIVCSplines, ProtPKPDDissolutionIVIVCGeneric
from pkpd.utils import parseOperation


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
        # self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def calculateAllIvIvC(self):
        self.parametersInVitro = []
        self.vesselNames = []
        self.profilesInVivo = []
        self.sampleNames = []
        idx=1
        for ptrProt in self.inputIVIVCs:
            parametersInVitro, vesselNames = ptrProt.get().getInVitroModels()
            profilesInVivo, sampleNames = ptrProt.get().getInVivoProfiles()
            self.parametersInVitro.append(parametersInVitro)
            self.vesselNames.append(vesselNames)
            self.profilesInVivo.append(profilesInVivo)
            self.sampleNames.append(sampleNames)

            if idx==1:
                self.experimentInVitro = ptrProt.get().experimentInVitro
                self.experimentInVivo = ptrProt.get().experimentInVivo
                self.varNameX = ptrProt.get().fitting.predictor.varName
                self.varNameY = ptrProt.get().fitting.predicted.varName
            idx+=1

        self.createOutputExperiments(set=1)

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

#    def createOutputStep(self):
#        self._defineOutputs(outputExperimentFabsSingle=self.outputExperimentFabsSingle)
#        for ptrExperiment in self.inputIVIVCs:
#            self._defineSourceRelation(ptrExperiment.get(), self.outputExperimentFabsSingle)
