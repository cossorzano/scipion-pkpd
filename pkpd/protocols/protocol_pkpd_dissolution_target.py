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

try:
    from itertools import izip
except ImportError:
    izip = zip

from math import sqrt
import numpy as np
import sys
from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.utils import uniqueFloatValues, computeXYmean, smoothPchip
from pkpd.pkpd_units import createUnit, PKPDUnit
from .protocol_pkpd import ProtPKPD

class ProtPKPDDissolutionTarget(ProtPKPD):
    """ Given an in-vivo absorption profile (assumed to be between 0 and 100), and the IVIVC parameters specified
        in this protocol, the protocol calculates the target in-vitro dissolution profile"""

    _label = 'dissol target'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputInVivo', params.PointerParam, label="Dissolution profiles in vivo",
                      pointerClass='ProtPKPDDeconvolve,ProtPKPDDeconvolutionWagnerNelson,ProtPKPDDeconvolutionLooRiegelman, PKPDExperiment',
                      help='Select an experiment with dissolution profiles')
        ts = form.addGroup("Time scaling (Delayed power scale (Fabs(t)=Adissol(k*(t-t0)^alpha))")
        ts.addParam('t0',params.FloatParam,label='t0',default=0,
                    help='Make sure it is in the same time units as the inputs')
        ts.addParam('k',params.FloatParam,label='k',default=1)
        ts.addParam('alpha',params.FloatParam,label='alpha',default=1)

        rs = form.addGroup("Response scaling (Affine transformation (Fabs(t)=A*Adissol(t)+B))")
        rs.addParam('A',params.FloatParam,label='A',default=1)
        rs.addParam('B',params.FloatParam,label='B',default=0)
        rs.addParam('saturate', params.BooleanParam, label='Saturate at 100%', default=True,
                    help='Saturate the calculated Adissol at 100%')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateTarget',self.inputInVivo.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def getInVivoProfiles(self):
        if hasattr(self.inputInVivo.get(),"outputExperiment"):
            fnPKPD = self.inputInVivo.get().outputExperiment.fnPKPD
        elif hasattr(self.inputInVivo.get(),"fnPKPD"):
            fnPKPD = self.inputInVivo.get().fnPKPD
        else:
            raise Exception("Cannot find a suitable filename for reading the experiment")
        experiment = self.readExperiment(fnPKPD)
        allY = []
        sampleNames = []
        for sampleName, sample in experiment.samples.items():
            t=sample.getValues("t")
            y=sample.getValues("A")
            allY.append((np.asarray(t,dtype=np.float64),np.asarray(y,dtype=np.float64)))
            sampleNames.append(sampleName)
        self.experimentInVivo = experiment
        return allY, sampleNames

    def addSample(self, sampleName):
        newSampleAdissol = PKPDSample()
        newSampleAdissol.sampleName = sampleName
        newSampleAdissol.variableDictPtr = self.outputExperiment.variables
        newSampleAdissol.descriptors = {}
        newSampleAdissol.addMeasurementColumn("tvitro", self.tvitroUnique)
        newSampleAdissol.addMeasurementColumn("Adissol",self.AdissolUnique)
        self.outputExperiment.samples[sampleName] = newSampleAdissol

    def produceAdissol(self,parameterInVitro,tmax):
        deltaT=np.min([(tmax+1)/1000,1.0])
        tvitro = np.arange(0,tmax+1,deltaT)
        if self.removeInVitroTlag:
            i=0
            for prmName in self.protFit.model.getParameterNames():
                if "tlag" in prmName:
                    parameterInVitro[i]=0.0
                i+=1
        self.protFit.model.x = tvitro
        Avitro = self.protFit.model.forwardModel(parameterInVitro)[0]
        return (tvitro, Avitro)

    def createOutputExperiment(self):
        tvitroVar = PKPDVariable()
        tvitroVar.varName = "tvitro"
        tvitroVar.varType = PKPDVariable.TYPE_NUMERIC
        tvitroVar.role = PKPDVariable.ROLE_TIME
        tvitroVar.units = createUnit(self.experimentInVivo.getTimeUnits().unit)
        tvitroVar.comment = "tvitro"

        AdissolVar = PKPDVariable()
        AdissolVar.varName = "Adissol"
        AdissolVar.varType = PKPDVariable.TYPE_NUMERIC
        AdissolVar.role = PKPDVariable.ROLE_MEASUREMENT
        AdissolVar.units = createUnit(PKPDUnit.UNIT_NONE)
        AdissolVar.comment = "Amount disolved in vitro"

        self.outputExperiment = PKPDExperiment()
        self.outputExperiment.variables[tvitroVar.varName] = tvitroVar
        self.outputExperiment.variables[AdissolVar.varName] = AdissolVar
        self.outputExperiment.general["title"] = "In-vitro target simulation"
        self.outputExperiment.general["comment"] = ""

    def calculateTarget(self, objId1):
        self.profilesInVivo, self.sampleNames=self.getInVivoProfiles()

        self.createOutputExperiment()

        # Compute all pairs
        for profileInVivo, sampleName in izip(self.profilesInVivo,self.sampleNames):
            tvivo=profileInVivo[0]
            Fabs=profileInVivo[1]
            FabsUnique, tvivoUnique = uniqueFloatValues(Fabs, tvivo)
            BFabs = InterpolatedUnivariateSpline(tvivoUnique, FabsUnique, k=1)

            tmax = np.max(tvivoUnique)
            deltaT = np.min([(tmax + 1) / 1000, 1.0])
            tvivop = np.arange(0, tmax + 1, deltaT)
            Fabsp = BFabs(tvivop)

            tvitro = self.k.get()*np.power(tvivop-self.t0.get(),self.alpha.get())
            Adissol = np.clip((Fabsp-self.B.get())/self.A.get(),0,None)
            if self.saturate:
                Adissol = np.clip(Adissol,None,100.0)

            self.AdissolUnique, self.tvitroUnique = uniqueFloatValues(Adissol, tvitro)

            self.addSample("target_%s" % sampleName)

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        self._defineSourceRelation(self.inputInVivo.get(), self.outputExperiment)
