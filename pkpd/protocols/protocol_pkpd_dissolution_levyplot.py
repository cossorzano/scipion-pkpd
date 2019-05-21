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
from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.utils import uniqueFloatValues
from pkpd.pkpd_units import createUnit
from .protocol_pkpd import ProtPKPD


# tested in test_workflow_levyplot

class ProtPKPDDissolutionLevyPlot(ProtPKPD):
    """ Calculate the Levy plot between two dissolution experiments. Each experiment may have
        several profiles and all vs all profiles are calculated"""

    _label = 'dissol Levy'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputInVitro', params.PointerParam, label="Dissolution profiles in vitro",
                      pointerClass='ProtPKPDDissolutionFit', help='Select an experiment with dissolution profiles')
        form.addParam('inputInVivo', params.PointerParam, label="Dissolution profiles in vivo",
                      pointerClass='ProtPKPDDeconvolve', help='Select an experiment with dissolution profiles')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllLevy',self.inputInVitro.get().getObjId(),self.inputInVivo.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def getInVivoProfiles(self):
        experiment = self.readExperiment(self.inputInVivo.get().outputExperiment.fnPKPD)
        allY = []
        for sampleName, sample in experiment.samples.iteritems():
            t=sample.getValues("t")
            y=sample.getValues("A")
            allY.append((np.asarray(t,dtype=np.float64),np.asarray(y,dtype=np.float64)))
        self.experimentInVivo = experiment
        return allY

    def getInVitroModels(self):
        allParameters = []
        self.protFit = self.inputInVitro.get()
        experiment = self.readExperiment(self.protFit.outputExperiment.fnPKPD)
        self.fitting = self.readFitting(self.protFit.outputFitting.fnFitting)
        self.varNameX = self.fitting.predictor.varName
        self.varNameY = self.fitting.predicted.varName

        self.protFit.model = self.protFit.createModel()
        self.protFit.model.setExperiment(experiment)
        self.protFit.model.setXVar(self.varNameX)
        self.protFit.model.setYVar(self.varNameY)
        self.protFit.setupFromFormParameters()
        self.protFit.experiment = experiment
        self.protFit.setupModel()

        parameterNames = self.protFit.model.getParameterNames()
        for sampleName, sample in experiment.samples.iteritems():
            parameters0 = []
            for parameterName in parameterNames:
                parameters0.append(float(sample.descriptors[parameterName]))
            allParameters.append(parameters0)
        self.experimentInVitro = experiment
        return allParameters

    def produceLevyPlot(self,tvitro,parameterInVitro,Avivo):
        Avivounique, tunique=uniqueFloatValues(Avivo,tvitro)
        B = InterpolatedUnivariateSpline(Avivounique, tunique, k=1)
        self.protFit.model.x=tvitro
        Avitro = self.protFit.model.forwardModel(parameterInVitro)[0]
        tvivo=[]
        for i in range(tvitro.shape[0]):
            tvivo.append(B(Avitro[i]))
            # print("tvitro=%f Avitro=%f Avivo=%f tvivoeqv=%f"%(tvitro[i],Avitro[i],Avivo[i],B(Avitro[i])))
        return (tvitro,np.asarray(tvivo,dtype=np.float64),Avitro)

    def addSample(self, sampleName, tvitro, tvivo):
        newSample = PKPDSample()
        newSample.sampleName = sampleName
        newSample.variableDictPtr = self.outputExperiment.variables
        newSample.descriptors = {}
        newSample.addMeasurementPattern(["tvivo"])
        newSample.addMeasurementColumn("tvitro", tvitro)
        newSample.addMeasurementColumn("tvivo",tvivo)
        self.outputExperiment.samples[sampleName] = newSample

    def calculateAllLevy(self, objId1, objId2):
        parametersInVitro=self.getInVitroModels()
        profilesInVivo=self.getInVivoProfiles()

        self.outputExperiment = PKPDExperiment()
        tvitrovar = PKPDVariable()
        tvitrovar.varName = "tvitro"
        tvitrovar.varType = PKPDVariable.TYPE_NUMERIC
        tvitrovar.role = PKPDVariable.ROLE_TIME
        tvitrovar.units = createUnit("min")

        tvivovar = PKPDVariable()
        tvivovar.varName = "tvivo"
        tvivovar.varType = PKPDVariable.TYPE_NUMERIC
        tvivovar.role = PKPDVariable.ROLE_MEASUREMENT
        tvivovar.units = createUnit("min")

        self.outputExperiment.variables[tvitrovar.varName] = tvitrovar
        self.outputExperiment.variables[tvivovar.varName] = tvivovar
        self.outputExperiment.general["title"] = "Levy plots"
        self.outputExperiment.general["comment"] = "Time in vivo vs time in vitro"

        i=1
        for parameterInVitro in parametersInVitro:
            for t,profileInVivo in profilesInVivo:
                tvitro, tvivo, _ = self.produceLevyPlot(t,parameterInVitro,profileInVivo)
                self.addSample("levy_%04d"%i,tvitro,tvivo)
                i+=1

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        self._defineSourceRelation(self.inputInVitro.get(), self.outputExperiment)
        self._defineSourceRelation(self.inputInVivo.get(), self.outputExperiment)

    def _validate(self):
        return []

    def _summary(self):
        retval = []
        return retval
