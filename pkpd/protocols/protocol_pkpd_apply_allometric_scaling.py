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
import math
import numpy as np

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDExperiment, PKPDAllometricScale
from pkpd.pkpd_units import strUnit
from pyworkflow.protocol.constants import LEVEL_ADVANCED


class ProtPKPDApplyAllometricScaling(ProtPKPD):
    """ Apply an allometric scaling previously calculated to an incoming experiment. The labels specified by the
        allometric scaling model will be rescaled to the target weight. Note that depending on the exponent of the
        fitting you may want to use a different predictor (weight*maximum lifespan potential, or weight*brain weight)
        see the rule of exponents (Mahmood and Balian 1996). """

    _label = 'apply allometric'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputPopulation', params.PointerParam, label="Input bootstrap population",
                      pointerClass='PKPDFitting',
                      help='The PK parameters of this experiment will be modified according to the allometric scale model.')
        form.addParam('inputAllometric', params.PointerParam, label="Allometric model", pointerClass='PKPDAllometricScale',
                      help='All variables specified by the allometric scale model will be adjusted')
        form.addParam('targetWeight', params.FloatParam, label="Target weight (kg)", default=25,
                      help='The PK parameters will be adjusted to this target weight')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runAdjust', self.targetWeight.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runAdjust(self, targetWeight):
        scaleModel = PKPDAllometricScale()
        scaleModel.load(self.inputAllometric.get().fnScale.get())

        self.population = self.readFitting(self.inputPopulation.get().fnFitting,cls="PKPDSampleFitBootstrap")
        self.experiment = PKPDExperiment()
        self.experiment.load(self.population.fnExperiment.get())

        for sampleFit in self.population.sampleFits:
            sample = self.experiment.samples[sampleFit.sampleName]
            sampleWeight = float(sample.getDescriptorValue(scaleModel.predictor))
            sample.setDescriptorValue(scaleModel.predictor,targetWeight)

            for varName, varUnits in scaleModel.averaged_vars:
                if varName in self.population.modelParameters:
                    idx = self.population.modelParameters.index(varName)
                    targetValue = scaleModel.models[varName][0]
                    sampleFit.parameters[0][idx] = targetValue

            for varName, varUnits in scaleModel.scaled_vars:
                if varName in self.population.modelParameters:
                    idx = self.population.modelParameters.index(varName)
                    k = scaleModel.models[varName][0]
                    a = scaleModel.models[varName][1]
                    targetValue = k*math.pow(targetWeight,a)
                    currentValue = k*math.pow(sampleWeight,a)
                    for j in range(sampleFit.parameters.shape[0]):
                        sampleFit.parameters[j][idx] *= targetValue/currentValue

        self.experiment.write(self._getPath("experiment.pkpd"))
        self.population.fnExperiment.set(self._getPath("experiment.pkpd"))
        self.population.write(self._getPath("bootstrapPopulation.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineOutputs(outputPopulation=self.population)
        self._defineSourceRelation(self.inputPopulation, self.experiment)
        self._defineSourceRelation(self.inputAllometric, self.experiment)
        self._defineSourceRelation(self.inputPopulation, self.population)
        self._defineSourceRelation(self.inputAllometric, self.population)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = ["Target weight: %f"%self.targetWeight.get()]
        return msg

    def _citations(self):
        return ['Sharma2009','Mahmood1996']