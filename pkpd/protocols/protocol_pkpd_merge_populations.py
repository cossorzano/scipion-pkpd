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

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDFitting, PKPDSampleFitBootstrap
import numpy as np
import copy


class ProtPKPDMergePopulations(ProtPKPD):
    """ Merge two populations. Both populations must have the same labels\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'merge populations'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputPopulation1', params.PointerParam, label="Population 1", important=True,
                      pointerClass='PKPDFitting', pointerCondition="isPopulation",
                      help='It must be a fitting coming from a bootstrap sample')
        form.addParam('inputPopulation2', params.PointerParam, label="Population 2", important=True,
                      pointerClass='PKPDFitting', pointerCondition="isPopulation",
                      help='It must be a fitting coming from a bootstrap sample')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runMerge',self.inputPopulation1.get().getObjId(), self.inputPopulation2.getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runMerge(self, objId1, objId2):
        self.population1 = self.readFitting(self.inputPopulation1.get().fnFitting,cls="PKPDSampleFitBootstrap")
        self.population2 = self.readFitting(self.inputPopulation2.get().fnFitting,cls="PKPDSampleFitBootstrap")
        self.printSection("Merging populations")

        self.fitting = PKPDFitting("PKPDSampleFitBootstrap")
        self.fitting.fnExperiment.set(self.population1.fnExperiment)
        self.fitting.predictor=self.population1.predictor
        self.fitting.predicted=self.population1.predicted
        self.fitting.modelParameterUnits = self.population1.modelParameterUnits
        self.fitting.modelParameters = self.population1.modelParameters
        self.fitting.modelDescription = self.population1.modelDescription

        newSampleFit = PKPDSampleFitBootstrap()
        newSampleFit.sampleName = "Merged population"
        newSampleFit.parameters = None
        for sampleFit in self.population1.sampleFits:
            if newSampleFit.parameters == None:
                newSampleFit.parameters = np.copy(sampleFit.parameters)
                newSampleFit.xB = copy.copy(sampleFit.xB)
                newSampleFit.yB = copy.copy(sampleFit.yB)
                newSampleFit.R2 = copy.copy(sampleFit.R2)
                newSampleFit.R2adj = copy.copy(sampleFit.R2adj)
                newSampleFit.AIC = copy.copy(sampleFit.AIC)
                newSampleFit.AICc = copy.copy(sampleFit.AICc)
                newSampleFit.BIC = copy.copy(sampleFit.BIC)
            else:
                newSampleFit.parameters = np.vstack([newSampleFit.parameters, sampleFit.parameters])
                newSampleFit.xB += sampleFit.xB
                newSampleFit.yB += sampleFit.yB
                newSampleFit.R2 += sampleFit.R2
                newSampleFit.R2adj += sampleFit.R2adj
                newSampleFit.AIC += sampleFit.AIC
                newSampleFit.AICc += sampleFit.AICc
                newSampleFit.BIC += sampleFit.BIC
        for sampleFit in self.population2.sampleFits:
            newSampleFit.parameters = np.vstack([newSampleFit.parameters, sampleFit.parameters])
            print(type(newSampleFit.xB))
            print(type(sampleFit.xB))
            newSampleFit.xB += sampleFit.xB
            newSampleFit.yB += sampleFit.yB
            newSampleFit.R2 += sampleFit.R2
            newSampleFit.R2adj += sampleFit.R2adj
            newSampleFit.AIC += sampleFit.AIC
            newSampleFit.AICc += sampleFit.AICc
            newSampleFit.BIC += sampleFit.BIC
        self.fitting.sampleFits.append(newSampleFit)

        self.fitting.write(self._getPath("bootstrapPopulation.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputPopulation=self.fitting)
        self._defineSourceRelation(self.inputPopulation1, self.fitting)
        self._defineSourceRelation(self.inputPopulation2, self.fitting)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=["Populations %s and %s were merged"%(self.getObjectTag(self.inputPopulation1.get()),
                                                  self.getObjectTag(self.inputPopulation2.get()))]
        return msg

    def _validate(self):
        msg=[]
        if not self.inputPopulation1.get().fnFitting.get().endswith("bootstrapPopulation.pkpd"):
            msg.append("Population 1 must be a bootstrap sample")
        if not self.inputPopulation2.get().fnFitting.get().endswith("bootstrapPopulation.pkpd"):
            msg.append("Population 2 must be a bootstrap sample")
        return msg