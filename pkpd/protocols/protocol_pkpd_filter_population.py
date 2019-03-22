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

# TESTED in test_workflow_gabrielsson_pk02.py

class ProtPKPDFilterPopulation(ProtPKPD):
    """ Filter a population by some criterion\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'filter population'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputPopulation', params.PointerParam, label="Population", important=True,
                      pointerClass='PKPDFitting',
                      help='It must be a fitting coming from a bootstrap sample')
        form.addParam('filterType', params.EnumParam, choices=["Exclude","Keep","Keep confidence interval"], label="Filter mode", default=0,
                      help='Exclude or keep samples meeting the following condition\n')
        form.addParam('condition', params.StringParam, label="Condition", default="$(AICc)>-10 or $(R2)<0.7",
                      help='Example: $(Cl)>0.25 and $(tlag)>0\n'\
                           'Example for the confidence interval: $(Cl) 95 (keep the central 95% of clearance).\n'\
                           'The variables R2, R2adj, AIC, AICc, and BIC can be used, \n'\
                           'e.g., $(AICc)>-10 selects those elements with suspicious fitting.\n'\
                           'Confidence intervals cannot be placed on the quality parameters (R2, R2adj, ...)')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runFilter',self.inputPopulation.get().getObjId(), self.filterType.get(), self.condition.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runFilter(self, objId, filterType, condition):
        self.population = self.readFitting(self.inputPopulation.get().fnFitting,cls="PKPDSampleFitBootstrap")
        self.printSection("Filtering population")

        self.fitting = PKPDFitting("PKPDSampleFitBootstrap")
        self.fitting.fnExperiment.set(self.population.fnExperiment)
        self.fitting.predictor=self.population.predictor
        self.fitting.predicted=self.population.predicted
        self.fitting.modelParameterUnits = self.population.modelParameterUnits
        self.fitting.modelParameters = self.population.modelParameters
        self.fitting.modelDescription = self.population.modelDescription

        newSampleFit = PKPDSampleFitBootstrap()
        newSampleFit.parameters = None
        filterType = self.filterType.get()
        for sampleFit in self.population.sampleFits:
            newSampleFit.sampleName = sampleFit.sampleName
            if newSampleFit.parameters == None:
                newSampleFit.parameters = np.empty((0,sampleFit.parameters.shape[1]))
                newSampleFit.xB = []
                newSampleFit.yB = []
                newSampleFit.R2 = []
                newSampleFit.R2adj = []
                newSampleFit.AIC = []
                newSampleFit.AICc = []
                newSampleFit.BIC = []

            if filterType<=1:
                conditionToEvaluate = self.condition.get()
            else:
                tokens=self.condition.get().strip().split(' ')
                variable = tokens[0][2:-1]
                confidenceLevel = float(tokens[1])
                columnIdx = -1
                for j in range(len(self.population.modelParameters)):
                    if self.population.modelParameters[j]==variable:
                        columnIdx = j
                        break
                if columnIdx==-1:
                    raise Exception("Cannot find %s amongst the model variables"%variable)
                values=sampleFit.parameters[:,columnIdx]
                alpha_2 = (100-confidenceLevel)/2
                limits = np.percentile(values,[alpha_2,100-alpha_2])
                conditionToEvaluate = "%s>=%f and %s<=%f"%(tokens[0],limits[0],tokens[0],limits[1])
                print("Condition to evaluate: %s"%conditionToEvaluate)

            for n in range(0,len(sampleFit.R2)):
                evaluatedCondition = copy.copy(conditionToEvaluate)
                evaluatedCondition = evaluatedCondition.replace('$(R2)',"(%f)"%sampleFit.R2[n])
                evaluatedCondition = evaluatedCondition.replace('$(R2adj)',"(%f)"%sampleFit.R2adj[n])
                evaluatedCondition = evaluatedCondition.replace('$(AIC)',"(%f)"%sampleFit.AIC[n])
                evaluatedCondition = evaluatedCondition.replace('$(AICc)',"(%f)"%sampleFit.AICc[n])
                evaluatedCondition = evaluatedCondition.replace('$(BIC)',"(%f)"%sampleFit.BIC[n])

                for j in range(len(self.population.modelParameters)):
                    evaluatedCondition = evaluatedCondition.replace('$(%s)'%self.population.modelParameters[j],"(%f)"%sampleFit.parameters[n,j])
                evaluatedCondition = eval(evaluatedCondition)
                if (filterType==0 and not evaluatedCondition) or (filterType>=1 and evaluatedCondition):
                    newSampleFit.parameters = np.vstack([newSampleFit.parameters, sampleFit.parameters[n,:]])
                    newSampleFit.xB += sampleFit.xB
                    newSampleFit.yB += sampleFit.yB
                    newSampleFit.R2.append(sampleFit.R2[n])
                    newSampleFit.R2adj.append(sampleFit.R2adj[n])
                    newSampleFit.AIC.append(sampleFit.AIC[n])
                    newSampleFit.AICc.append(sampleFit.AICc[n])
                    newSampleFit.BIC.append(sampleFit.BIC[n])

        self.fitting.sampleFits.append(newSampleFit)
        self.fitting.write(self._getPath("bootstrapPopulation.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputPopulation=self.fitting)
        self._defineSourceRelation(self.inputPopulation, self.fitting)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        if self.filterType.get()==0:
            filterType = "Exclude"
        elif self.filterType.get()==1:
            filterType = "Keep"
        elif self.filterType.get()==2:
            filterType = "Keep confidence interval"
        msg=["Filter type: %s"%filterType]
        msg.append("Condition: %s"%self.condition.get())
        return msg

    def _validate(self):
        msg=[]
        if not self.inputPopulation.get().fnFitting.get().endswith("bootstrapPopulation.pkpd"):
            msg.append("Population must be a bootstrap sample")
        return msg
