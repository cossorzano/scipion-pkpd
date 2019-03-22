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
from pkpd.objects import PKPDExperiment
from pkpd.pkpd_units import PKPDUnit


class ProtPKPDCreateLabel2Exps(ProtPKPD):
    """ Create label by performing calculations on already existing labels from two different experiments.\n
        The protocol assumes that the same samples are present in both experiments and calculations are performed
        only on those samples with the same name in both experiments.
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'create label 2 Exps'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment1', params.PointerParam, label="Input experiment1", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('inputExperiment2', params.PointerParam, label="Input experiment2", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('labelToAdd', params.StringParam, label="Label to add", default="",
                      help='Name of the variable to add')
        form.addParam('expression', params.StringParam, label="Expression to calculate", default="",
                      help='For example, to estimate the Mean Absorption Time you may subtract the MRT of an intravenous experiment from the MRT of an oral experiment: $Exp1(MRT)-$Exp2(MRT)')
        form.addParam('units', params.StringParam, label="Units", default="None",
                      help='For example, L/kg')
        form.addParam('comment', params.StringParam, label="Label comment", default="",
                      help='For example, apparent volume of distribution per kilogram')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runCreate',self.inputExperiment1.get().getObjId(), self.inputExperiment2.get().getObjId(),
                                 self.labelToAdd.get(), self.expression.get(), self.comment.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runCreate(self, objId1, objId2, labelToAdd, expression, comment):
        self.experiment1 = self.readExperiment(self.inputExperiment1.get().fnPKPD)
        self.experiment2 = self.readExperiment(self.inputExperiment2.get().fnPKPD)

        labelToAdd = self.labelToAdd.get().replace(' ',"_")
        units = PKPDUnit(self.units.get())
        for sampleName, sample1 in self.experiment1.samples.iteritems():
            if sampleName in self.experiment2.samples:
                expression1 = sample1.substituteValuesInExpression(self.expression.get(),"Exp1")
                sample2 = self.experiment2.samples[sampleName]
                expression2 = sample2.substituteValuesInExpression(expression1,"Exp2")
                varValue = eval(expression2, {"__builtins__" : {"True": True, "False": False, "None": None} }, {})
                self.experiment1.addParameterToSample(sampleName, labelToAdd, units.unit, self.comment.get(), varValue)
            else:
                # Remove this sample from the experiment because we cannot perform this calculation
                self.experiment1.samples.pop(sampleName)

        self.writeExperiment(self.experiment1,self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment1)
        self._defineSourceRelation(self.inputExperiment1, self.experiment1)
        self._defineSourceRelation(self.inputExperiment2, self.experiment1)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=["%s created as %s"%(self.labelToAdd.get(),self.expression.get())]
        return msg
