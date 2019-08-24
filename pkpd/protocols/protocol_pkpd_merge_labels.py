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
from pkpd.objects import PKPDExperiment, PKPDVariable
from pkpd.pkpd_units import PKPDUnit


class ProtPKPDMergeLabels(ProtPKPD):
    """ Merge the labels of Experiment 2 into the samples of Experiment 1.\n
        If a label is in Experiment 1 and in Experiment 2, it remains the value from Experiment 1.
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'merge labels'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment1', params.PointerParam, label="Experiment 1",
                      pointerClass='PKPDExperiment',
                      help='Labels from Experiment 2 will be merged into Experiment 1. If a label already exists in Experiment 1, this has preference.')
        form.addParam('inputExperiment2', params.PointerParam, label="Experiment 2",
                      pointerClass='PKPDExperiment',
                      help='Labels from Experiment 2 will be merged into Experiment 1. If a label already exists in Experiment 1, this has preference.')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runMerge',self.inputExperiment1.get().getObjId(), self.inputExperiment2.getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runMerge(self, objId1, objId2):
        self.experiment = self.readExperiment(self.inputExperiment1.get().fnPKPD)
        self.experiment2 = self.readExperiment(self.inputExperiment2.get().fnPKPD)
        self.printSection("Merging labels")

        labelsToAdd = []
        for varName, variable in self.experiment2.variables.iteritems():
            if variable.role == PKPDVariable.ROLE_LABEL:
                if not varName in self.experiment.variables:
                    labelsToAdd.append(variable)

        for sampleName, sample in self.experiment.samples.iteritems():
            if sampleName in self.experiment2.samples:
                sample2 = self.experiment2.samples[sampleName]
                for variable in labelsToAdd:
                    varValue = sample2.descriptors[variable.varName]
                    self.experiment.addParameterToSample(sampleName, variable.varName, variable.units.unit, variable.comment,
                                                         varValue)

        self.writeExperiment(self.experiment,self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment1, self.experiment)
        self._defineSourceRelation(self.inputExperiment2, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=["Labels from Experiment %s were merged into %s"%(self.getObjectTag(self.inputExperiment1.get()),
                                                              self.getObjectTag(self.inputExperiment2.get()))]
        return msg
