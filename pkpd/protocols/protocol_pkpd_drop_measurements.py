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
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable


class ProtPKPDDropMeasurements(ProtPKPD):
    """ Filter measurements.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'drop measurements'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('varsToDrop', params.StringParam, label="Variables to drop", default="",
                      help='List of variable names to drop separated by commas')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runDrop',self.inputExperiment.get().getObjId(), self.varsToDrop.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runDrop(self, objId, varsToDrop):
        import copy
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        self.printSection("Dropping variables")
        varsToDrop = []
        for varName in self.varsToDrop.get().split(','):
            varsToDrop.append(varName.strip())

        filteredExperiment = PKPDExperiment()
        filteredExperiment.general = copy.copy(experiment.general)
        filteredExperiment.variables = {}
        for varName, variable in experiment.variables.iteritems():
            if not varName in varsToDrop:
                filteredExperiment.variables[varName] = copy.copy(variable)
        filteredExperiment.samples = {}
        filteredExperiment.doses = copy.copy(experiment.doses)

        for sampleKey, sample in experiment.samples.iteritems():
            candidateSample = PKPDSample()
            candidateSample.variableDictPtr    = filteredExperiment.variables
            candidateSample.doseDictPtr        = filteredExperiment.doses
            candidateSample.sampleName          = copy.copy(sample.sampleName)
            candidateSample.doseList           = copy.copy(sample.doseList)
            candidateSample.descriptors        = copy.copy(sample.descriptors)
            candidateSample.measurementPattern = []
            for varName in sample.measurementPattern:
                if not varName in varsToDrop:
                    candidateSample.measurementPattern.append(varName)
                    exec("candidateSample.measurement_%s = copy.copy(sample.measurement_%s)"%(varName,varName))
            filteredExperiment.samples[candidateSample.varName] = candidateSample

        self.writeExperiment(filteredExperiment,self._getPath("experiment.pkpd"))
        self.experiment = filteredExperiment

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=["Variables %s dropped"%self.varsToDrop.get()]
        return msg

    def filterVarForWizard(self, v):
        """ Define the type of variables required (used in wizard). """
        return v.isMeasurement()
