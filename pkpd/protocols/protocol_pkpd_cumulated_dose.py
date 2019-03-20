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


class ProtPKPDCumulatedDose(ProtPKPD):
    """ Create label with the cumulated dose between two time points.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'cumulated dose'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('labelToAdd', params.StringParam, label="Label to add", default="CumulatedDose",
                      help='Name of the variable to add')
        form.addParam('fromTime', params.FloatParam, label="From (h)", default=0)
        form.addParam('toTime', params.FloatParam, label="To (h)", default=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runCreate',self.inputExperiment.get().getObjId(), self.labelToAdd.get(),
                                 self.fromTime.get(), self.toTime.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runCreate(self, objId, labelToAdd, fromTime, toTime):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        labelToAdd = self.labelToAdd.get().replace(' ',"_")
        fromTime*=60
        toTime*=60
        comment = "Cumulated dose between %s and %s"%(fromTime,toTime)
        for sampleName, sample in self.experiment.samples.iteritems():
            sample.interpretDose()
            Dunits = sample.getDoseUnits()
            varValue = sample.getCumulatedDose(fromTime,toTime)
            self.experiment.addParameterToSample(sampleName, labelToAdd, Dunits, comment, varValue)

        self.writeExperiment(self.experiment,self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=["Cumulated dose between %s and %s"%(self.fromTime.get(),self.toTime.get())]
        return msg
