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
from itertools import izip_longest


class ProtPKPDCreateLabel(ProtPKPD):
    """ Create label by performing calculations on already existing labels.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'create label'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('labelToAdd', params.StringParam, label="Label(s) to add", default="",
                      help='Name of the variable to add. If several names are given, separated by semicolons.')
        form.addParam('rewrite', params.BooleanParam, label="Rewrite labels", default=False,
                      help='Set this flag to true if you want to rewrite the content of existing labels.')
        form.addParam('expression', params.StringParam, label="Expression(s) to calculate", default="",
                      help='For example, to normalize the apparent volume of distribution by the animal weight use $(Vd)/$(weight), a literal as "T1", or a constant 1.5. '\
                           'If several labels are created, separated by semicolons.')
        form.addParam('units', params.StringParam, label="Units", default="None",
                      help='For example, L/kg. If several labels are created, separated by semicolons. '\
                           'Set to None for no units.')
        form.addParam('comment', params.StringParam, label="Label comment(s)", default="",
                      help='For example, apparent volume of distribution per kilogram. '\
                            'If several labels are created, separated by semicolons.')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runCreate',self.inputExperiment.get().getObjId(), self.labelToAdd.get(),
                                 self.expression.get(), self.comment.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runCreate(self, objId, labelToAdd, expression, comment):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        labels = self.labelToAdd.get().split(';')
        expressions = self.expression.get().split(';')
        units = self.units.get().split(';')
        comments = self.comment.get().split(';')

        for label, expression, unit, comment in izip_longest(labels,expressions,units,comments,fillvalue=""):
            labelToAdd = label.strip().replace(' ',"_")
            units = PKPDUnit(unit.strip())
            for sampleName, sample in self.experiment.samples.iteritems():
                varValue = sample.evaluateExpression(expression.strip())
                self.experiment.addParameterToSample(sampleName, labelToAdd, units.unit, comment.strip(), varValue,
                                                     self.rewrite.get())

        self.writeExperiment(self.experiment,self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        labels = self.labelToAdd.get().split(';')
        expressions = self.expression.get().split(';')
        for label, expression in izip_longest(labels,expressions,fillvalue=""):
            msg.append("%s created as %s"%(label.strip(),expression.strip()))
        return msg
