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
import sys

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.data import PKPDExperiment
from pyworkflow.protocol.constants import LEVEL_ADVANCED

class ProtPKPDJoinSamples(ProtPKPD):
    """ Join samples.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'join samples'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment1', params.PointerParam, label="Input experiment 1",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('prefix1', params.StringParam, label="Prefix 1", expertLevel=LEVEL_ADVANCED,
                      help='Select prefix like Exp1_')
        form.addParam('inputExperiment2', params.PointerParam, label="Input experiment 2",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('prefix2', params.StringParam, label="Prefix 2", expertLevel=LEVEL_ADVANCED,
                      help='Select prefix like Exp2_')
        form.addParam('title', params.StringParam, label="Join experiment title",
                      help='If empty, the title of the first experiment is taken')
        form.addParam('comment', params.StringParam, label="Join experiment comment",
                      help='If empty, the comment of the first experiment is taken')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runJoin',self.inputExperiment1.get().getObjId(), self.inputExperiment1.get().getObjId(), \
                                 self.prefix1.get(), self.prefix2.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runJoin(self, objId1, objId2, prefix1, prefix2):
        experiment1 = self.readExperiment(self.inputExperiment1.get().fnPKPD)
        experiment2 = self.readExperiment(self.inputExperiment2.get().fnPKPD)
        self.printSection("Joining")

        self.experiment = PKPDExperiment()

        # General
        if self.title!=None:
            self.experiment.general["title"]=self.title.get()
        if self.comment!=None:
            self.experiment.general["comment"]=self.comment.get()
        for key, value in experiment1.general.iteritems():
            if not (key in self.experiment.general):
                self.experiment.general[key] = copy.copy(value)
        for key, value in experiment2.general.iteritems():
            if not (key in self.experiment.general):
                self.experiment.general[key] = copy.copy(value)

        # Variables
        for key, value in experiment1.variables.iteritems():
            if not (key in self.experiment.variables):
                self.experiment.variables[key] = copy.copy(value)
        for key, value in experiment2.variables.iteritems():
            if not (key in self.experiment.variables):
                self.experiment.variables[key] = copy.copy(value)

        # Doses
        for key, value in experiment1.doses.iteritems():
            dose = copy.copy(value)
            dose.varName = "%s%s"%(self.prefix1.get(),key)
            self.experiment.doses[dose.varName] = dose
        for key, value in experiment2.doses.iteritems():
            dose = copy.copy(value)
            dose.varName = "%s%s"%(self.prefix2.get(),key)
            self.experiment.doses[dose.varName] = dose

        # Samples
        for key, value in experiment1.samples.iteritems():
            sample = copy.copy(value)
            sample.sampleName = "%s%s"%(self.prefix1.get(),key)
            sample.doseList = ["%s%s"%(self.prefix1.get(),doseName) for doseName in sample.doseList]
            self.experiment.samples[sample.sampleName] = sample
        for key, value in experiment2.samples.iteritems():
            sample = copy.copy(value)
            sample.sampleName = "%s%s"%(self.prefix2.get(),key)
            sample.doseList = ["%s%s"%(self.prefix2.get(),doseName) for doseName in sample.doseList]
            self.experiment.samples[sample.sampleName] = sample

        # Print and save
        self.writeExperiment(self.experiment,self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment1, self.experiment)
        self._defineSourceRelation(self.inputExperiment2, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        import re
        errors = []
        if self.prefix1.get()==self.prefix2.get():
            experiment1 = PKPDExperiment()
            experiment1.load(self.inputExperiment1.get().fnPKPD)
            experiment2 = PKPDExperiment()
            experiment2.load(self.inputExperiment2.get().fnPKPD)

            # Check if there are repeated doses
            for doseName1 in experiment1.doses:
                if doseName1 in experiment2.doses:
                    errors.append("Dose %s is repeated in both experiments"%doseName1)

            # Check if there are repeated samples
            for sampleName1 in experiment1.samples:
                if sampleName1 in experiment2.samples:
                    errors.append("Sample %s is repeated in both experiments"%sampleName1)
        # if self.prefix1.get()!="" and not re.match("[_A-Za-z][_a-zA-Z0-9]*$",self.prefix1.get()):
        #     errors.append("Prefix1 is not well formatted")
        # if self.prefix2.get()!="" and not re.match("[_A-Za-z][_a-zA-Z0-9]*$",self.prefix2.get()):
        #     errors.append("Prefix2 is not well formatted")
        if len(errors)>0:
            errors.append("Use the prefixes in the Advanced options")
        return errors
