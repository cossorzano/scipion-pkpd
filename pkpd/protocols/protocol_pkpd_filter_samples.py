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

import sys

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDExperiment, PKPDVariable


class ProtPKPDFilterSamples(ProtPKPD):
    """ Filter samples.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'filter samples'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('filterType', params.EnumParam, choices=["Exclude","Keep","Remove NA"], label="Filter mode", default=0,
                      help='Exclude or keep samples meeting the following condition\n'\
                           "NA values are excluded in keep filters, and kept in exclude filters")
        form.addParam('condition', params.StringParam, label="Condition", condition="filterType!=2",
                      help='Example: $(weight)<200 and $(sex)=="female"\n'
                           'You may use any of the variables of the experiment or the variable $(sampleName)')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runFilter',self.inputExperiment.get().getObjId(), self.filterType.get(), \
                                 self.condition.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runFilter(self, objId, filterType, condition):
        import copy
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        self.printSection("Filtering")
        if self.filterType.get()==0:
            filterType="exclude"
        elif self.filterType.get()==1:
            filterType="keep"
        else:
            filterType="rmNA"

        filteredExperiment = PKPDExperiment()
        filteredExperiment.general = copy.copy(experiment.general)
        filteredExperiment.variables = copy.copy(experiment.variables)
        filteredExperiment.groups = copy.copy(experiment.groups)
        filteredExperiment.samples = {}
        filteredExperiment.doses = {}

        # http://stackoverflow.com/questions/701802/how-do-i-execute-a-string-containing-python-code-in-python
        safe_list = ['descriptors']
        safe_dict = dict([ (k, locals().get(k, None)) for k in safe_list ])
        usedDoses = []
        for sampleKey, sample in experiment.samples.iteritems():
            ok = filterType=="rmNA"
            try:
                if filterType == "rmNA":
                    conditionPython = "True"
                else:
                    conditionPython = copy.copy(condition)
                conditionPython = conditionPython.replace('$(sampleName)','"%s"'%sample.sampleName)
                for key, variable in experiment.variables.iteritems():
                    if key in sample.descriptors:
                        value = sample.descriptors[key]
                        if value=="NA":
                            conditionPython="False"
                        else:
                            if filterType!="rmNA":
                                if variable.varType == PKPDVariable.TYPE_NUMERIC:
                                    conditionPython = conditionPython.replace("$(%s)"%key,"%f"%float(sample.descriptors[key]))
                                else:
                                    conditionPython = conditionPython.replace("$(%s)"%key,"'%s'"%sample.descriptors[key])
                ok=eval(conditionPython, {"__builtins__" : {"True": True, "False": False} }, {})
            except:
                print sys.exc_info()[0]
                pass
            if (ok and (filterType=="keep" or filterType=="rmNA")) or (not ok and filterType=="exclude"):
                filteredExperiment.samples[sampleKey] = copy.copy(sample)
                for doseName in sample.doseList:
                    usedDoses.append(doseName)

        if len(usedDoses)>0:
            for doseName in usedDoses:
                filteredExperiment.doses[doseName] = copy.copy(experiment.doses[doseName])

        self.writeExperiment(filteredExperiment,self._getPath("experiment.pkpd"))
        self.experiment = filteredExperiment

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        if self.filterType.get()==0:
            msg.append("Exclude %s"%self.condition.get())
        elif self.filterType.get()==1:
            msg.append("Keep %s"%self.condition.get())
        elif self.filterType.get()==2:
            msg.append("Remove NA")
        return msg