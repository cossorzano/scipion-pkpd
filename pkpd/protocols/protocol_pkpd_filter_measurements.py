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

# TESTED in test_workflow_gabrielsson_pk02.py
# TESTED in test_workflow_gabrielsson_pk04.py
# TESTED in test_workflow_gabrielsson_pk06.py
# TESTED in test_workflow_gabrielsson_pk07.py
# TESTED in test_workflow_gabrielsson_pk11.py

class ProtPKPDFilterMeasurements(ProtPKPD):
    """ Filter measurements.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'filter measurements'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('filterType', params.EnumParam, choices=["Exclude","Keep","Remove NA","Remove ULOQ & LLOQ","Substitute LLOQ","Substitute ULOQ"],
                      label="Filter mode", default=0,
                      help='Exclude or keep measurements meeting the following condition.\n"\
                           "NA values are excluded in keep filters, and kept in exclude filters')
        form.addParam('condition', params.StringParam, label="Condition", condition="filterType<=1",
                      help='Example: $(t)<200\n $(Cp)>=1000 and $(Cp)<=2000"')
        form.addParam('substitute', params.FloatParam, label="Substitute by", condition="filterType==4 or filterType==5")

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
        elif self.filterType.get()==2:
            filterType="rmNA"
        elif self.filterType.get()==3:
            filterType="rmLL"
        elif self.filterType.get()==4:
            filterType="subsLL"
        elif self.filterType.get()==5:
            filterType="subsUL"

        filteredExperiment = PKPDExperiment()
        filteredExperiment.general = copy.copy(experiment.general)
        filteredExperiment.variables = copy.copy(experiment.variables)
        filteredExperiment.samples = {}
        filteredExperiment.doses = {}
        filteredExperiment.vias = {}

        usedDoses = []
        for sampleKey, sample in experiment.samples.iteritems():
            candidateSample = PKPDSample()
            candidateSample.variableDictPtr    = copy.copy(sample.variableDictPtr)
            candidateSample.doseDictPtr        = copy.copy(sample.doseDictPtr)
            candidateSample.sampleName         = copy.copy(sample.sampleName)
            candidateSample.doseList           = copy.copy(sample.doseList)
            candidateSample.descriptors        = copy.copy(sample.descriptors)
            candidateSample.measurementPattern = copy.copy(sample.measurementPattern)

            N = 0 # Number of initial measurements
            if len(sample.measurementPattern)>0:
                aux=getattr(sample,"measurement_%s"%sample.measurementPattern[0])
                N = len(aux)
            if N==0:
                continue

            # Create empty output variables
            Nvar = len(sample.measurementPattern)
            convertToFloat = []
            for i in range(0,Nvar):
                exec("candidateSample.measurement_%s = []"%sample.measurementPattern[i])
                convertToFloat.append(sample.variableDictPtr[sample.measurementPattern[i]].varType == PKPDVariable.TYPE_NUMERIC)

            for n in range(0,N):
                toAdd = []
                okToAddTimePoint = True
                conditionPython = copy.copy(condition)
                for i in range(0,Nvar):
                    exec("aux=sample.measurement_%s[%d]"%(sample.measurementPattern[i],n))
                    if filterType=="rmNA":
                        if aux=="NA" or aux=="None":
                            okToAddTimePoint = False
                        else:
                            toAdd.append(aux)
                    elif filterType=="rmLL":
                        if aux=="LLOQ" or aux=="ULOQ":
                            okToAddTimePoint = False
                        else:
                            toAdd.append(aux)
                    elif filterType=="subsLL":
                        okToAddTimePoint = True
                        if aux=="LLOQ":
                            toAdd.append(str(self.substitute.get()))
                        else:
                            toAdd.append(aux)
                    elif filterType=="subsUL":
                        okToAddTimePoint = True
                        if aux=="ULOQ":
                            toAdd.append(str(self.substitute.get()))
                        else:
                            toAdd.append(aux)
                    else:
                        # Keep or exclude
                        toAdd.append(aux)
                        varString = "$(%s)"%sample.measurementPattern[i]
                        if varString in conditionPython:
                            if aux=="NA":
                                okToAddTimePoint = False
                            else:
                                if convertToFloat[i]:
                                    conditionPython = conditionPython.replace(varString,"%f"%float(aux))
                                else:
                                    conditionPython = conditionPython.replace(varString,"'%s'"%aux)
                if (filterType=="exclude" or filterType=="keep") and okToAddTimePoint:
                    okToAddTimePoint = eval(conditionPython, {"__builtins__" : None }, {})
                    if filterType=="exclude":
                        okToAddTimePoint = not okToAddTimePoint
                if okToAddTimePoint:
                    for i in range(0,Nvar):
                        exec("candidateSample.measurement_%s.append('%s')"%(sample.measurementPattern[i],toAdd[i]))

            N = len(getattr(sample,"measurement_%s"%sample.measurementPattern[0])) # Number of final measurements
            if N!=0:
                filteredExperiment.samples[candidateSample.sampleName] = candidateSample
                for doseName in candidateSample.doseList:
                    if not doseName in usedDoses:
                        usedDoses.append(doseName)

        if len(usedDoses)>0:
            for doseName in usedDoses:
                filteredExperiment.doses[doseName] = copy.copy(experiment.doses[doseName])
                viaName=experiment.doses[doseName].via.viaName
                if not viaName in filteredExperiment.vias:
                    filteredExperiment.vias[viaName] = copy.copy(experiment.vias[viaName])

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
