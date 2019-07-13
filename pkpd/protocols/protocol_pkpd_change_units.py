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
from pkpd.pkpd_units import  unitFromString, convertUnits, strUnit, PKPDUnit

# TESTED in test_workflow_gabrielsson_pk01.py
# TESTED in test_workflow_gabrielsson_pk03.py
# TESTED in test_workflow_gabrielsson_pk04.py
# TESTED in test_workflow_gabrielsson_pk05.py
# TESTED in test_workflow_gabrielsson_pk06.py
# TESTED in test_workflow_gabrielsson_pk08.py
# TESTED in test_workflow_gabrielsson_pk11.py
# TESTED in test_workflow_gabrielsson_pk12.py
# TESTED in test_workflow_gabrielsson_pk14.py
# TESTED in test_workflow_gabrielsson_pk16.py
# TESTED in test_workflow_gabrielsson_pk20.py
# TESTED in test_workflow_gabrielsson_pk25.py
# TESTED in test_workflow_gabrielsson_pk39.py
# TESTED in test_workflow_gabrielsson_pk43.py
# TESTED in test_workflow_levyplot


class ProtPKPDChangeUnits(ProtPKPD):
    """ Change units of a given variable.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'change units'

    choicesTime=["h","min","sec"]
    choicesInvTime=["1/h","1/min","1/sec"]
    choicesWeight=["kg","g","mg","ug","ng"]
    choicesVolume=["L","mL","uL","nL"]
    choicesConc=["g/L","mg/L","ug/L","ng/L","g/mL","mg/mL","ug/mL","g/uL"]
    choicesAUC=["g*h/L","mg*h/L","ug*h/L","ng*h/L","g*h/mL","g*h/uL","g*min/L","mg*min/L",\
                 "ug*min/L","ng*min/L","g*min/mL","g*min/uL"]
    choicesAUMC=["mg*h^2/L","ug*h^2/L","ng*h^2/L","g*h^2/mL","g*h^2/uL","g*min^2/L","mg*min^2/L",\
                 "ug*min^2/L","ng*min^2/L","g*min^2/mL","g*min^2/uL"]
    choicesCl=["L/h","mL/h","uL/h","nL/h","L/min","mL/min","uL/min","nL/min","L/s","mL/s","uL/s","nL/s"]
    choicesVnorm=["L/kg","L/g"]
    choicesWeightInvTime=["kg/h","g/h","mg/h","ug/h","ng/h","kg/min","g/min","mg/min",\
                  "ug/min","ng/min","kg/s","g/s","mg/s","ug/s","ng/s"]

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('labelToChange', params.StringParam, label="Label to change", default="",
                      help='Name of the variable to change its units')
        form.addParam('newUnitsCategory', params.EnumParam, choices=["Time","1/Time","Weight","Volume",\
                                                                     "Weight/Volume","Weight*Time/Volume",\
                                                                     "Weight*Time^2/Volume","Volume/Time",\
                                                                     "Volume/Weight","Weight/Time"], label="Unit category", default=0)
        form.addParam('newUnitsCategoryTime', params.EnumParam, choices=ProtPKPDChangeUnits.choicesTime, label="Change time to",
                      condition="newUnitsCategory==0", default=1)
        form.addParam('newUnitsCategoryInvTime', params.EnumParam, choices=ProtPKPDChangeUnits.choicesInvTime, label="Change 1/time to",
                      condition="newUnitsCategory==1", default=1)
        form.addParam('newUnitsCategoryWeight', params.EnumParam, choices=ProtPKPDChangeUnits.choicesWeight, label="Change weight to",
                      condition="newUnitsCategory==2", default=1)
        form.addParam('newUnitsCategoryVolume', params.EnumParam, choices=ProtPKPDChangeUnits.choicesVolume, label="Change volume to",
                      condition="newUnitsCategory==3", default=1)
        form.addParam('newUnitsCategoryConc', params.EnumParam, choices=ProtPKPDChangeUnits.choicesConc, label="Change weight/volume [concentration] to",
                      condition="newUnitsCategory==4", default=1)
        form.addParam('newUnitsCategoryAUC', params.EnumParam, choices=ProtPKPDChangeUnits.choicesAUC, label="Change weight*time/volume [AUC] to",
                      condition="newUnitsCategory==5", default=1)
        form.addParam('newUnitsCategoryAUMC', params.EnumParam, choices=ProtPKPDChangeUnits.choicesAUMC, label="Change weight*time^2/volume [AUMC] to",
                      condition="newUnitsCategory==6", default=1)
        form.addParam('newUnitsCategoryCl', params.EnumParam, choices=ProtPKPDChangeUnits.choicesCl, label="Change volume/time [clearance] to",
                      condition="newUnitsCategory==7", default=1)
        form.addParam('newUnitsCategoryVnorm', params.EnumParam, choices=ProtPKPDChangeUnits.choicesVnorm, label="Change volume/weight [normalized vol] to",
                      condition="newUnitsCategory==8", default=1)
        form.addParam('newUnitsCategoryWeightInvTime', params.EnumParam, choices=ProtPKPDChangeUnits.choicesWeightInvTime, label="Change weight/time to",
                      condition="newUnitsCategory==9", default=1)


    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runChange',self.inputExperiment.get().getObjId(), self.labelToChange.get(),self.newUnitsCategory.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runChange(self, objId, labelToChange, category):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        currentUnit = self.experiment.getVarUnits(self.labelToChange.get())
        newUnit = unitFromString(self._getNewUnit())
        K=convertUnits(1.0,currentUnit,newUnit)

        variable = self.experiment.variables[self.labelToChange.get()]

        for sampleName, sample in self.experiment.samples.iteritems():
            if variable.varName == "dose":
                pass
            else:
                if variable.role == PKPDVariable.ROLE_LABEL:
                    varValue = float(sample.descriptors[variable.varName])
                    sample.descriptors[variable.varName] = K*varValue
                elif variable.role == PKPDVariable.ROLE_MEASUREMENT:
                    newValues = []
                    for x in sample.getValues(variable.varName):
                        if x=="NA" or x=="LLOQ" or x=="ULOQ" or x=="None":
                            newValues.append(x)
                        else:
                            newValues.append(str(K*float(x)))
                    sample.setValues(variable.varName,newValues)
                elif variable.role == PKPDVariable.ROLE_TIME:
                    newValues = []
                    for x in sample.getValues(variable.varName):
                        newValues.append(str(K*float(x)))
                    sample.setValues(variable.varName,newValues)
        variable.units = PKPDUnit()
        variable.units.unit = newUnit

        self.writeExperiment(self.experiment,self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _getNewUnit(self):
        newUnit = "none"
        if self.newUnitsCategory==0:
            newUnit = ProtPKPDChangeUnits.choicesTime[self.newUnitsCategoryTime.get()]
        elif self.newUnitsCategory==1:
            newUnit = ProtPKPDChangeUnits.choicesInvTime[self.newUnitsCategoryInvTime.get()]
        elif self.newUnitsCategory==2:
            newUnit = ProtPKPDChangeUnits.choicesWeight[self.newUnitsCategoryWeight.get()]
        elif self.newUnitsCategory==3:
            newUnit = ProtPKPDChangeUnits.choicesVolume[self.newUnitsCategoryVolume.get()]
        elif self.newUnitsCategory==4:
            newUnit = ProtPKPDChangeUnits.choicesConc[self.newUnitsCategoryConc.get()]
        elif self.newUnitsCategory==5:
            newUnit = ProtPKPDChangeUnits.choicesAUC[self.newUnitsCategoryAUC.get()]
        elif self.newUnitsCategory==6:
            newUnit = ProtPKPDChangeUnits.choicesAUMC[self.newUnitsCategoryAUMC.get()]
        elif self.newUnitsCategory==7:
            newUnit = ProtPKPDChangeUnits.choicesCl[self.newUnitsCategoryCl.get()]
        elif self.newUnitsCategory==8:
            newUnit = ProtPKPDChangeUnits.choicesVnorm[self.newUnitsCategoryVnorm.get()]
        elif self.newUnitsCategory==9:
            newUnit = ProtPKPDChangeUnits.choicesWeightInvTime[self.newUnitsCategoryWeightInvTime.get()]
        return newUnit

    def _summary(self):
        msg=["%s changed to %s"%(self.labelToChange.get(),self._getNewUnit())]
        return msg

    def _validate(self):
        errors=[]
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD, False)
        if not self.labelToChange.get() in experiment.variables:
            errors.append("Cannot find %s as variable"%self.labelToChange)
        else:
            variable = experiment.variables[self.labelToChange.get()]
            if variable.varType!=PKPDVariable.TYPE_NUMERIC:
                errors.append("%s is not a numeric variable"%variable.varName)
            else:
                currentUnit = experiment.getVarUnits(self.labelToChange.get())
                newUnit = unitFromString(self._getNewUnit())
                try:
                    K=convertUnits(1,currentUnit,newUnit)
                    if K==None:
                        errors.append("Unknown conversion from %s to %s. If it makes sense and it is not implemented you may contact info@kinestat.com"%\
                                      (strUnit(currentUnit),strUnit(newUnit)))
                except:
                    errors.append("Unknown conversion from %s to %s. If it makes sense and it is not implemented you may contact info@kinestat.com"%(currentUnit,newUnit))
        return errors

    def filterVarForWizard(self, v):
        """ Define the type of variables required (used in wizard). """
        return v.isNumeric()
