# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Carlos Oscar Sorzano (info@kinestat.com)
# *              
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os
import Tkinter as tk
import ttk
import pyworkflow.object as pwobj
from pyworkflow.wizard import Wizard
import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.tree import TreeProvider, BoundTree
from pyworkflow.gui.dialog import ListDialog
import pyworkflow.gui as gui

from pkpd.protocols import *
from pkpd.viewers.tk_ode import PKPDODEDialog, PKPDFitDialog


class VariablesProvider(dialog.Dialog):
    """
    This class create a frame with the pkpd variables
    """
    def __init__(self, parent, title, **kwargs):
        """ From kwargs:
                message: message tooltip to show when browsing.
                selected: the item that should be selected.
                validateSelectionCallback: a callback function to validate selected items.
        """
        self.targetProtocol = kwargs['targetProtocol']
        self.params = kwargs['params']
        self.values = []
        self.oldValues = {}
        self.experiment = kwargs['experiment']
        self.provider = kwargs['provider']

        dialog.Dialog.__init__(self, parent, title,
                               buttons=[('Select', dialog.RESULT_YES),
                                        ('Cancel', dialog.RESULT_CANCEL)])

    def body(self, bodyFrame):
        bodyFrame.config(bg='white')
        gui.configureWeigths(bodyFrame)
        self._createContent(bodyFrame)

    def _createContent(self, content):
        """
        Create the wizard content
        """
        # Create the variables
        frame = tk.Frame(content, bg='white')
        variables = tk.LabelFrame(frame, text="Parameters", bg='white')

        variables.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        self.samplesTree = self._addBoundTree(variables,
                                              self.provider, 4)

        # Create the content
        parameteres = tk.LabelFrame(frame, text="Variables", bg='white')
        parameteres.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        self._insertParameters(parameteres)

        frame.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(frame)

    def _insertParameters(self, content):
        row = 0
        for param in self.params:
            paramValue = self._createVariablesValue(content, param)
            paramValue.grid(row=row, column=0, sticky='news', padx=5, pady=5)
            row += 1

    def _createVariablesValue(self, frame, paramName):
        """
        Create the wizard variables with the values
        """
        def _fillVariableValues(paramName, varValue, cbox):
            self._objects = self.provider.getObjects()
            index = 0
            valueIndex = -1
            valuesList = []
            for obj in self._objects:
                objDict = self.provider.getObjectInfo(obj)
                if objDict is not None:
                    key = objDict.get('key')
                    # text = objDict.get('text', key)
                    # parent = objDict.get('parent', None)
                    valuesList.append(key)
                    if key == varValue:
                        valueIndex = index
                    index += 1
            cbox['values'] = valuesList
            if valueIndex != -1:
                cbox.current(valueIndex)
            else:
                _selection_changed(None, paramName, valuesList[0])
                cbox.current(0)

        def _selection_changed(event, param, newValue):
            for value in self.values:
                if value[0].lower() == param:
                    value[1].set(newValue)

        paramValueFrame = tk.Frame(frame, bg='white')
        gui.configureWeigths(paramValueFrame)
        label = tk.Label(paramValueFrame, text=(str(self.targetProtocol.getParam(paramName).label) + ': '))
        label.grid(row=0, column=0, sticky='news')
        combobox = ttk.Combobox(paramValueFrame, name=paramName.lower(), state="readonly")
        combobox.grid(row=0, column=1, sticky='news')
        combobox.bind("<<ComboboxSelected>>", lambda event: _selection_changed(event,
                                                                               combobox._name,
                                                                               combobox['value'][combobox.current()]))
        varValue = getattr(self.targetProtocol, paramName, None)
        self.values.append((paramName, varValue))
        self.oldValues[paramName] = varValue.get()
        _fillVariableValues(paramName, varValue, combobox)

        return paramValueFrame

    def _addBoundTree(self, parent, provider, height):
        bt = BoundTree(parent, provider, height=height)
        bt.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(parent)
        return bt

    def restartValues(self, restartAll=False):
        if restartAll:
            for value in self.values:
                self.oldValues[value[0]] = value[1].get()
        else:
            for key, oldValue in self.oldValues.items():
                for value in self.values:
                    if value[0] == key:
                        value[1].set(pwobj.String(oldValue))


class FilterVariablesTreeProvider(TreeProvider):
    """ Simplified view of VariablesTreeProvider with less columns.
    Additionally, we can filter by a given function. """

    def __init__(self, experiment, filter=None):
        self.params = params
        self.experiment = experiment
        self.filter = filter or self._noFilter

    def _noFilter(self, v):
        return True

    def getColumns(self):
        return [('Name', 100), ('Unit', 100), ('Comment', 300)]

    def getObjects(self):
        sortedVars = []
        for key in sorted(self.experiment.variables.keys()):
            v = self.experiment.variables[key]
            if self.filter(v):
                sortedVars.append(v)
        return sortedVars

    def getObjectInfo(self, obj):
        key = obj.varName
        return {'key': key, 'text': key,
                'values': (obj.getUnitsString(),
                           obj.comment)}


class PKPDChooseVariableWizard(Wizard):
    _targets = [(ProtPKPDRegressionLabel, ['labelX', 'labelY']),
                (ProtPKPDChangeUnits, ['labelToChange']),
                (ProtPKPDStatsExp1Subgroups2Mean, ['labelToCompare']),
                (ProtPKPDExponentialFit, ['predictor', 'predicted']),
                (ProtPKPDDissolutionFit, ['predictor', 'predicted']),
                (ProtPKPDExportToCSV, ['tVar', 'xVar']),
                (ProtPKPDEliminationRate, ['predictor', 'predicted']),
                (ProtPKPDMonoCompartment, ['predictor', 'predicted']),
                (ProtPKPDMonoCompartmentConv, ['predictor', 'predicted']),
                (ProtPKPDMonoCompartmentUrine, ['predictor', 'predicted',
                                                'Au']),
                (ProtPKPDMonoCompartmentClint,['predictor', 'predicted']),
                (ProtPKPDMonoCompartmentLinkPD, ['predictor', 'predicted',
                                                 'E']),
                (ProtPKPDMonoCompartmentPD, ['predictor', 'predicted', 'E']),
                (ProtPKPDNCANumeric, ['xVar']),
                (ProtPKPDTwoCompartments, ['predictor', 'predicted']),
                (ProtPKPDTwoCompartmentsConv, ['predictor', 'predicted']),
                (ProtPKPDTwoCompartmentsUrine, ['predictor', 'predicted', 'Au']),
                (ProtPKPDTwoCompartmentsClint, ['predictor', 'predicted']),
                (ProtPKPDTwoCompartmentsClintMetabolite, ['predictor',
                                                          'predicted',
                                                          'metabolite']),
                (ProtPKPDTwoCompartmentsAutoinduction, ['predictor',
                                                        'predicted']),
                (ProtPKPDTwoCompartmentsBoth, ['predictor', 'predicted',
                                               'Cperipheral']),
                (ProtPKPDTwoCompartmentsBothPD, ['predictor', 'predicted',
                                                 'Cperipheral', 'E']),
                (ProtPKPDMonoCompartmentUrine, ['predicted']),
                (ProtPKPDSimulateGenericPD, ['predictor']),
                ]

    def show(self, form, *params):
        protocol = form.protocol
        fullParams = self._getAllParams(protocol)
        experiment = protocol.getAttributeValue('inputExperiment', None)

        if experiment is None:
            form.showError("Select input experiment first.")
        else:
            experiment.load()
            filterFunc = getattr(protocol, 'filterVarForWizard', None)
            provider = FilterVariablesTreeProvider(experiment)
            dlg = VariablesProvider(form.root, self.getTitle(),
                                    targetProtocol=protocol,
                                    params=fullParams,
                                    experiment=experiment,
                                    provider=provider)
            if dlg.resultYes():
                for label in fullParams:
                    for value in dlg.values:
                        if value[0] == label:
                            self.setFormValues(form, label, value[1])
                            break
                dlg.restartValues(restartAll=True)
            else:
                dlg.restartValues()

    def getTitle(self):
        return "Choose variable"

    def setFormValues(self, form, label, values):
        form.setVar(label, values)

    def getSelectMode(self):
        return "browse"

    def _getAllParams(self, protocol):
        """ Return the list of all target parameters associated with this
        protocol. This function is useful when more than one parameter is
        associated to the wizard in the same protocol. """

        for k, v in self._targets:
            if k.__name__ == protocol.getClassName():
                return v

        return []


class DoseTreeProvider(TreeProvider):
    def __init__(self, experiment):
        self.params = params
        self.experiment = experiment

    def getColumns(self):
        return [('Name', 100)]

    def getObjects(self):
        sortedDoses = []
        for key in sorted(self.experiment.doses.keys()):
            d = self.experiment.doses[key]
            sortedDoses.append(d)
        return sortedDoses

    def getObjectInfo(self, obj):
        key = obj.doseName
        return {'key': key, 'text': key}


class PKPDChooseDoseWizard(Wizard):
    _targets = [(ProtPKPDChangeVia, ['doseName', 'doseVia']),
                ]

    def show(self, form, *params):
        protocol = form.protocol
        fullParams = self._getAllParams(protocol)
        experiment = protocol.getAttributeValue('inputExperiment', None)

        if experiment is None:
            form.showError("Select input experiment first.")
        else:
            experiment.load()
            provider = DoseTreeProvider(experiment)

            dlg = VariablesProvider(form.root, "Choose dose name",
                                    targetProtocol=protocol,
                                    params=fullParams,
                                    experiment=experiment,
                                    provider=provider,
                                    selectmode=self.getSelectMode())

            if dlg.resultYes():
                for label in fullParams:
                    for value in dlg.values:
                        if value[0] == label:
                            self.setFormValues(form, label, value[1])
                            break
                dlg.restartValues(restartAll=True)
            else:
                dlg.restartValues()

    def setFormValues(self, form, label, values):
        form.setVar(label, values)

    def getSelectMode(self):
        return "browse"

    def _getAllParams(self, protocol):
        """ Return the list of all target parameters associated with this
        protocol. This function is useful when more than one parameter is
        associated to the wizard in the same protocol. """

        for k, v in self._targets:
            if k.__name__ == protocol.getClassName():
                return v

        return []


class PKPDChooseSeveralVariableWizard(PKPDChooseVariableWizard):
    _targets = [(ProtPKPDDropMeasurements, ['varsToDrop']),
                (ProtPKPDScaleToCommonDose, ['measurementsToChange'])]

    def getTitle(self):
        return "Choose variable(s)"

    def setFormValues(self, form, label, values):
        form.setVar(label, ', '.join([var.varName for var in values]))

    def getSelectMode(self):
        return "extended"


class SimpleListTreeProvider(TreeProvider):
    """ A simple TreeProvider over the elements of a string list """

    def __init__(self, strList, name='Name', width=100):
        self.objects = [pwobj.String(value=s) for s in strList]
        self.columns = [(name, width)]

    def getColumns(self):
        return self.columns

    def getObjects(self):
        return self.objects

    def getObjectInfo(self, obj):
        key = obj.get()
        return {'key': key, 'text': key,
                'values': ()}


class PKPDVariableTemplateWizard(Wizard):
    _targets = [(ProtPKPDImportFromText, ['variables'])]

    def show(self, form, *params):

        protocol = form.protocol
        labels = self._getAllParams(protocol)
        fnCSV = protocol.getAttributeValue('inputFile', "")
        if not os.path.exists(fnCSV):
            form.showError("Select a valid CSV input file first.")
        else:
            varNames = getVarNamesFromCSVfile(fnCSV)
            tp = SimpleListTreeProvider(varNames, name="Variables")
            dlg = dialog.ListDialog(form.root, "Choose variable(s)", tp,
                                    selectmode='extended')
            if dlg.resultYes():
                strToAdd = ""
                for value in dlg.values:
                    strToAdd += ("\n%s ; [Units/none] ; [numeric/text] ; "
                                 "[time/label/measurement] ; [Comment]"
                                 % (value.get()))
                if strToAdd != "":
                    for label in labels:
                        currentValue = protocol.getAttributeValue(label, "")
                        form.setVar(label, currentValue+strToAdd)

    def _getAllParams(self, protocol):
        """ Return the list of all target parameters associated with this
        protocol. This function is useful when more than one parameter is
        associated to the wizard in the same protocol. """

        for k, v in self._targets:
            if k.__name__ == protocol.getClassName():
                return v

        return []


class PKPDDoseTemplateWizard(Wizard):
    _targets = [(ProtPKPDImportFromText, ['doses']),
                (ProtPKPDODESimulate, ['doses'])
                ]

    def show(self, form, *params):
        label = params[0]
        protocol = form.protocol
        currentValue = protocol.getAttributeValue(label, "")
        template = ("\nInfusion0 ; via=Intravenous; infusion; t=0.5:0.75 h; "
                    "d=60*weight/1000 mg\n Bolus1 ; via=Oral; bolus; t=2 h; "
                    "d=100 mg\n Treatment ; via=Oral; repeated_bolus; "
                    "t=0:8:48 h; d=100 mg")
        form.setVar(label, currentValue+template)


class PKPDDosesToSamplesTemplateWizard(Wizard):
    _targets = [(ProtPKPDImportFromText, ['dosesToSamples'])
                ]
    def show(self, form, *params):
        label = params
        protocol = form.protocol
        fnCSV = protocol.getAttributeValue('inputFile', "")
        if not os.path.exists(fnCSV):
            form.showError("Select a valid CSV input file first.")
        else:
            doseNames = []
            for line in protocol.doses.get().replace('\n',';;').split(';;'):
                tokens = line.split(';')
                if len(tokens)==4:
                    doseNames.append(tokens[0].strip())

            sampleNamesInCSV = getSampleNamesFromCSVfile(fnCSV)
            currentValue = protocol.getAttributeValue(label, "")
            sampleNamesAssigned = []
            for line in currentValue.replace('\n',';;').split(';;'):
                tokens = line.split(';')
                if len(tokens)==2:
                    sampleNamesAssigned.append(tokens[0].strip())

            dlg = MultiListDialog(form.root, "Test",
                                         [SimpleListTreeProvider(list(set(sampleNamesInCSV)-set(sampleNamesAssigned)),
                                                                 name="Samples"),
                                          SimpleListTreeProvider(doseNames,
                                                                 name="Doses")],
                             selectmode='extended')
            if dlg.resultYes():
                sampleList = dlg.values[0]
                doseList = dlg.values[1]
                for sample in sampleList:
                    if currentValue!="":
                        currentValue+="\n"
                    currentValue+="%s; "%sample.get().strip()
                    doseNo = 1
                    for dose in doseList:
                        if doseNo!=1:
                            currentValue+=", %s"%dose.get()
                        else:
                            currentValue+="%s"%dose.get()
                        doseNo += 1
                if currentValue!="":
                    form.setVar(label, currentValue)


class PKPDODEWizard(Wizard):
    _targets = [(ProtPKPDMonoCompartment, ['bounds']),
                (ProtPKPDMonoCompartmentConv, ['bounds']),
                (ProtPKPDMonoCompartmentUrine, ['bounds']),
                (ProtPKPDMonoCompartmentClint, ['bounds']),
                (ProtPKPDMonoCompartmentLinkPD, ['bounds']),
                (ProtPKPDMonoCompartmentPD, ['bounds']),
                (ProtPKPDTwoCompartments, ['bounds']),
                (ProtPKPDTwoCompartmentsConv, ['bounds']),
                (ProtPKPDTwoCompartmentsUrine, ['bounds']),
                (ProtPKPDTwoCompartmentsBoth, ['bounds']),
                (ProtPKPDTwoCompartmentsBothPD, ['bounds']),
                (ProtPKPDTwoCompartmentsClint, ['bounds']),
                (ProtPKPDTwoCompartmentsClintMetabolite, ['bounds']),
                (ProtPKPDTwoCompartmentsAutoinduction, ['bounds']),
                (ProtPKPDGenericFit, ['bounds'])]

    _nonODE = [ProtPKPDGenericFit]

    def show(self, form, *params):
        protocol = form.protocol
        labels = self._getAllParams(protocol)
        experiment = protocol.getAttributeValue('inputExperiment', None)

        if experiment is None:
            form.showError("Select the input experiment first.")
        else:
            # try:
                protocol.setInputExperiment() # this load the experiment
                if not type(protocol) in PKPDODEWizard._nonODE:
                    protocol.clearGroupParameters()
                protocol.getXYvars()
                if not type(protocol) in PKPDODEWizard._nonODE:
                    protocol.configureSource(protocol.createDrugSource())
                protocol.setupModel()
                if not type(protocol) in PKPDODEWizard._nonODE:
                    protocol.model.drugSource = protocol.drugSource

                i = 0
                for sampleName, sample in protocol.experiment.samples.iteritems():
                    sample.interpretDose()
                    if i == 0:
                        if not type(protocol) in PKPDODEWizard._nonODE:
                            protocol.model.drugSource.setDoses(sample.parsedDoseList,0,1) # Needed to correctly identify the parameters of the source
                        protocol.model.setSample(sample)
                        protocol.calculateParameterUnits(sample)
                    i += 1

                if not type(protocol) in PKPDODEWizard._nonODE:
                    dlg = PKPDODEDialog(form.root, "Select Parameter Bounds",
                                        targetProtocol=protocol,
                                        varNameX=protocol.predictor.get(),
                                        varNameY=protocol.predicted.get())
                else:
                    dlg = PKPDFitDialog(form.root, "Select Parameter Bounds",
                                        targetProtocol=protocol,
                                        varNameX=protocol.predictor.get(),
                                        varNameY=protocol.predicted.get())

                if dlg.resultYes():
                    boundStr = ""
                    i = 1
                    for bound in dlg.getBoundsList():
                        if i > 1:
                            boundStr += "; "
                        boundStr += str(bound)
                        i += 1
                    if boundStr != "":
                        for label in labels:
                            form.setVar(label, boundStr)
            # except Exception as e:
                # pass
                # form.showError("Error: %s" % str(e))

    def _getAllParams(self, protocol):
        """ Return the list of all target parameters associated with this
        protocol. This function is useful when more than one parameter is
        associated to the wizard in the same protocol. """

        for k, v in self._targets:
            if k.__name__ == protocol.getClassName():
                return v

        return []


class MultiListDialog(ListDialog):
    """
    Select elements among several lists.
    """
    def _createTree(self, parent):
        treeFrame = tk.Frame(parent)
        treeFrame.grid(row=0, column=0, sticky='news')
        treeFrame.rowconfigure(0, weight=1)
        self.trees = []
        self.tree = None

        for i, p in enumerate(self.provider):
            t = BoundTree(treeFrame, p,
                          selectmode=self._selectmode)
            t.grid(row=0, column=i, padx=5, pady=5, sticky='news')
            treeFrame.columnconfigure(i, weight=1)
            self.trees.append(t)

    def apply(self):
        self.values = []
        for t in self.trees:
            self.values.append(t.getSelectedObjects())