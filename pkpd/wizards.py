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

import pyworkflow.object as pwobj
from pyworkflow.wizard import Wizard
import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.tree import TreeProvider

from pkpd.protocols.protocol_pkpd_regression_labels import ProtPKPDRegressionLabel
from pkpd.protocols.protocol_pkpd_drop_measurements import ProtPKPDDropMeasurements
from pkpd.protocols.protocol_pkpd_change_units import ProtPKPDChangeUnits
from pkpd.protocols.protocol_pkpd_scale_to_common_dose import ProtPKPDScaleToCommonDose
from pkpd.protocols.protocol_pkpd_stats_oneExperiment_twoSubgroups_mean import ProtPKPDStatsExp1Subgroups2Mean
from pkpd.protocols.protocol_pkpd_exponential_fit import ProtPKPDExponentialFit
from pkpd.protocols.protocol_pkpd_elimination_rate import ProtPKPDEliminationRate
from pkpd.protocols.protocol_pkpd_monocompartment import ProtPKPDMonoCompartment
from pkpd.protocols.protocol_pkpd_two_compartments import ProtPKPDTwoCompartments
from pkpd.protocols.protocol_pkpd_simulate_generic_pd import ProtPKPDSimulateGenericPD
from pkpd.protocols.protocol_pkpd_stats_twoExperiments_twoSubgroups_mean import ProtPKPDStatsExp2Subgroups2Mean
from pkpd.protocols.protocol_pkpd_import_from_csv import ProtPKPDImportFromText, getSampleNamesFromCSVfile, getVarNamesFromCSVfile
from pkpd.protocols.protocol_pkpd_bootstrap_simulate import ProtPKPDODESimulate
from pkpd.protocols.protocol_pkpd_monocompartment_urine import ProtPKPDMonoCompartmentUrine
from pkpd.protocols.protocol_pkpd_pdgeneric_fit import ProtPKPDGenericFit
from pkpd.protocols.protocol_pkpd_change_via import ProtPKPDChangeVia
from pkpd.protocols.protocol_pkpd_twocompartments_urine import ProtPKPDTwoCompartmentsUrine
from pkpd.protocols.protocol_pkpd_monocompartment_clint import ProtPKPDMonoCompartmentClint
from pkpd.protocols.protocol_pkpd_two_compartments_clint import ProtPKPDTwoCompartmentsClint
from pkpd.protocols.protocol_pkpd_two_compartments_metabolite import ProtPKPDTwoCompartmentsClintMetabolite
from pkpd.protocols.protocol_pkpd_two_compartments_autoinduction import ProtPKPDTwoCompartmentsAutoinduction
from pkpd.protocols.protocol_pkpd_twocompartments_both import ProtPKPDTwoCompartmentsBoth
from pkpd.protocols.protocol_pkpd_monocompartment_linkpd import ProtPKPDMonoCompartmentLinkPD
from pkpd.protocols.protocol_pkpd_monocompartment_pd import ProtPKPDMonoCompartmentPD
from pkpd.protocols.protocol_pkpd_twocompartments_both_pd import ProtPKPDTwoCompartmentsBothPD

from pkpd.tk_ode import PKPDODEDialog, PKPDFitDialog


class FilterVariablesTreeProvider(TreeProvider):
    """ Simplified view of VariablesTreeProvider with less columns.
    Additionally, we can filter by a given function. """

    def __init__(self, experiment, filter=None):
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
    _targets = [(ProtPKPDRegressionLabel, ['labelX']),
                (ProtPKPDRegressionLabel, ['labelY']),
                (ProtPKPDChangeUnits, ['labelToChange']),
                (ProtPKPDStatsExp1Subgroups2Mean, ['labelToCompare']),
                (ProtPKPDExponentialFit, ['predictor']),
                (ProtPKPDExponentialFit, ['predicted']),
                (ProtPKPDEliminationRate, ['predictor']),
                (ProtPKPDEliminationRate, ['predicted']),
                (ProtPKPDMonoCompartment, ['predictor']),
                (ProtPKPDMonoCompartment, ['predicted']),
                (ProtPKPDMonoCompartmentUrine, ['predictor']),
                (ProtPKPDMonoCompartmentUrine, ['predicted']),
                (ProtPKPDMonoCompartmentUrine, ['Au']),
                (ProtPKPDMonoCompartmentLinkPD, ['predictor']),
                (ProtPKPDMonoCompartmentLinkPD, ['predicted']),
                (ProtPKPDMonoCompartmentLinkPD, ['E']),
                (ProtPKPDMonoCompartmentPD, ['predictor']),
                (ProtPKPDMonoCompartmentPD, ['predicted']),
                (ProtPKPDMonoCompartmentPD, ['E']),
                (ProtPKPDTwoCompartments, ['predictor']),
                (ProtPKPDTwoCompartments, ['predicted']),
                (ProtPKPDTwoCompartmentsClint, ['predictor']),
                (ProtPKPDTwoCompartmentsClint, ['predicted']),
                (ProtPKPDTwoCompartmentsClintMetabolite, ['predictor']),
                (ProtPKPDTwoCompartmentsClintMetabolite, ['predicted']),
                (ProtPKPDTwoCompartmentsClintMetabolite, ['metabolite']),
                (ProtPKPDTwoCompartmentsAutoinduction, ['predictor']),
                (ProtPKPDTwoCompartmentsAutoinduction, ['predicted']),
                (ProtPKPDTwoCompartmentsBoth, ['predictor']),
                (ProtPKPDTwoCompartmentsBoth, ['predicted']),
                (ProtPKPDTwoCompartmentsBoth, ['Cperipheral']),
                (ProtPKPDTwoCompartmentsBothPD, ['predictor']),
                (ProtPKPDTwoCompartmentsBothPD, ['predicted']),
                (ProtPKPDTwoCompartmentsBothPD, ['Cperipheral']),
                (ProtPKPDTwoCompartmentsBothPD, ['E']),
                (ProtPKPDMonoCompartmentUrine, ['predicted']),
                (ProtPKPDSimulateGenericPD, ['predictor']),
                (ProtPKPDStatsExp2Subgroups2Mean, ['label1',
                                                   'inputExperiment1']),
                (ProtPKPDStatsExp2Subgroups2Mean, ['label2',
                                                   'inputExperiment2'])
                ]

    def show(self, form, *params):
        protocol = form.protocol
        label = params[0]

        fullParams = self._getAllParams(protocol)
        experiment = None
        if len(fullParams)>1:
            experiment = protocol.getAttributeValue(fullParams[1], None)
        else:
            experiment = protocol.getAttributeValue('inputExperiment', None)

        if experiment is None:
            form.showError("Select input experiment first.")
        else:
            experiment.load()
            filterFunc = getattr(protocol, 'filterVarForWizard', None)
            tp = FilterVariablesTreeProvider(experiment, filter=filterFunc)
            dlg = dialog.ListDialog(form.root, self.getTitle(), tp,
                                    selectmode=self.getSelectMode())
            if dlg.resultYes():
                self.setFormValues(form, label, dlg.values)

    def getTitle(self):
        return "Choose variable"

    def setFormValues(self, form, label, values):
        form.setVar(label, values[0].varName)

    def getSelectMode(self):
        return "browse"


class DoseTreeProvider(TreeProvider):
    def __init__(self, experiment):
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
    _targets = [(ProtPKPDChangeVia, ['doseName','doseVia']),
                ]

    def show(self, form, *params):
        protocol = form.protocol
        label = params[0]

        experiment = protocol.getAttributeValue('inputExperiment', None)

        if experiment is None:
            form.showError("Select input experiment first.")
        else:
            experiment.load()
            dlg = dialog.ListDialog(form.root, "Choose dose name",
                                    DoseTreeProvider(experiment),
                                    selectmode=self.getSelectMode())
            if dlg.resultYes():
                self.setFormValues(form, label, dlg.values)

    def setFormValues(self, form, label, values):
        form.setVar('doseName', values[0].doseName)
        form.setVar('doseVia', values[0].via)

    def getSelectMode(self):
        return "browse"


class PKPDChooseSeveralVariableWizard(PKPDChooseVariableWizard):
    _targets = [(ProtPKPDDropMeasurements, ['varsToDrop']),
                (ProtPKPDScaleToCommonDose,['measurementsToChange'])]

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
    _targets = [(ProtPKPDImportFromText, ['variables'])
                ]

    def show(self, form, *params):
        label = params[0]
        protocol = form.protocol
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
                if strToAdd!="":
                    currentValue = protocol.getAttributeValue(label, "")
                    form.setVar(label, currentValue+strToAdd)


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
        label = params[0]
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

            dlg = dialog.MultiListDialog(form.root, "Test",
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
                (ProtPKPDMonoCompartmentUrine, ['bounds']),
                (ProtPKPDMonoCompartmentClint, ['bounds']),
                (ProtPKPDMonoCompartmentLinkPD, ['bounds']),
                (ProtPKPDMonoCompartmentPD, ['bounds']),
                (ProtPKPDTwoCompartments, ['bounds']),
                (ProtPKPDTwoCompartmentsUrine, ['bounds']),
                (ProtPKPDTwoCompartmentsBoth, ['bounds']),
                (ProtPKPDTwoCompartmentsBothPD, ['bounds']),
                (ProtPKPDTwoCompartmentsClint, ['bounds']),
                (ProtPKPDTwoCompartmentsClintMetabolite, ['bounds']),
                (ProtPKPDTwoCompartmentsAutoinduction, ['bounds']),
                (ProtPKPDGenericFit, ['bounds']),
                ]

    _nonODE = [ProtPKPDGenericFit]

    def show(self, form, *params):
        label = params[0]
        protocol = form.protocol
        experiment = protocol.getAttributeValue('inputExperiment', None)

        if experiment is None:
            form.showError("Select the input experiment first.")
        else:
            # try:
                protocol.setInputExperiment() # this load the experiment
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
                    if i==0:
                        if not type(protocol) in PKPDODEWizard._nonODE:
                            protocol.model.drugSource.setDoses(sample.parsedDoseList,0,1) # Needed to correctly identify the parameters of the source
                        protocol.model.setSample(sample)
                        protocol.calculateParameterUnits(sample)
                    i+=1

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
                        if i>1:
                            boundStr+="; "
                        boundStr += str(bound)
                        i += 1
                    if boundStr!="":
                        form.setVar(label, boundStr)
            # except Exception as e:
                # pass
                # form.showError("Error: %s" % str(e))
