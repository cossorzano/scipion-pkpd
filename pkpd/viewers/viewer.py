# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@gmail.com)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

from os.path import exists
import numpy as np

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pkpd.objects import PKPDExperiment, PKPDAllometricScale
from pyworkflow.gui.text import openTextFileEditor
from pyworkflow.em.viewers.plotter import EmPlotter
from pkpd.objects import PKPDFitting, PKPDSignalAnalysis
from pkpd.pkpd_units import strUnit

from pkpd.protocols.protocol_batch_create_experiment import BatchProtCreateExperiment
from pkpd.protocols.protocol_pkpd_export_to_csv import ProtPKPDExportToCSV
from pkpd.protocols.protocol_pkpd_statistics_labels import ProtPKPDStatisticsLabel
from pkpd.protocols.protocol_pkpd_regression_labels import ProtPKPDRegressionLabel
from pkpd.protocols.protocol_pkpd_ode_bootstrap import ProtPKPDODEBootstrap
from pkpd.protocols.protocol_pkpd_filter_population import ProtPKPDFilterPopulation
from pkpd.protocols.protocol_pkpd_merge_populations import ProtPKPDMergePopulations

from pkpd.viewers.tk_experiment import ExperimentWindow
from pkpd.viewers.tk_populations import PopulationWindow

class PKPDViewer(Viewer):
    def _createExperiment(self):
        """ Create a new experiment after manipulation of
        the currently displayed experiment and register the action.
        """
        sampleKeys = self.windowDisplayed.samplesTree.selection()
        samples = ';'.join([self.windowDisplayed.experiment.samples[k].sampleName
                            for k in sampleKeys])

        prot = self.protocol
        project = prot.getProject()
        newProt = project.newProtocol(BatchProtCreateExperiment)
        newProt.inputExperiment.set(self.windowDisplayed.experiment)
        newProt.listOfSamples.set(samples)
        newProt.newTitle.set(self.windowDisplayed._titleVar.get())
        newProt.newComment.set(self.windowDisplayed._commentText.getText())
        project.launchProtocol(newProt)


class PKPDExperimentViewer(PKPDViewer):
    """ Visualization of a given PKPDExperiment
    """
    _targets = [PKPDExperiment]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        obj.load()
        self.windowDisplayed = self.tkWindow(ExperimentWindow,
                                           title='Experiment Viewer',
                                           experiment=obj,
                                           callback=self._createExperiment)
        self.windowDisplayed.show()


class PKPDFittingViewer(PKPDViewer):
    """ Visualization of a given PKPDFitting
    """
    _targets = [PKPDFitting]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        fitting = obj
        if fitting.isPopulation():
            fitting.sampleFittingClass="PKPDSampleFitBootstrap"
            fitting.load()
            self.windowDisplayed = self.tkWindow(PopulationWindow,
                                                  title='Population Viewer',
                                                  population=fitting)
            self.windowDisplayed.show()
        else:
            fitting.load()
            if hasattr(self.protocol,"outputExperiment"):
                experiment = self.protocol.outputExperiment
                experiment.load()
            else:
                experiment = fitting.loadExperiment()

            self.windowDisplayed = self.tkWindow(ExperimentWindow,
                                               title='Fitting Viewer',
                                               experiment=experiment,
                                               fitting=fitting,
                                               callback=self._createExperiment)
            self.windowDisplayed.show()


class PKPDCSVViewer(Viewer):
    """ Wrapper to visualize CSV files
    """
    _label = 'viewer csv'
    _targets = [ProtPKPDExportToCSV]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        fnCSV = self.protocol.getFilenameOut()

        if exists(fnCSV):
            openTextFileEditor(fnCSV)


class PKPDAnalysisViewer(Viewer):
    """ Wrapper to visualize analysis results
    """
    _label = 'viewer analysis'
    _targets = [PKPDSignalAnalysis]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        fnAnalysis = obj.fnAnalysis.get()

        if exists(fnAnalysis):
            openTextFileEditor(fnAnalysis)


class PKPDStatisticsLabelViewer(Viewer):
    """ Wrapper to visualize statistics
    """
    _label = 'viewer statistics'
    _targets = [ProtPKPDStatisticsLabel]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        fnStatistics = self.protocol._getPath("statistics.txt")

        if exists(fnStatistics):
            openTextFileEditor(fnStatistics)


class PKPDRegressionLabelsViewer(Viewer):
    """ Wrapper to visualize regression
    """
    _label = 'viewer regression'
    _targets = [ProtPKPDRegressionLabel]
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **kwargs):
        fnResults = self.protocol._getPath("results.txt")

        if exists(fnResults):
            X, Y, _, _ = self.protocol.getXYValues(False)

            minX = min(X)
            maxX = max(X)
            step = (maxX - minX) / 50
            xValues = np.arange(minX, maxX+step, step)
            yValues = self.protocol.evalFunction(xValues)

            plotter = EmPlotter(style='seaborn-whitegrid')
            varNameX = self.protocol.labelX.get()
            varNameY = self.protocol.labelY.get()
            ax = plotter.createSubPlot("Regression Plot",
                                       "%s [%s]"%(varNameX,strUnit(self.protocol.experiment.variables[varNameX].units.unit)),
                                       "%s [%s]"%(varNameY,strUnit(self.protocol.experiment.variables[varNameY].units.unit)))
            ax.plot(xValues, yValues)
            ax.plot(X, Y, 'o')

            return [plotter]
        else:
            return [self.errorMessage("Result file '%s' not produced yet. "
                                      % fnResults)]


class PKPDPopulationViewer(Viewer):
    """ Visualization of a given Population
    """
    _targets = [ProtPKPDODEBootstrap,
                ProtPKPDFilterPopulation,
                ProtPKPDMergePopulations]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        if hasattr(obj,"outputPopulation"):
            population = PKPDFitting("PKPDSampleFitBootstrap")
            population.load(obj.outputPopulation.fnFitting)

            self.populationWindow = self.tkWindow(PopulationWindow,
                                                  title='Population Viewer',
                                                  population=population)
            self.populationWindow.show()

class PKPDAllometricScalingViewer(Viewer):
    """ Visualization of an allometric scaling
    """
    _targets = [PKPDAllometricScale]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        model = PKPDAllometricScale()
        model.load(obj.fnScale.get())

        x = np.log10(np.asarray(model.X))
        xlabel = "%s [%s]" % (model.predictor, model.predictorUnits)
        for varName, varUnits in model.scaled_vars:
            plotter = EmPlotter(style='seaborn-whitegrid')
            y = np.log10(np.asarray(model.Y[varName]))
            ylabel = "%s [%s]" % (varName, varUnits)
            ax = plotter.createSubPlot(varName, xlabel, ylabel)
            ax.plot(x, y, '.', label='Species')
            ax.plot(x, np.log10(model.models[varName][0])+x*model.models[varName][1],'r', label='R2=%f'%model.qualifiers[varName][0])
            leg = ax.legend(loc='upper right')
            if leg:
                leg.draggable()
            plotter.show()

        for varName, varUnits in model.averaged_vars:
            plotter = EmPlotter(style='seaborn-whitegrid')
            y = np.asarray(model.Y[varName])
            ylabel = "%s [%s]" % (varName, varUnits)
            ax = plotter.createSubPlot("Scatter Plot", xlabel, ylabel)
            ax.plot(x, y, '.', label='Species')
            ax.plot(x, model.models[varName][0]*np.ones(x.shape),'r',label='Std=%f'%model.qualifiers[varName][0])
            leg = ax.legend(loc='upper right')
            if leg:
                leg.draggable()
            plotter.show()

