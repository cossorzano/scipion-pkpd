# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@gmail.com)
# *              Carlos Oscar Sorzano (info@kinestat.com)
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

import math
import numpy as np
from itertools import izip
from datetime import datetime
import Tkinter as tk

import pyworkflow.gui.dialog as dialog
import pyworkflow.gui as gui
from pyworkflow.gui.tree import TreeProvider, BoundTree
from pyworkflow.em.viewers.plotter import EmPlotter
from pkpd.pkpd_units import strUnit


class SamplesTreeProvider(TreeProvider):
    def __init__(self, experiment, fitting=None):
        self.experiment = experiment

    def getColumns(self):
        return [('Name', 60), ('Dose', 60)]

    def getObjects(self):
        sortedSamples = []
        for key in sorted(self.experiment.samples.keys()):
            sample = self.experiment.samples[key]
            sortedSamples.append(sample)
        return sortedSamples

    def getObjectInfo(self, obj):
        key = obj.sampleName
        values = [','.join(obj.doseList)]

        return {'key': key, 'text': key,
                'values': tuple(values)
                }


class MinMaxSlider(tk.Frame):
    """
    Create a personalized frame that contains:
        label, min entry, slider and max entry
    It also keeps a variable with the value
    """
    def __init__(self, master, label, from_=0, to=100, callback=None,
                 numberOfSteps=25):

        self.callback = callback
        self.numberOfSteps = numberOfSteps
        step = (to - from_) / numberOfSteps
        value = (from_ + to) * 0.5
        self.var = tk.DoubleVar()
        self.var.set(float(value))
        self.varMin = tk.DoubleVar()
        self.varMin.set(float(from_))
        self.varMax = tk.DoubleVar()
        self.varMax.set(float(to))

        tk.Frame.__init__(self, master)
        tk.Label(self, text=label).pack(side=tk.LEFT, padx=2, pady=2,
                                        anchor='s')

        def _entry(var):
            entry = tk.Entry(self, textvariable=var, bg='white', width=10)
            entry.pack(side=tk.LEFT, padx=2, pady=2, anchor='s')
            entry.bind('<Return>', self._onBoundChanged)

        _entry(self.varMin)

        self.slider = tk.Scale(self, from_=from_, to=to, variable=self.var,
                               bigincrement=step, resolution=step,
                               orient=tk.HORIZONTAL)
        self.slider.pack(side=tk.LEFT, padx=2)
        self.slider.bind('<ButtonRelease-1>', self._onButtonRelease)

        _entry(self.varMax)

    def getMinMax(self):
        return (self.varMin.get(), self.varMax.get())

    def getValue(self):
        return self.var.get()

    def _onBoundChanged(self, e=None):
        v = self.getValue()
        minV, maxV = self.getMinMax()
        step = (maxV - minV) / self.numberOfSteps
        self.slider.config(from_=minV, to=maxV,
                           bigincrement=step, resolution=step)
        if v > maxV or v < minV:
            self.var.set(0.5 * (minV + maxV))

    def _onButtonRelease(self, e=None):
        if self.callback:
            self.callback(e)


class PKPDResponsiveDialog(dialog.Dialog):
    def __init__(self, parent, title, **kwargs):
        """ From kwargs:
                message: message tooltip to show when browsing.
                selected: the item that should be selected.
                validateSelectionCallback: a callback function to validate selected items.
        """
        self.values = []
        self.plotter = None
        self.targetProtocol = kwargs['targetProtocol']
        self.experiment = self.targetProtocol.experiment
        self.varNameX = kwargs['varNameX']
        self.varNameY = kwargs['varNameY']
        self.provider = SamplesTreeProvider(self.experiment)
        self.model = self.loadModel()
        self.validateSelectionCallback = kwargs.get('validateSelectionCallback', None)
        self.setLabels()

        dialog.Dialog.__init__(self, parent, title,
                        buttons=[('Select', dialog.RESULT_YES),
                                 ('Cancel', dialog.RESULT_CANCEL)])

    def setLabels(self):
        pass

    def loadModel(self):
        model = self.targetProtocol.createModel()
        model.setExperiment(self.experiment)
        # if hasattr(self.protODE, "deltaT"):
        #     model.deltaT = self.protODE.deltaT.get()
        model.setXVar(self.varNameX)
        model.setYVar(self.varNameY)
        return model

    def body(self, bodyFrame):
        bodyFrame.config(bg='white')
        gui.configureWeigths(bodyFrame)
        self._createSamplesFrame(bodyFrame)
        self._createSlidersFrame(bodyFrame)
        self._createLogsFrame(bodyFrame)

    def _createSamplesFrame(self, content):
        frame = tk.Frame(content, bg='white')
        #frame = tk.LabelFrame(content, text='General')
        lfSamples = tk.LabelFrame(frame, text="Samples", bg='white')
        gui.configureWeigths(frame)
        lfSamples.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        self.samplesTree = self._addBoundTree(lfSamples,
                                                  self.provider, 10)
        self.samplesTree.itemClick = self._onSampleChanged

        frame.grid(row=0, column=0, sticky='news', padx=5, pady=5)

    def _createSlidersFrame(self, content):
        frame = tk.Frame(content, bg='white')
        lfBounds = tk.LabelFrame(frame, text="Parameter Bounds", bg='white')
        gui.configureWeigths(frame)

        i = 0
        self.sliders = {}
        paramUnits = self.targetProtocol.parameterUnits
        for paramName, bounds in self.targetProtocol.getParameterBounds().iteritems():
            bounds = bounds or (0, 1)
            slider = MinMaxSlider(lfBounds, "%s [%s]"%(paramName,strUnit(paramUnits[i])),
                                  bounds[0], bounds[1],
                                  callback=self._onVarChanged)
            slider.grid(row=i, column=0, padx=5, pady=5)
            self.sliders[paramName] = slider
            i += 1

        lfBounds.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        frame.grid(row=0, column=1, sticky='news', padx=5, pady=5)

    def _createLogsFrame(self, content):
        frame = tk.Frame(content)

        def addVar(text, col, varName):
            varFrame = tk.Frame(frame)
            varFrame.grid(row=0, column=col, sticky='new')
            label = tk.Label(varFrame, text=text)#, font=self.fontBold)
            label.grid(row=0, column=0, padx=5, pady=2, sticky='nw')
            combo = tk.Label(varFrame, text=varName, width=10)
            combo.grid(row=0, column=1, sticky='nw', padx=5, pady=5)
            radioVar = tk.IntVar()
            radio = tk.Checkbutton(varFrame, text='Log10', variable=radioVar)
            radio.grid(row=0, column=2, sticky='nw', padx=5, pady=5)
            return combo, radio, radioVar

        self.timeWidget = addVar('Time variable', 0, self.varNameX)
        self.timeWidget[2].trace('w', self._onLogChanged)
        self.measureWidget = addVar('Measure variable', 1, self.varNameY)
        measureVar = self.measureWidget[2]
        measureVar.set(True)
        measureVar.trace('w', self._onLogChanged)
        frame.grid(row=1, column=0, columnspan=2, sticky='news', padx=5, pady=5)

    def _addBoundTree(self, parent, provider, height):
        bt = BoundTree(parent, provider, height=height)
        bt.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(parent)
        return bt

    def apply(self):
        self.values = []

    def _onVarChanged(self, *args):
        sampleKeys = self.samplesTree.selection()

        if sampleKeys:
            self.computeFit()
            self.plotResults()
        else:
            dialog.showInfo("Warning","Please select some sample(s) to plot.",self)

    def computeFit(self):
        currentParams = []
        for paramName in self.targetProtocol.getParameterNames():
            currentParams.append(self.sliders[paramName].getValue())

        self.targetProtocol.setParameters(currentParams)
        self.ypValues = self.targetProtocol.forwardModel(currentParams, self.xpValues)

    def getBoundsList(self):
        boundList = []
        for paramName in self.targetProtocol.getParameterNames():
            boundList.append(self.sliders[paramName].getMinMax())
        return boundList

    def useTimeLog(self):
        return self.timeWidget[2].get()

    def useMeasureLog(self):
        return self.measureWidget[2].get()

    def getUnits(self, varName):
        return self.experiment.variables[varName].getUnitsString()

    def getLabel(self, varName, useLog):
        varLabel = '%s [%s]' % (varName, self.getUnits(varName))
        if useLog:
            varLabel = "log10(%s)" % varLabel
        return varLabel

    def getTimeLabel(self):
        return self.getLabel(self.varNameX, self.useTimeLog())

    def getMeasureLabel(self):
        return self.getLabel(self.varNameY, self.useMeasureLog())

    def computePlotValues(self, xValues, yValues):
        useMeasureLog = self.useMeasureLog()
        useTimeLog = self.useTimeLog()

        if not (useMeasureLog or useTimeLog):
            newXValues = xValues
            newYValues = yValues
        else:
            # If log will be used either for time or measure var
            # we need to filter elements larger than 0
            newXValues = []
            newYValues = []

            def _value(v, useLog):
                if useLog:
                    return math.log10(v) if v > 0 else None
                return v

            for x, y in izip(xValues, yValues):
                x = _value(x, useTimeLog)
                y = _value(y, useMeasureLog)

                if x is not None and y is not None:
                    newXValues.append(x)
                    newYValues.append(y)

        return newXValues, newYValues

    def _updateModel(self):
        """ This function should be called whenever the sample changes """
        pass

    def _onLogChanged(self, *args):
        # We will treat a log change as a sample change to plot
        self._onSampleChanged()

    def _onSampleChanged(self, e=None):
        sampleKeys = self.samplesTree.selection()

        if sampleKeys:
            # When the sample is changed we need to re-compute (with log or not)
            # the x, y values
            self.sample = self.experiment.samples[sampleKeys[0]]
            self.xValues, self.yValues = self.sample.getXYValues(self.varNameX,
                                                                 self.varNameY)
            self.newXValues, self.newYValues = self.computePlotValues(self.xValues[0],
                                                                      self.yValues[0])
            self._updateModel()
            self.computeFit()
            self.plotResults()
        else:
            dialog.showInfo("Warning","Please select some sample(s) to plot.",self)

    def plotResults(self):
        if self.plotter is None or self.plotter.isClosed():
            self.plotter = EmPlotter(style='seaborn-whitegrid')
            doShow = True
        else:
            doShow = False
            ax = self.plotter.getLastSubPlot()
            self.plotter.clear()

        ax = self.plotter.createSubPlot("Sample: %s" % self.sample.sampleName,
                                        self.getTimeLabel(),
                                        self.getMeasureLabel())
        self.newXPValues, self.newYPValues = self.computePlotValues(self.xpValues[0],
                                                                    self.ypValues[0])
        ax.plot(self.newXValues, self.newYValues, 'x', label="Observations")
        ax.plot(self.newXPValues, self.newYPValues, label="Fit")
        ax.legend()

        if doShow:
            self.plotter.show()
        else:
            self.plotter.draw()

    def destroy(self):
        """Destroy the window"""
        if not (self.plotter is None or self.plotter.isClosed()):
            self.plotter.close()
        dialog.Dialog.destroy(self)


class PKPDODEDialog(PKPDResponsiveDialog):
    def _updateModel(self):
        self.xpValues = [np.asarray([x for x in np.arange(0,np.max(self.xValues),4)])]
        self.targetProtocol.model.t0 = 0
        self.targetProtocol.model.tF = np.max(self.xValues)
        self.targetProtocol.drugSource.setDoses(self.sample.parsedDoseList,
                                                self.targetProtocol.model.t0,
                                                self.targetProtocol.model.tF)
        self.targetProtocol.configureSource(self.targetProtocol.drugSource)
        self.targetProtocol.model.drugSource = self.targetProtocol.drugSource
        # Necessary to count the number of source and PK parameters
        self.targetProtocol.getParameterNames()

class PKPDFitDialog(PKPDResponsiveDialog):
    def _updateModel(self):
        xLength=np.max(self.xValues)-np.min(self.xValues)
        self.xpValues = np.asarray([x for x in np.arange(np.min(self.xValues),np.max(self.xValues),xLength/25)])
        self.targetProtocol.getParameterNames()
