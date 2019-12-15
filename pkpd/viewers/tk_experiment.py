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
import Tkinter as tk
import ttk
from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.object as pwobj
from pkpd.utils import uniqueFloatValues
import pyworkflow.gui as gui
from pyworkflow.gui.widgets import Button, HotButton, ComboBox
from pyworkflow.gui.text import TaggedText
from pyworkflow.gui.tree import TreeProvider, BoundTree
from pyworkflow.em.viewers.plotter import EmPlotter


class VariablesTreeProvider(TreeProvider):
    def __init__(self, experiment):
        self.experiment = experiment

    def getColumns(self):
        return [('Name', 60), ('Unit', 60), ('Format', 60),
                ('Type', 60), ('Role', 60), ('Comment', 100)]

    def getObjects(self):
        sortedVars = []
        for key in sorted(self.experiment.variables.keys()):
            sortedVars.append(self.experiment.variables[key])
        return sortedVars

    def getObjectInfo(self, obj):
        key = obj.varName
        return {'key': key, 'text': key,
                'values': (obj.getUnitsString(),
                           obj.displayString,
                           obj.getTypeString(),
                           obj.getRoleString(),
                           obj.comment)}


class DosesTreeProvider(TreeProvider):
    def __init__(self, experiment):
        self.experiment = experiment

    def getColumns(self):
        return [('Name', 60), ('Via',40), ('Dose', 100), ('Amount', 80)]

    def getObjects(self):
        sortedDoses = []
        for key in sorted(self.experiment.doses.keys()):
            sortedDoses.append(self.experiment.doses[key])
        return sortedDoses

    def getObjectInfo(self, obj):
        key = obj.doseName
        return {'key': key, 'text': key,
                'values': (obj.via.viaName,
                           obj.getDoseString(),
                           obj.doseAmount+" "+obj.getDUnitsString()
                          )
                }


class GroupsTreeProvider(TreeProvider):
    def __init__(self, experiment):
        self.experiment = experiment

    def getColumns(self):
        return [('Name', 60), ('Samples', 60)]

    def getObjects(self):
        sortedGroups = []
        for key in sorted(self.experiment.groups.keys()):
            sortedGroups.append(self.experiment.groups[key])
        return sortedGroups

    def getObjectInfo(self, obj):
        key = obj.groupName
        return {'key': key, 'text': key,
                'values': (obj.getSamplesString()
                          )
                }

class ViasTreeProvider(TreeProvider):
    def __init__(self, experiment):
        self.experiment = experiment

    def getColumns(self):
        return [('Name', 60), ('Tlag', 60), ('Bioavailability',60)]

    def getObjects(self):
        sortedVias = []
        for key in sorted(self.experiment.vias.keys()):
            sortedVias.append(self.experiment.vias[key])
        return sortedVias

    def getObjectInfo(self, obj):
        key = obj.viaName
        if 'tlag' in obj.paramsToOptimize:
            tlagStr="To be optimized"
        else:
            tlagStr=str(obj.tlag)+" "+obj.tunits._toString()
        if 'bioavailability' in obj.paramsToOptimize:
            bioavailabilityStr="To be optimized"
        else:
            bioavailabilityStr=str(obj.bioavailability)
        return {'key': key, 'text': key,
                'values': (tlagStr,
                           bioavailabilityStr
                          )
                }


class SamplesTreeProvider(TreeProvider):
    FIT_COLUMNS = ['R2', 'R2adj', 'AIC', 'AICc', 'BIC']

    def __init__(self, experiment, fitting=None):
        self.experiment = experiment
        self.fitting = fitting

        numberOfSamples = len(experiment.samples)
        if numberOfSamples > 0:
            sample = self.experiment.samples.values()[0]
            if sample.descriptors:
                self.columns = [(key, 60)
                                for key in sorted(sample.descriptors.keys())]
            else:
                self.columns = []
            if self.fitting:
                self.columns.extend([(key, 60) for key in self.FIT_COLUMNS])
        else:
            self.columns = []

    def getColumns(self):
        return [('Name', 60), ('Dose', 60),('Groups',60)] + self.columns

    def getObjects(self):
        sortedSamples = []
        for key in sorted(self.experiment.samples.keys()):
            sample = self.experiment.samples[key]
            if self.fitting:
                sample.fit = self.fitting.getSampleFit(key)
            sortedSamples.append(sample)
        return sortedSamples

    def getObjectInfo(self, obj):
        key = obj.sampleName
        values = [','.join(obj.doseList)]
        values +=[','.join(obj.groupList)]
        if obj.descriptors:
            values += [obj.descriptors[k] for k, _ in self.columns if k in obj.descriptors]

        if self.fitting:
            if obj.fit!=None:
                values.extend([getattr(obj.fit, attr) for attr in self.FIT_COLUMNS])

        return {'key': key, 'text': key,
                'values': tuple(values)
                }


class MeasurementTreeProvider(TreeProvider):
    def __init__(self, sample):
        self.sample = sample
        self.columns = [(key, 60) for key in self.sample.measurementPattern]

    def getColumns(self):
        return [('Sample #',60)]+self.columns

    def getObjects(self):
        return self.sample.getSampleMeasurements()

    def getObjectInfo(self, obj):
        key = obj.n
        return {'key': key, 'text': key,
                'values': tuple(obj.getValues())
                }


class FitValuesTreeProvider(TreeProvider):
    def __init__(self, fitting, sampleFit):
        self.sampleFit = sampleFit
        if type(fitting.predicted)==list:
            self.columns = [('n',30), (fitting.predictor.getLabel(), 60)]+[(y.getLabel(),60) for y in fitting.predicted]+\
                           [("Fitted Value", 60)]
        else:
            self.columns = [('n',30), (fitting.predictor.getLabel(), 60),
                            (fitting.predicted.getLabel(), 60),
                            ("Fitted Value", 60)]

    def getColumns(self):
        return self.columns

    def getObjects(self):
        fit = self.sampleFit

        def _object(n, x, y, yp):
            obj = pwobj.Object()
            obj.n = n
            obj.x = x
            obj.y = y
            obj.yp = yp
            return obj

        if type(fit.y[0])==list:
            retval = []
            Ncols=len(fit.y)
            n=0
            i=0
            for xi, yi, ypi in izip(fit.x,fit.y,fit.yp):
                for x,y,yp in izip(xi,yi,ypi):
                    yList=[]
                    for ii in range(Ncols):
                        if ii==i:
                            yList.append(y)
                        else:
                            yList.append("NA")
                    retval.append(_object(n, x, yList, yp))
                    n+=1
                i+=1
            return retval
        else:
            nrange=[x for x in range(len(fit.x))]
            return [_object(n, x, y, yp) for n, x, y, yp in izip(nrange, fit.x, fit.y, fit.yp)]

    def getObjectInfo(self, obj):
        key = obj.n
        if type(obj.y)==list:
            return {'key': key, 'text': key,
                    'values': tuple([obj.x]+obj.y+[obj.yp])
                    }
        else:
            return {'key': key, 'text': key,
                    'values': (obj.x, obj.y, obj.yp)
                    }


class ExperimentWindow(gui.Window):
    """ This class creates a Window that will display some Point's
    contained in a Data object.
    It will allow to launch 1D, 2D and 3D plots by selecting any
    combination of the x1, x2...xn from the Point dimension.
    Points can be selected by either Click and Drag in the Scatter plot or..
    by creating an Expression.
    Finally, there is a button 'Create Cluster' that will call a callback
    fuction to take care of it.
    """
    def __init__(self, **kwargs):
        gui.Window.__init__(self,  minsize=(420, 200), **kwargs)
        self.experiment = kwargs.get('experiment')
        self.fitting = kwargs.get('fitting', None)
        self.callback = kwargs.get('callback', None)
        content = tk.Frame(self.root)
        self._createContent(content)
        content.grid(row=0, column=0, sticky='news')
        content.columnconfigure(0, weight=1)
        content.rowconfigure(0, weight=1)
        self.plotter = None

    def _createContent(self, content):
        # Create and fill the frame containing the Experiment
        # info, variables and doses
        p = tk.PanedWindow(content, orient=tk.VERTICAL, bg='white')

        p.grid(row=0, column=0, sticky='news', padx=5, pady=5)

        self._createTopFrame(p)
        # Create the middle frame containing the Samples Box
        self._createSamplesFrame(p)
        # Create the last frame with the buttons
        self._createButtonsFrame(content)
        #self._updateSelectionLabel()

    def _createTopFrame(self, content):
        frame = tk.Frame(content)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        tab = ttk.Notebook(frame)
        tab.grid(row=0, column=0, sticky='news', padx=5, pady=5)

        def addTab(label):
            lf = tk.Frame(tab)
            tab.add(lf, text=label)
            return lf

        lfGeneral = addTab('General')
        self._addLabel(lfGeneral, 'Title', 0, 0)
        self._titleVar = tk.StringVar()
        self._titleVar.set(self.experiment.general['title'])
        titleEntry = tk.Entry(lfGeneral, width=26, textvariable=self._titleVar,
                              bg='white')
        titleEntry.grid(row=0, column=1, sticky='nw', padx=5, pady=(5, 0))
        self._addLabel(lfGeneral, 'Comment', 1, 0)
        commentText = gui.text.Text(lfGeneral, width=30, height=3, bg='white')
        commentText.setText(self.experiment.general['comment'])
        commentText.grid(row=1, column=1, sticky='nw', padx=5, pady=(5, 0))
        self._commentText = commentText

        lfVars = addTab('Variables')
        self.varsTree = self._addBoundTree(lfVars, VariablesTreeProvider, 5)

        lfVias = addTab('Vias')
        self.viasTree = self._addBoundTree(lfVias, ViasTreeProvider, 5)

        lfDoses = addTab('Doses')
        self.dosesTree = self._addBoundTree(lfDoses, DosesTreeProvider, 5)

        lfGroups = addTab('Groups')
        self.groupsTree = self._addBoundTree(lfGroups, GroupsTreeProvider, 5)

        # Fitting tab
        if self.fitting:
            tabFitting = addTab('Fitting')
            t = TaggedText(tabFitting, width=80, height=10, bg='white')
            t.grid(row=0, column=0, sticky='nw', padx=10, pady=10)
            t.setText('Select one of the samples to see more information.')
            t.setReadOnly(True)
            self.fittingText = t
            # Buttons for fitting
            buttonsFrame = self._createFittingButtonsFrame(tabFitting)
            buttonsFrame.grid(row=0, column=1, sticky='news', padx=5, pady=5)

        #frame.grid(row=0, column=0, sticky='news', padx=5, pady=(10, 5))
        content.add(frame, sticky='news', padx=5, pady=5)

    def _createSamplesFrame(self, content):
        frame = tk.Frame(content)
        #frame = tk.LabelFrame(content, text='General')
        lfSamples = tk.LabelFrame(frame, text='Samples')
        gui.configureWeigths(frame)
        lfSamples.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        self.samplesTree = self._addBoundTree(lfSamples, SamplesTreeProvider,
                                              10, fitting=self.fitting)
        self.samplesTree.itemDoubleClick = self._onSampleDoubleClick
        self.samplesTree.itemClick = self._onSampleClick

        plotFrame = tk.Frame(lfSamples)
        plotFrame.grid(row=1, column=0, sticky='ws', padx=5, pady=5)

        # Add a combobox with the variable for time
        timeVars = [v.varName for v in self.experiment.variables.values()
                    if v.role == v.ROLE_TIME]
        timeVars += [v.varName for v in self.experiment.variables.values()
                    if v.role == v.ROLE_MEASUREMENT]

        measureVars = [v.varName for v in self.experiment.variables.values()
                    if v.role == v.ROLE_MEASUREMENT]

        def addVar(text, col, choices):
            varFrame = tk.Frame(plotFrame)
            varFrame.grid(row=0, column=col, sticky='new')
            label = tk.Label(varFrame, text=text, font=self.fontBold)
            label.grid(row=0, column=0, padx=5, pady=2, sticky='nw')
            combo = ComboBox(varFrame, choices, width=10)
            combo.grid(row=0, column=1, sticky='nw', padx=5, pady=5)
            radioVar = tk.IntVar()
            radio = tk.Checkbutton(varFrame, text='Log10', variable=radioVar)
            radio.grid(row=0, column=2, sticky='nw', padx=5, pady=5)
            return combo, radio, radioVar

        self.timeWidget = addVar('Time variable', 0, timeVars)
        self.measureWidget = addVar('Measure variable', 1, list(set(measureVars+timeVars)))
        self.measureWidget[2].set(True)

        self.plotButton = Button(plotFrame, '   Plot   ', font=self.fontBold,
                                 command=self._onPlotClick,
                                 tooltip='Select one or more samples to plot '
                                         'their measures of the selected '
                                         'variables (optionally in log).')

        self.plotButton.grid(row=0, column=2, sticky='ne', padx=5)

        self.useCurrentPlotVar = tk.IntVar()
        self.useCurrentPlot = tk.Checkbutton(plotFrame, text='Use current plot', variable=self.useCurrentPlotVar)
        self.useCurrentPlot.grid(row=1, column=0, sticky='nw', padx=5, pady=5)

        self.plotSummaryButton = Button(plotFrame, ' Summary Plot ',
                                        font=self.fontBold,
                                        command=self._onPlotSummaryClick,
                                        tooltip='Select several samples to plot'
                                                ' their statistics'
                                                ' (min, max, avg and std).')

        self.plotSummaryButton.grid(row=1, column=2, sticky='ne',
                                    padx=5, pady=5)

        #frame.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        content.add(frame, sticky='news', padx=5, pady=5)

    def _createButtonsFrame(self, content):
        frame = tk.Frame(content)
        gui.configureWeigths(frame)
        buttonsFrame = tk.Frame(frame)
        buttonsFrame.grid(row=0, column=0, sticky='ne')
        closeButton = Button(buttonsFrame, 'Close', command=self.close,
                             imagePath='fa-times.png')
        closeButton.grid(row=0, column=0, sticky='ne', padx=5)
        self.newButton = HotButton(buttonsFrame, '   New Experiment   ',
                                   command=self._onCreateClick,
                                   tooltip='Create a new experiment with the '
                                            'selected samples. You can also edit'
                                            'title and comment.')

        self.newButton.grid(row=0, column=1, sticky='ne', padx=5)

        frame.grid(row=1, column=0, sticky='news', padx=5, pady=5)

    def _createFittingButtonsFrame(self, content):
        buttonsFrame = tk.Frame(content)
        buttonsFrame.grid(row=0, column=0, sticky='ne')

        plotValuesButton = Button(buttonsFrame, 'Plot Fit',
                                  command=self._onPlotFitClick)
        plotValuesButton.grid(row=0, column=0, sticky='sew', padx=5, pady=5)

        openValuesButton = Button(buttonsFrame, 'Open Fit',
                                  command=self._onOpenFitClick)
        openValuesButton.grid(row=1, column=0, sticky='sew', padx=5, pady=5)

        return buttonsFrame


    def _addLabel(self, parent, text, r, c):
        label = tk.Label(parent, text=text, font=self.fontBold)
        label.grid(row=r, column=c, padx=5, pady=5, sticky='ne')
        return label

    def _addBoundTree(self, parent, ProviderClass, height, **kwargs):
        bt = BoundTree(parent, ProviderClass(self.experiment, **kwargs),
                       height=height)
        bt.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(parent)
        return bt

    def getUnits(self, varName):
        return self.experiment.variables[varName].getUnitsString()

    def getLabel(self, varName, useLog):
        varLabel = '%s [%s]' % (varName, self.getUnits(varName))
        if useLog:
            varLabel = "log10(%s)" % varLabel
        return varLabel

    def getTimeVarName(self):
        return self.timeWidget[0].getText()

    def useTimeLog(self):
        return self.timeWidget[2].get()

    def getTimeLabel(self):
        return self.getLabel(self.getTimeVarName(), self.useTimeLog())

    def getMeasureVarName(self):
        return self.measureWidget[0].getText()

    def useMeasureLog(self):
        return self.measureWidget[2].get()

    def getMeasureLabel(self):
        return self.getLabel(self.getMeasureVarName(), self.useMeasureLog())

    def getPlotValues(self, sample):
        xValues, yValues = sample.getXYValues(self.getTimeVarName(),
                                  self.getMeasureVarName())
        xValues=xValues[0] # From [array(...)] to array(...)
        yValues=yValues[0]

        useMeasureLog = self.useMeasureLog()
        useTimeLog = self.useTimeLog()

        if not (useMeasureLog or useTimeLog):
            return xValues, yValues

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

    def reuseFigure(self):
        return None if not self.useCurrentPlotVar.get() else 'active'

    def _onPlotClick(self, e=None):
        sampleKeys = self.samplesTree.selection()

        if sampleKeys:
            if self.plotter is None or self.plotter.isClosed():
                self.plotter = EmPlotter(style='seaborn-whitegrid', figure=self.reuseFigure())
                doShow = True
                ax = self.plotter.createSubPlot("Plot", self.getTimeLabel(),
                                                self.getMeasureLabel())
                self.plotDict = {}
            else:
                doShow = False
                ax = self.plotter.getLastSubPlot()

            samples = [self.experiment.samples[k] for k in sampleKeys]
            for s in samples:
                if not s.sampleName in self.plotDict:
                    x, y = self.getPlotValues(s)
                    ax.plot(x, y, label=s.sampleName)
                    self.plotDict[s.sampleName] = True
            leg = ax.legend()
            if leg:
                leg.draggable()

            if doShow:
                self.plotter.show()
            else:
                self.plotter.draw()
        else:
            self.showInfo("Please select some sample(s) to plot.")

    def _onPlotSummaryClick(self, e=None):
        sampleKeys = self.samplesTree.selection()
        n = len(sampleKeys)

        if n == 1:
            self.showInfo("Please select several samples to plot.")
        else:
            if n > 1:
                samples = [self.experiment.samples[k] for k in sampleKeys]
            else:
                samples = self.experiment.samples.values()

            xmin=1e38
            xmax=-1e38
            for s in samples:
                xValues, _ = self.getPlotValues(s)
                xmin = min(xmin, min(xValues))
                xmax = max(xmax, max(xValues))

            dataDict = {} # key will be time values
            xrange=np.arange(xmin,xmax,(xmax-xmin)/300.0)
            for s in samples:
                xValues, yValues = self.getPlotValues(s)
                xValuesUnique,yValuesUnique = uniqueFloatValues(xValues,yValues)
                B = InterpolatedUnivariateSpline(xValuesUnique,yValuesUnique, k=1)
                yrange = B(xrange)
                for x, y in izip(xrange, yrange):
                    if x in dataDict:
                        dataDict[x].append(y)
                    else:
                        dataDict[x] = [y]

            sortedTime = sorted(dataDict.keys())
            # We will store five values (min, 25%, 50%, 75%, max)
            # for each of the time entries computed
            percentileList = [0, 25, 50, 75, 100]
            Y = np.zeros((len(sortedTime), 5))
            for i, t in enumerate(sortedTime):
                Y[i,:] = np.percentile(dataDict[t], percentileList)

            plotter = EmPlotter(style='seaborn-whitegrid', figure=self.reuseFigure())
            # *** Pending
            ax = plotter.createSubPlot("Summary Plot", self.getTimeLabel(),
                                       self.getMeasureLabel())
            ax.plot(sortedTime, Y[:, 0], 'r--', label="Minimum")
            ax.plot(sortedTime, Y[:, 1], 'b--', label="25%")
            ax.plot(sortedTime, Y[:, 2], 'g', label="50% (Median)")
            ax.plot(sortedTime, Y[:, 3], 'b--', label="75%")
            ax.plot(sortedTime, Y[:, 4], 'r--', label="Maximum")
            ax.grid(True)

            ax.legend()
            plotter.show()

    def _onCreateClick(self, e=None):
        sampleKeys = self.samplesTree.selection()
        if sampleKeys and self.callback:
            self.callback()
        else:
            self.showInfo("Please select some sample(s) to create a new experiment.")

    def _onSampleClick(self, obj):
        # In fitting mode we need to display some information for
        # the given sample fit
        if self.fitting:
            selection = self.samplesTree.selection()
            if selection:
                sampleKey = selection[0]
                sampleFit = self.fitting.getSampleFit(sampleKey)
                textMsg = sampleFit.getBasicInfo()
            else:
                textMsg = 'Select one of the samples to see more information.'

            self.fittingText.setReadOnly(False)
            self.fittingText.setText(textMsg)
            self.fittingText.setReadOnly(True)

        if not (self.plotter is None or self.plotter.isClosed()):
            self._onPlotClick()

    def _onSampleDoubleClick(self, obj):
        sampleKeys = self.samplesTree.selection()
        if sampleKeys:
            MeasureWindow(masterWindow=self,
                          sample=self.experiment.samples[sampleKeys[0]]).show()

    def _onPlotFitClick(self, e=None):
        sampleKeys = self.samplesTree.selection()

        if sampleKeys:
            def _value(v, useLog):
                if useLog:
                    return math.log10(v) if v > 0 else None
                return v

            def _values(values, useLog=self.useMeasureLog()):
                return [_value(float(x), useLog) for x in values]

            def _plot(varLabelX, varLabelY, x, y, yp):
                plotter = EmPlotter(style='seaborn-whitegrid', figure=self.reuseFigure())
                ax = plotter.createSubPlot("Plot", varLabelX, varLabelY)
                xValues = _values(x, useLog=self.useTimeLog())
                ax.plot(xValues, _values(y), 'x', label="Observations")
                idx = np.argsort(xValues)
                ypValues = _values(yp)
                ax.plot(np.asarray(xValues)[idx], np.asarray(ypValues)[idx], 'g', label="Fit")
                leg = ax.legend()
                if leg:
                    leg.draggable()

                plotter.show()

            # Get first selected element
            fit = self.fitting.getSampleFit(sampleKeys[0])
            if fit==None:
                return
            varLabelX = self.getLabel(self.fitting.predictor.varName,
                                      self.useTimeLog())
            if type(self.fitting.predicted)==list:
                i=0
                for v in self.fitting.predicted:
                    varLabelY = self.getLabel(v.varName, self.useMeasureLog())
                    _plot(varLabelX, varLabelY, fit.x[i], fit.y[i], fit.yp[i])
                    i+=1
            else:
                varLabelY = self.getLabel(self.fitting.predicted.varName,
                                          self.useMeasureLog())
                _plot(varLabelX, varLabelY, fit.x[0], fit.y[0], fit.yp[0])

        else:
            self.showInfo("Please select some sample(s) to plot.")

    def _onOpenFitClick(self, e=None):
        sampleKeys = self.samplesTree.selection()

        if sampleKeys:
            # Get first selected element
            fit = self.fitting.getSampleFit(sampleKeys[0])
            FitWindow(masterWindow=self,
                      fitting=self.fitting, sampleFit=fit).show()
        else:
            self.showInfo("Please select some sample(s) to plot.")

    def _onClosing(self):
        if self.plotter:
            self.plotter.close()
        gui.Window._onClosing(self)


class MeasureWindow(gui.Window):
    def __init__(self, **kwargs):
        gui.Window.__init__(self,  minsize=(420, 200), **kwargs)
        self.title = kwargs.get('title', 'Measurements')
        self.sample = kwargs.get('sample', None)
        self.provider = kwargs.get('provider', None)
        if not self.provider:
            self.provider = MeasurementTreeProvider(self.sample)
        content = tk.Frame(self.root)
        self._createContent(content)
        content.grid(row=0, column=0, sticky='news')
        content.columnconfigure(0, weight=1)

    def _createContent(self, content):
        self._createMeasurementFrame(content)
        self._createButtonsFrame(content)

    def _createMeasurementFrame(self, content):
        frame = tk.Frame(content)
        #frame = tk.LabelFrame(content, text='General')
        lfSamples = tk.LabelFrame(frame, text=self.title)
        gui.configureWeigths(frame)
        lfSamples.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        self.measurementTree = self._addBoundTree(lfSamples,
                                                  self.provider, 10)

        frame.grid(row=0, column=0, sticky='news', padx=5, pady=5)

    def _createButtonsFrame(self, content):
        frame = tk.Frame(content)
        gui.configureWeigths(frame)
        buttonsFrame = tk.Frame(frame)
        buttonsFrame.grid(row=0, column=0, sticky='ne')
        closeButton = Button(buttonsFrame, 'Close', command=self.close,
                             imagePath='fa-times.png')

        frame.grid(row=2, column=0, sticky='news', padx=5, pady=5)

    def _addBoundTree(self, parent, provider, height):
        bt = BoundTree(parent, provider, height=height)
        bt.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(parent)
        return bt


class FitWindow(MeasureWindow):
    def __init__(self, **kwargs):
        self.fitting = kwargs['fitting']
        self.sampleFit = kwargs['sampleFit']

        MeasureWindow.__init__(self, title="Fitted Values",
                               provider=FitValuesTreeProvider(self.fitting,
                                                              self.sampleFit),
                               **kwargs)
