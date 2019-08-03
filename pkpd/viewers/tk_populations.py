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

import numpy as np
from collections import OrderedDict
import Tkinter as tk

import pyworkflow.object as pwobj
import pyworkflow.gui as gui
from pyworkflow.gui.widgets import HotButton
from pyworkflow.gui.tree import TreeProvider, BoundTree
from pyworkflow.gui.text import TaggedText
from pyworkflow.em.viewers.plotter import EmPlotter

from pkpd.pkpd_units import PKPDUnit


class PopulationVar(pwobj.Object):
    def __init__(self, **kwargs):
        self.data = kwargs
        for k, v in kwargs.iteritems():
            setattr(self, k, v)


class PopulationVariablesTreeProvider(TreeProvider):
    def __init__(self, variables):
        #self.population = population
        self._variables = variables

    def getColumns(self):
        w = 60
        return [('Name', w), ('Unit', w),
                ('Mean', w), ('Std', w),
                ('Minimum', w), ('2.5%', w), ('25%', w), ('50%', w),
                ('75%', w), ('97.5%', w), ('Maximum', w)]

    def getObjects(self):
        return self._variables

    def getObjectInfo(self, obj):
        key = obj.varName
        p = obj.percentiles
        return {'key': key, 'text': key,
                'values': (obj.unitStr,
                           obj.mu, obj.sigma,
                           p[0], p[1], p[2], p[3], p[4], p[5], p[6])}


class PopulationWindow(gui.Window):
    """
    """
    def __init__(self, **kwargs):
        gui.Window.__init__(self,  minsize=(620, 200), **kwargs)
        self.population = kwargs.get('population')
        self.callback = kwargs.get('callback', None)
        self.plotter = None
        self._variables = OrderedDict()
        self._variablesDict = {}
        self._loadVariables(self.population)

        content = tk.Frame(self.root)
        self._createContent(content)
        content.grid(row=0, column=0, sticky='news')
        gui.configureWeigths(content)

    def _loadVariables(self, population):
        self.observations = population.getAllParameters()
        mu, sigma, self.R, percentiles = population.getStats(self.observations)

        for i, varName in enumerate(population.modelParameters):
            v = PopulationVar(index=i,
                              varName=varName,
                              unitStr=PKPDUnit.codeToString(population.modelParameterUnits[i]),
                              mu=mu[i], sigma=sigma[i],
                              percentiles=[p[i] for p in percentiles])
            self._variables[varName] = v

    def _createContent(self, content):
        # Create top frame with the variables list
        self._createTopFrame(content)
        # Create button frame with Plot button
        self._createButtonsFrame(content)
        # Create an info frame
        self._createBottomFrame(content)

    def _createTopFrame(self, content):
        frame = tk.Frame(content)
        lfSamples = tk.LabelFrame(frame, text='Plot')
        gui.configureWeigths(frame)
        lfSamples.grid(row=0, column=0, sticky='news', padx=5, pady=5)

        # Create tree for Population Variables
        tp = PopulationVariablesTreeProvider(self._variables.values())
        self.tree = BoundTree(frame, tp, height=10)
        self.tree.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        gui.configureWeigths(frame)

        frame.grid(row=0, column=0, sticky='news', padx=5, pady=(10, 5))

    def _createButtonsFrame(self, content):
        frame = tk.Frame(content)
        self.plotButton = HotButton(frame, '   Plot   ', font=self.fontBold,
                                 command=self._onPlotClick,
                                 tooltip='Select one or two variables to plot ')

        self.plotButton.grid(row=0, column=0, sticky='se', padx=5)
        frame.grid(row=1, column=0, sticky='sew', padx=5, pady=5)
        gui.configureWeigths(frame)

    def _createBottomFrame(self, content):
        frame = tk.Frame(content)
        t = TaggedText(frame, height=10)
        t.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        t.addLine("*Correlation Matrix*")
        t.addText(np.array2string(self.R))
        frame.grid(row=2, column=0, sticky='sew', padx=5, pady=5)
        gui.configureWeigths(frame)

    def _onPlotClick(self, e=None):
        selection = self.tree.selection()
        n = len(selection)

        if n < 1 or n > 2:
            self.showError("Select one or two variables to plot.")
        else:
            plotter = EmPlotter(style='seaborn-whitegrid')
            varX = self._variables[selection[0]]
            xValues = self.observations[:, varX.index]

            def _label(var):
                return "%s [%s]" % (var.varName, var.unitStr)

            if n == 1:
                plotter.createSubPlot("Histogram", _label(varX), "Count")
                plotter.plotHist(xValues, 50)
            else: # n == 2
                varY = self._variables[selection[1]]
                yValues = self.observations[:, varY.index]
                ax = plotter.createSubPlot("Scatter Plot",
                                           _label(varX), _label(varY))
                ax.plot(xValues, yValues, '.')
            plotter.show()

