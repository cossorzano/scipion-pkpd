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

from itertools import izip
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pyworkflow.em.viewers.plotter import EmPlotter
from pkpd.objects import PKPDExperiment
from pkpd.utils import uniqueFloatValues

from pkpd.protocols import ProtPKPDIVIVCInternalValidity

class PKPDDissolutionIVIVCInternalValidityViewer(Viewer):
    _targets = [ProtPKPDIVIVCInternalValidity]
    _environments = [DESKTOP_TKINTER]

    def getPlotValues(self,sample):
        xValues, yValues = sample.getXYValues(self.timeVarName, self.CVarName)
        return xValues[0], yValues[0] # From [array(...)] to array(...)

    def getSummary(self,fnPKPD):
        experiment = PKPDExperiment()
        experiment.load(fnPKPD)
        self.timeVarName = experiment.getTimeVariable()
        self.CVarName = experiment.getMeasurementVariables()[0] # The first one

        xmin = 1e38
        xmax = -1e38
        for sampleName, sample in experiment.samples.iteritems():
            xValues, _ = self.getPlotValues(sample)
            xmin = min(xmin, min(xValues))
            xmax = max(xmax, max(xValues))

        dataDict = {}  # key will be time values
        xrange = np.arange(xmin, xmax, (xmax - xmin) / 300.0)
        for sampleName, sample in experiment.samples.iteritems():
            xValues, yValues = self.getPlotValues(sample)
            xValuesUnique, yValuesUnique = uniqueFloatValues(xValues, yValues)
            B = InterpolatedUnivariateSpline(xValuesUnique, yValuesUnique, k=1)
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
            Y[i, :] = np.percentile(dataDict[t], percentileList)
        return sortedTime, Y


    def visualize(self, obj, **kwargs):
        prot = obj
        auc = np.genfromtxt(prot._getExtraPath('errorAUC.txt'))
        cmax = np.genfromtxt(prot._getExtraPath('errorCmax.txt'))

        plotter = EmPlotter(style='seaborn-whitegrid')
        plotter.createSubPlot("Histogram of error AUC0t", "Error AUC0t", "Count")
        plotter.plotHist(auc[~np.isnan(auc)], 50)
        plotter.show()

        plotter = EmPlotter(style='seaborn-whitegrid')
        plotter.createSubPlot("Histogram of error Cmax", "Error Cmax", "Count")
        plotter.plotHist(cmax[~np.isnan(cmax)], 50)
        plotter.show()

        sortedTimeTrue, Ytrue = self.getSummary(prot.inputExperiment.get().fnPKPD)
        sortedTimeSimulated, Ysimulated = self.getSummary(prot.inputSimulated.get().fnPKPD)

        plotter = EmPlotter(style='seaborn-whitegrid')
        ax = plotter.createSubPlot("Summary Plot", self.timeVarName, self.CVarName)
        ax.plot(sortedTimeTrue, Ytrue[:, 0], 'r--', label="Minimum In-vivo", linewidth=2)
        ax.plot(sortedTimeTrue, Ytrue[:, 1], 'b--', label="25%  In-vivo", linewidth=2)
        ax.plot(sortedTimeTrue, Ytrue[:, 2], 'g', label="50% (Median)  In-vivo", linewidth=2)
        ax.plot(sortedTimeTrue, Ytrue[:, 3], 'b--', label="75%  In-vivo", linewidth=2)
        ax.plot(sortedTimeTrue, Ytrue[:, 4], 'r--', label="Maximum  In-vivo", linewidth=2)
        ax.plot(sortedTimeSimulated, Ysimulated[:, 0], 'r--', label="Minimum Simulated")
        ax.plot(sortedTimeSimulated, Ysimulated[:, 1], 'b--', label="25% Simulated")
        ax.plot(sortedTimeSimulated, Ysimulated[:, 2], 'g', label="50% (Median) Simulated")
        ax.plot(sortedTimeSimulated, Ysimulated[:, 3], 'b--', label="75% Simulated")
        ax.plot(sortedTimeSimulated, Ysimulated[:, 4], 'r--', label="Maximum Simulated")
        ax.grid(True)
        ax.legend()
        plotter.show()

        plotter = EmPlotter(style='seaborn-whitegrid')
        ax = plotter.createSubPlot("Mean Plot", self.timeVarName, self.CVarName)
        ax.plot(sortedTimeTrue, Ytrue[:, 2], 'g', label="50% (Median)  In-vivo", linewidth=2)
        ax.plot(sortedTimeSimulated, Ysimulated[:, 2], 'g', label="50% (Median) Simulated")
        ax.grid(True)
        ax.legend()
        plotter.show()
