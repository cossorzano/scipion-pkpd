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

from pkpd.protocols import ProtPKPDCompareExperiments

class PKPDCompareExperimentsViewer(Viewer):
    _targets = [ProtPKPDCompareExperiments]
    _environments = [DESKTOP_TKINTER]

    def getPlotValues(self,sample):
        xValues, yValues = sample.getXYValues(self.timeVarName, self.CVarName)
        return xValues[0], yValues[0] # From [array(...)] to array(...)

    def getSummary(self,fnPKPD,X,Y):
        experiment = PKPDExperiment()
        experiment.load(fnPKPD)
        self.timeVarName = X
        self.CVarName = Y

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
        X2=prot.X1.get() if prot.X2.get()=="" else prot.X2.get()
        Y2=prot.Y1.get() if prot.Y2.get()=="" else prot.Y2.get()

        sortedX1, Y1 = self.getSummary(prot.inputExperiment1.get().fnPKPD, prot.X1.get(), prot.Y1.get())
        sortedX2, Y2 = self.getSummary(prot.inputExperiment2.get().fnPKPD, X2, Y2)

        plotter = EmPlotter(style='seaborn-whitegrid')
        ax = plotter.createSubPlot("Summary Plot", self.timeVarName, self.CVarName)
        ax.plot(sortedX1, Y1[:, 0], 'r--', label="Minimum Exp1", linewidth=2)
        ax.plot(sortedX1, Y1[:, 1], 'b--', label="25% Exp1", linewidth=2)
        ax.plot(sortedX1, Y1[:, 2], 'g', label="50% (Median) Exp1", linewidth=2)
        ax.plot(sortedX1, Y1[:, 3], 'b--', label="75% Exp1", linewidth=2)
        ax.plot(sortedX1, Y1[:, 4], 'r--', label="Maximum Exp1", linewidth=2)

        ax.plot(sortedX2, Y2[:, 0], 'r--', label="Minimum Exp2")
        ax.plot(sortedX2, Y2[:, 1], 'b--', label="25% Exp2")
        ax.plot(sortedX2, Y2[:, 2], 'g', label="50% (Median) Exp2")
        ax.plot(sortedX2, Y2[:, 3], 'b--', label="75% Exp2")
        ax.plot(sortedX2, Y2[:, 4], 'r--', label="Maximum Exp2")

        ax.grid(True)
        ax.legend()
        plotter.show()
