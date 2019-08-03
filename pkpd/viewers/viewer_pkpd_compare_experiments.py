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

    def getData(self,fnPKPD,XvarName,YvarName):
        experiment = PKPDExperiment()
        experiment.load(fnPKPD)
        self.timeVarName = XvarName
        self.CVarName = YvarName

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
        if self.prot.showSummary.get():
            # We will store five values (min, 25%, 50%, 75%, max)
            # for each of the time entries computed
            percentileList = [0, 25, 50, 75, 100]
            Y = np.zeros((len(sortedTime), 5))
            for i, t in enumerate(sortedTime):
                Y[i, :] = np.percentile(dataDict[t], percentileList)
        else:
            Y=np.zeros((len(sortedTime),len(experiment.samples)))
            for i, t in enumerate(sortedTime):
                for j, C in enumerate(dataDict[t]):
                    Y[i, j] = C

        return sortedTime, Y


    def visualize(self, obj, **kwargs):
        self.prot = obj

        plotter = EmPlotter(style='seaborn-whitegrid')
        title = "Summary plot" if self.prot.showSummary else "Individual plots"
        ax = plotter.createSubPlot("Summary Plot", self.prot.X1.get(), self.prot.Y1.get())
        if self.prot.twoExperiments.get()==0:
            X2 = self.prot.X1.get() if self.prot.X2.get() == "" else self.prot.X2.get()
            Y2 = self.prot.Y1.get() if self.prot.Y2.get() == "" else self.prot.Y2.get()

            sortedX1, Y1 = self.getData(self.prot.inputExperiment1.get().fnPKPD, self.prot.X1.get(), self.prot.Y1.get())
            sortedX2, Y2 = self.getData(self.prot.inputExperiment2.get().fnPKPD, X2, Y2)

            if self.prot.showSummary:
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
            else:
                ax.plot(sortedX1, Y1, linewidth=2, label="Exp1")
                ax.plot(sortedX2, Y2, label="Exp2")

        else:
            listColors = ['b','g','r','k','c','m','y']
            for idx,experimentPtr in enumerate(self.prot.inputExperiments):
                sortedX, Y = self.getData(experimentPtr.get().fnPKPD, self.prot.X1.get(), self.prot.Y1.get())
                if self.prot.showSummary:
                    ax.plot(sortedX, Y[:, 0], linestyle='--', color=listColors[idx], label="Minimum Exp%d"%idx)
                    ax.plot(sortedX, Y[:, 1], linestyle='--', color=listColors[idx], label="25%% Exp%d"%idx)
                    ax.plot(sortedX, Y[:, 2], color=listColors[idx], label="50%% (Median) Exp%d"%idx)
                    ax.plot(sortedX, Y[:, 3], linestyle='--', color=listColors[idx], label="75%% Exp%d"%idx)
                    ax.plot(sortedX, Y[:, 4], linestyle='--', color=listColors[idx], label="Maximum Exp%d"%idx)
                else:
                    ax.plot(sortedX, Y[:,0], color=listColors[idx],label="Exp%d"%idx)
                    ax.plot(sortedX, Y[:,1:], color=listColors[idx])

        ax.legend()
        ax.grid(True)
        plotter.show()
