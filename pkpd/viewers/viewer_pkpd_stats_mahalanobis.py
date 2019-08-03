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

import os
from numpy import genfromtxt, isnan

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pyworkflow.em.viewers.plotter import EmPlotter

from pkpd.protocols import ProtPKPDStatsMahalanobis

class PKPDStatsMahalanobisViewer(Viewer):
    _targets = [ProtPKPDStatsMahalanobis]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        prot = obj
        fn = prot._getExtraPath('D11.txt')
        if os.path.exists(fn):
            D11 = genfromtxt(fn)
            plotter = EmPlotter(style='seaborn-whitegrid')
            plotter.createSubPlot("Histogram of D11", "D11", "Count")
            plotter.plotHist(D11[~isnan(D11)], 50)
            plotter.show()

        fn = prot._getExtraPath('D12.txt')
        if os.path.exists(fn):
            D12 = genfromtxt(fn)
            plotter = EmPlotter(style='seaborn-whitegrid')
            plotter.createSubPlot("Histogram of D12", "D12", "Count")
            plotter.plotHist(D12[~isnan(D12)], 50)
            plotter.show()
