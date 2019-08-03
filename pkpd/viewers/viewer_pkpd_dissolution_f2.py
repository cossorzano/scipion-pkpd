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

from numpy import genfromtxt, isnan

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pyworkflow.em.viewers.plotter import EmPlotter

from pkpd.protocols import ProtPKPDDissolutionF2

class PKPDDissolutionF2Viewer(Viewer):
    _targets = [ProtPKPDDissolutionF2]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        prot = obj
        f1 = genfromtxt(prot._getExtraPath('f1.txt'))
        f2 = genfromtxt(prot._getExtraPath('f2.txt'))

        plotter = EmPlotter(style='seaborn-whitegrid')
        plotter.createSubPlot("Histogram of f1", "f1", "Count")
        plotter.plotHist(f1[~isnan(f1)], 50)
        plotter.show()

        plotter = EmPlotter(style='seaborn-whitegrid')
        plotter.createSubPlot("Histogram of f2", "f2", "Count")
        plotter.plotHist(f2[~isnan(f2)], 50)
        plotter.show()
