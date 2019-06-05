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

from pkpd.protocols import ProtPKPDIVIVCInternalValidity

class PKPDDissolutionIVIVCInternalValidityViewer(Viewer):
    _targets = [ProtPKPDIVIVCInternalValidity]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        prot = obj
        auc = genfromtxt(prot._getExtraPath('errorAUC.txt'))
        cmax = genfromtxt(prot._getExtraPath('errorCmax.txt'))

        plotter = EmPlotter()
        plotter.createSubPlot("Histogram of error AUC0t", "Error AUC0t", "Count")
        plotter.plotHist(auc[~isnan(auc)], 50)
        plotter.show()

        plotter = EmPlotter()
        plotter.createSubPlot("Histogram of error Cmax", "Error Cmax", "Count")
        plotter.plotHist(cmax[~isnan(cmax)], 50)
        plotter.show()
