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

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pyworkflow.gui.text import openTextFileEditor

from pkpd.protocols import ProtPKPDBEPowerAnalysis

class PKPDBEPowerAnalysisViewer(Viewer):
    _targets = [ProtPKPDBEPowerAnalysis]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        prot = obj
        fnPlot = prot._getPath('plot.png')
        if os.path.exists(fnPlot):
            openTextFileEditor(fnPlot)