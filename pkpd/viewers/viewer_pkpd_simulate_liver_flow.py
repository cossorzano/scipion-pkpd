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

import numpy as np

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pyworkflow.em.viewers.plotter import EmPlotter

from pkpd.protocols import ProtPKPDSimulateLiverFlow

class PKPDSimulateLiverFlowViewer(Viewer):
    _targets = [ProtPKPDSimulateLiverFlow]
    _environments = [DESKTOP_TKINTER]


    def visualize(self, obj, **kwargs):
        prot = obj
        fnProfiles = prot._getPath("profiles.txt")
        fh = open(fnProfiles,"r")
        state = 0
        legends = ['Iliver', 'Iinlet', 'Isys']
        for line in fh:
            if state==0:
                tokens = line.split("::")
                title = tokens[1].strip()
                I3=[]
                state=1
            elif state==1:
                tokens=line.strip().split()
                if len(tokens)==0:
                    plotter = EmPlotter(style='seaborn-whitegrid')
                    ax = plotter.createSubPlot("Simulation", "t [h]", "[I] [mg/mL]")
                    t = np.asarray(I3[0],dtype=np.float64)/60
                    for n in range(1,len(I3)):
                        y = np.asarray(I3[n],dtype=np.float64)
                        ax.plot(t, y, label=legends[n-1])
                    ax.legend()
                    plotter.show()
                    state=0
                else:
                    if len(I3)==0:
                        for n in range(len(tokens)):
                            I3.append([])
                for n in range(len(tokens)):
                    I3[n].append(tokens[n])
        fh.close()
