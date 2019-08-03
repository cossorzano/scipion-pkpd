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

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pyworkflow.em.viewers.plotter import EmPlotter

from pkpd.protocols import ProtPKPDSimulateDrugInteractions


class PKPDSimulateDrugInteractionsViewer(Viewer):
    _targets = [ProtPKPDSimulateDrugInteractions]
    _environments = [DESKTOP_TKINTER]

    def addLimits(self,plotter,previousType,minX,maxX):
        if previousType=="ReversibleLiver":
            ax = plotter.getLastSubPlot()
            ax.plot([minX, maxX],[1.02,1.02],'b--', label="EMA")
            ax.plot([minX, maxX],[1.1,1.1],'r-.', label="FDA")
            leg=ax.legend()
            if leg:
                leg.draggable()
            plotter.draw()
        elif previousType=="ReversibleGut" or previousType=="TimeDependentGut":
            ax = plotter.getLastSubPlot()
            ax.plot([minX, maxX],[11,11],'r--', label="EMA/FDA")
            leg=ax.legend()
            if leg:
                leg.draggable()
            plotter.draw()
        elif previousType=="TimeDependentLiver":
            ax = plotter.getLastSubPlot()
            ax.plot([minX, maxX],[1.25,1.25],'b--', label="EMA")
            ax.plot([minX, maxX],[1.1,1.1],'r-.', label="FDA")
            leg=ax.legend()
            if leg:
                leg.draggable()
            plotter.draw()
        elif previousType=="InductionLiver":
            ax = plotter.getLastSubPlot()
            ax.plot([minX, maxX],[0.9,0.9],'r--', label="FDA")
            leg=ax.legend()
            if leg:
                leg.draggable()
            plotter.draw()
        elif previousType=="TransporterGut":
            ax = plotter.getLastSubPlot()
            ax.plot([minX, maxX],[11,11],'r--', label="EMA/FDA")
            leg=ax.legend()
            if leg:
                leg.draggable()
            plotter.draw()
        elif previousType=="TransporterLiver":
            ax = plotter.getLastSubPlot()
            ax.plot([minX, maxX],[1.04,1.04],'b--', label="EMA")
            ax.plot([minX, maxX],[1.1,1.1],'r--', label="FDA")
            leg=ax.legend()
            if leg:
                leg.draggable()
            plotter.draw()
        elif previousType=="TransporterRenal":
            ax = plotter.getLastSubPlot()
            ax.plot([minX, maxX],[1.02,1.02],'b--', label="EMA")
            ax.plot([minX, maxX],[1.1,1.1],'r--', label="FDA")
            leg=ax.legend()
            if leg:
                leg.draggable()
            plotter.draw()

    def visualize(self, obj, **kwargs):
        prot = obj
        fnProfiles = prot._getPath("profiles.txt")
        fh = open(fnProfiles,"r")
        Rtype = []
        Rlegends = []
        R = []
        state = 0
        for line in fh:
            if state==0:
                tokens = line.split("::")
                Rtype.append(tokens[0])
                Rlegends.append(tokens[1].strip())
                Ri=[]
                state=1
            elif state==1:
                tokens=line.strip().split()
                if len(tokens)==0:
                    R.append(Ri)
                    state=0
                else:
                    if len(Ri)==0:
                        for n in range(len(tokens)):
                            Ri.append([])
                for n in range(len(tokens)):
                    Ri[n].append(tokens[n])
        fh.close()

        plotter = None
        previousType = ""
        for legend, Ri, Rtypei  in izip(Rlegends, R, Rtype):
            if plotter is None or Rtypei!=previousType:
                if previousType!="":
                    self.addLimits(plotter,previousType,minX,maxX)
                plotter = EmPlotter(style='seaborn-whitegrid')
                doShow = True
                if Rtypei=="ReversibleLiver" or Rtypei=="TimeDependentLiver" or Rtypei=="InductionLiver" or Rtypei=="StaticLiver" or Rtypei=="TransporterLiver":
                    Ilabel="[Ih] [uM]"
                elif Rtypei=="ReversibleGut" or Rtypei=="TimeDependentGut" or Rtypei=="InductionGut" or Rtypei=="StaticGut" or Rtypei=="TransporterGut":
                    Ilabel="[Ig] [uM]"
                elif Rtypei=="TransporterRenal":
                    Ilabel="[Cmax] [uM]"
                ax = plotter.createSubPlot("Plot", Ilabel, "R")
                previousType = Rtypei
                minX = None
                maxX = None
            else:
                doShow = False
                ax = plotter.getLastSubPlot()

            x = np.asarray(Ri[0],dtype=np.float64)
            if len(Ri)==2:
                y=np.asarray(Ri[1],dtype=np.float64)
            else:
                y=np.asarray(Ri[2],dtype=np.float64)
            ax.plot(x, y, label=legend)

            minXi = np.min(x)
            maxXi = np.max(x)
            if minX==None:
                minX=minXi
                maxX=maxXi
            minX=min(minXi,minX)
            maxX=min(maxXi,maxX)

            leg = ax.legend()
            if leg:
                leg.draggable()

            if doShow:
                plotter.show()
            else:
                plotter.draw()
        self.addLimits(plotter,previousType,minX,maxX)
