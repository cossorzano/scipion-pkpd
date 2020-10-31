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
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pyworkflow.em.viewers.plotter import EmPlotter

from pkpd.protocols import ProtPKPDParticleSize

class PKPDParticleSizeViewer(Viewer):
    _targets = [ProtPKPDParticleSize]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        self.prot = obj

        # Experimental
        x = np.asarray([float(xi.strip()) for xi in self.prot.x.get().split(',')])
        p = np.asarray([float(xi.strip()) for xi in self.prot.p.get().split(',')])
        if self.prot.descending:
            x=np.flip(x,0)
            p=np.flip(p,0)
        logx = np.log(x)
        logx = np.insert(logx, 0, np.log(np.min(x)/2))
        p=p/np.sum(p)

        barx = []
        bary = p
        widths = []
        locx = []
        labels = []

        for i in range (0,p.shape[0]):
            barx.append(0.5*(logx[i]+logx[i+1]))
            widths.append(logx[i+1]-logx[i])
            locx.append(np.log(x[i]))
            labels.append("log(%4.2f)"%x[i])

        plotter =EmPlotter(style='seaborn-whitegrid')
        ax = plotter.createSubPlot("Particle size distribution", "log(Particle size)", "Fraction")
        ax.bar(barx, bary, width=widths, linewidth=1, label="Experimental", edgecolor="black")
        plt.xticks(locx, labels)

        # Theoretical
        fhSummary = open(self.prot._getPath("summary.txt"))
        lineno=0
        for line in fhSummary.readlines():
            if lineno==0:
                mu = float((line.split()[2]).split('=')[1])
            elif lineno==1:
                sigma = float((line.split()[2]).split('=')[1])
            lineno+=1
        fhSummary.close()

        logx = np.arange(np.min(logx),np.max(logx),(np.max(logx)-np.min(logx))/100)
        theox = norm.pdf(logx,mu,sigma)
        ax.plot(logx,theox/np.max(theox)*np.max(bary), color='red', label='Theoretical (log-normal)')

        # General
        ax.legend()
        ax.grid(True)
        plotter.show()
