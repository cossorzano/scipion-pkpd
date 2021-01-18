# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (info@kinestat.com)
# *
# * Kinestat Pharma
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
# *  e-mail address 'info@kinestat.com'
# *
# **************************************************************************

from itertools import izip
import math
import numpy as np
from scipy.stats import norm

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD


# tested in test_workflow_dissolution_f2.py

class ProtPKPDParticleSize(ProtPKPD):
    """ Estimate the particle size distribution parameters. The size is supposed to be log-normal
        and the mean and standard deviation of the log-normal distribution are estimated. Given a set
        of particle sizes {x1, x2, ..., xN} and a set of bin occupancies (they can be in absolute or
        relative counts or mass) {p1, p2, ..., pN}, this protocol estimates the parameters of the
        log-normal distribution that is compatible with these observations. It also reports
        the fitting quality of the distribution.

        If it is descending order, then it is assumed that p1=Prob{x1<=x<x2}, p2=Prob{x2<=x<x3}, ..., pN=Prob{x<=xN}
        If it is ascending order, then it is assumed that p1=Prob{x<=x1}, p2=Prob{x1<x<=x2}, ..., pN=Prob{xN<=x}

        The program calculates the mu and sigma of the log-normal distribution and their exponentials.
        """

    _label = 'particle size distribution'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('x', params.StringParam, label="Particle size (x)", default="",
                      help='Separated by commas, example: x1, x2, x3')
        form.addParam('p', params.StringParam, label="Mass, counts or probability (p)", default="",
                      help='The probability of each bin is pi/sum(pi). Separated by commas, example: p1, p2, p3')
        form.addParam('descending', params.BooleanParam, label="Descending", default=True,
                      help='If it is descending order, then it is assumed that p1=Prob{x1<=x<x2}, p2=Prob{x2<=x<x3}, ..., pN=Prob{x<=xN}\n'
                           'If it is ascending order, then it is assumed that p1=Prob{x<=x1}, p2=Prob{x1<x<=x2}, ..., pN=Prob{xN<=x}')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculate')

    #--------------------------- STEPS functions --------------------------------------------
    def calculate(self):
        x = np.asarray([float(xi.strip()) for xi in self.x.get().split(',')])
        p = np.asarray([float(xi.strip()) for xi in self.p.get().split(',')])
        if self.descending:
            x=np.flip(x,0)
            p=np.flip(p,0)

        P = np.cumsum(p/np.sum(p))
        z = norm.ppf(P)

        A = np.vstack([np.ones(len(z[:-1])), z[:-1]]).T
        b = np.log(x[:-1])
        mu, sigma = np.linalg.lstsq(A, b, rcond=None)[0]
        fhSummary = open(self._getPath("summary.txt"),"w")
        self.doublePrint(fhSummary, "Log-normal mean mu=%f exp(mu)=%f"%(mu,math.exp(mu)))
        self.doublePrint(fhSummary, "Log-normal mean sigma=%f exp(sigma)=%f"%(sigma,math.exp(sigma)))

        bp = A.dot(np.array([mu, sigma]))
        residuals = b-bp
        R2=1-np.var(residuals)/np.var(b)
        self.doublePrint(fhSummary, "R2=%f" % R2)
        fhSummary.close()
        print(" ")


        residuals = np.append(residuals,np.NaN)
        i=1
        for xi, pi, Pi, zi, ei in izip(x, p, P, z, residuals):
            print("x%d=%f (log=%f) p%d=%f Cumulated fraction (P%d)=%f (z=%f, residual=%f)"%\
                  (i,xi,math.log(xi),i,pi,i,Pi,zi,ei))
            i+=1

    def _summary(self):
        retval = []
        self.addFileContentToMessage(retval,self._getPath("summary.txt"))
        return retval
