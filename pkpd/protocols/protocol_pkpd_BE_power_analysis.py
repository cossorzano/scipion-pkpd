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

import os

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd import Plugin as pkpdPlugin


class ProtPKPDBEPowerAnalysis(ProtPKPD):
    """ Power analysis for Bioequivalence studies.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'BE power analysis'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('labelToAdd', params.StringParam, label="Label to add", default="CumulatedDose",
                      help='Name of the variable to add')
        form.addParam('fromTime', params.FloatParam, label="From (h)", default=0)
        form.addParam('toTime', params.FloatParam, label="To (h)", default=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runPA')

    #--------------------------- STEPS functions --------------------------------------------
    def runPA(self):
        fnPlot = self._getPath('plot.png')

        rScript = """
        library(PowerTOST)
        x=pa.ABE(CV = 0.20, theta0 = 1)
        print(x)
        png(file="%s")
        plot(x)
        dummy=dev.off()
        """%fnPlot
        response=pkpdPlugin.runRscript(rScript)
        print("Running in R ======")
        print(rScript)
        print("R response ========")
        print(response)

        fh = open(self._getPath('summary.txt'),'w')
        fh.write(response)
        fh.close()

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        retval = []
        if pkpdPlugin.getRscript() is None:
            retval.append("Cannot find Rscript in the system PATH")
        else:
            if not pkpdPlugin.checkRPackage('PowerTOST'):
                retval.append("PowerTOST package is not installed in R")
        return retval

    def _summary(self):
        msg=[]
        fnSummary = self._getPath('summary.txt')
        if os.path.exists(fnSummary):
            fh = open(fnSummary,'r')
            msg = [x.strip() for x in fh.readlines()]
            fh.close()
        return msg
