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
from pyworkflow.protocol.constants import LEVEL_ADVANCED

class ProtPKPDBEPowerAnalysis(ProtPKPD):
    """ Power analysis for Bioequivalence studies.\n
        For further help, see https://cran.r-project.org/web/packages/PowerTOST/vignettes/PA.html
        https://cran.r-project.org/web/packages/PowerTOST/PowerTOST.pdf
        https://cran.r-project.org/web/packages/PowerTOST/index.html
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'BE power analysis'

    DESIGNSABE   = ['2x2x2','2x2x3','2x2x4','2x3x3', '2x4x2', '2x4x4', '3x6x3', '4x4']
    DESIGNSscABE = ['2x2x3','2x2x4','2x3x3']
    DESIGNSNTID  = ['2x2x3','2x2x4']
    REGULATORS = ['EMA', "HC", "FDA", "GCC"]

    METHOD_ABE = 0
    METHOD_scABE = 1
    METHOD_NTIDFDA = 2

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('method', params.EnumParam, label="Method",
                      choices=['Average bioequivalence', 'Scaled average bioequivalence',
                               'Narrow Therapeutic Index Drugs (FDA)'], default=self.METHOD_ABE,
                      help='See pa.ABE, pa.scABE, and pa.NTIDFDA in '
                           'https://cran.r-project.org/web/packages/PowerTOST/PowerTOST.pdf')
        form.addParam('regulator', params.EnumParam, label="Regulator", choices=self.REGULATORS, default=0,
                      condition='method==1', help='See pa.scABE in '
                           'https://cran.r-project.org/web/packages/PowerTOST/PowerTOST.pdf')
        form.addParam('designABE', params.EnumParam, label="Design",
                      choices=self.DESIGNSABE, default=0, condition='method==0',
                      help='Treatments x Periods x Sequences.\n'
                           'For scABE the valid designs are 2x3x3, 2x2x4, and 2x2x3.\n'
                           'See design in https://cran.r-project.org/web/packages/PowerTOST/vignettes/vignette.html')
        form.addParam('designscABE', params.EnumParam, label="Design",
                      choices=self.DESIGNSscABE, default=0, condition='method==1',
                      help='Treatments x Periods x Sequences.\n'
                           'See design in https://cran.r-project.org/web/packages/PowerTOST/vignettes/vignette.html')
        form.addParam('designNTID', params.EnumParam, label="Design",
                      choices=self.DESIGNSNTID, default=0, condition='method==2',
                      help='Treatments x Periods x Sequences.\n'
                           'See design in https://cran.r-project.org/web/packages/PowerTOST/vignettes/vignette.html')
        form.addParam('CV', params.FloatParam, label="Coefficient of variation", default=0.2,
                      help='See https://cran.r-project.org/web/packages/PowerTOST/vignettes/vignette.html')
        form.addParam('theta0', params.FloatParam, label="True or assumed T/R ratio", default=0.95,
                      help='See theta0 in https://cran.r-project.org/web/packages/PowerTOST/vignettes/vignette.html')
        form.addParam('minPower', params.FloatParam, label="Minimum power", default=0.7, expertLevel=LEVEL_ADVANCED,
                      help='See power in https://cran.r-project.org/web/packages/PowerTOST/vignettes/vignette.html')
        form.addParam('targetPower', params.FloatParam, label="Target power", default=0.8, expertLevel=LEVEL_ADVANCED,
                      help='See power in https://cran.r-project.org/web/packages/PowerTOST/vignettes/vignette.html')
        form.addParam('theta1', params.FloatParam, label="Lower limit", default=0.8, expertLevel=LEVEL_ADVANCED,
                      help='See theta1 in https://cran.r-project.org/web/packages/PowerTOST/vignettes/vignette.html')
        form.addParam('theta2', params.FloatParam, label="Upper limit", default=1.25, expertLevel=LEVEL_ADVANCED,
                      help='See theta2 in https://cran.r-project.org/web/packages/PowerTOST/vignettes/vignette.html')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runPA')

    #--------------------------- STEPS functions --------------------------------------------
    def runPA(self):
        fnPlot = self._getPath('plot.png')

        rScript = "library(PowerTOST)\n"
        extraArgs = ""
        design = ""
        if self.method==self.METHOD_ABE:
            rfunc = "pa.ABE"
            design = self.DESIGNSABE[self.designABE.get()]
        elif self.method==self.METHOD_scABE:
            rfunc = "pa.scABE"
            design = self.DESIGNSscABE[self.designscABE.get()]
            extraArgs = ', regulator="%s"'%self.REGULATORS[self.regulator.get()]
        else:
            rfunc = "pa.NTIDFDA"
            design = self.DESIGNSNTID[self.designNTID.get()]
        rScript+='x=%s(CV=%f, theta0=%f, theta1=%f, theta2=%f, minpower=%f, targetpower=%f, design="%s"%s)\n'%\
                 (rfunc, self.CV, self.theta0, self.theta1, self.theta2, self.minPower, self.targetPower,
                  design, extraArgs)
        rScript+="print(x)\n"
        rScript+='png(file="%s")\n'%fnPlot
        rScript+='plot(x)\n'
        rScript+='dummy=dev.off()\n'
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
