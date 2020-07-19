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

import copy
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.utils import computeXYmean, twoWayUniqueFloatValues
from .protocol_pkpd import ProtPKPD


class ProtPKPDDissolutionLevyPlotJoin(ProtPKPD):
    """ Join several Levy plots into a single one. The strategy is to compute the average of all the plots involved in the
        Levy plot process: 1) tvivo -> tvitro """

    _label = 'dissol levyplot join avg'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputLevyplots', params.MultiPointerParam, label="Levy plots",
                      pointerClass='PKPDExperiment', help='Choose experiments with IVIV correlations (only the Fabs experiments)')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllLevy')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def calculateAllLevy(self):
        L1=[]
        for ptrExperiment in self.inputLevyplots:
            experiment=PKPDExperiment()
            experiment.load(ptrExperiment.get().fnPKPD.get())
            x,y = experiment.getXYMeanValues("tvivo","tvitro")
            L1.append((x,y))

        tvivo, tvitro = computeXYmean(L1, common=True)

        x, y = twoWayUniqueFloatValues(tvivo, tvitro)
        Bt=InterpolatedUnivariateSpline(x, y, k=1)

        vtvitroReinterpolated = np.zeros(len(tvivo))
        for i in range(len(tvivo)):
            tvivoi = tvivo[i]
            vtvitroReinterpolated[i]=Bt(tvivoi)

        tvitroReinterpolatedVar=experiment.variables["tvitro"]
        tvivoVar=experiment.variables["tvivo"]

        self.outputExperimentSingle = PKPDExperiment()
        self.outputExperimentSingle.variables[tvitroReinterpolatedVar.varName] = tvitroReinterpolatedVar
        self.outputExperimentSingle.variables[tvivoVar.varName] = tvivoVar
        self.outputExperimentSingle.general["title"] = "Levy plot"
        self.outputExperimentSingle.general["comment"] = "tvitro vs tvivo"

        sampleName="jointLevyPlot"
        newSampleSingle = PKPDSample()
        newSampleSingle.sampleName = sampleName
        newSampleSingle.variableDictPtr = self.outputExperimentSingle.variables
        newSampleSingle.descriptors = {}
        newSampleSingle.addMeasurementColumn("tvitro", vtvitroReinterpolated)
        newSampleSingle.addMeasurementColumn("tvivo", tvivo)

        self.outputExperimentSingle.samples[sampleName] = newSampleSingle

        self.outputExperimentSingle.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperimentSingle)
        for ptrExperiment in self.inputLevyplots:
            self._defineSourceRelation(ptrExperiment.get(), self.outputExperimentSingle)

    def _summary(self):
        return []