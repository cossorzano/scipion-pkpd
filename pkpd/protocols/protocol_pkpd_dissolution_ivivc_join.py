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

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.utils import computeXYmean, twoWayUniqueFloatValues
from .protocol_pkpd import ProtPKPD


class ProtPKPDDissolutionIVIVCJoin(ProtPKPD):
    """ Join several IVIVCs into a single one."""

    _label = 'dissol ivivc join'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputIVIVCs', params.MultiPointerParam, label="IVIVCs Fabs",
                      pointerClass='PKPDExperiment', help='Choose experiments with IVIV correlations (only the Fabs experiments)')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllIvIvC')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def calculateAllIvIvC(self):
        L1=[]
        L2=[]
        L3=[]
        for ptrExperiment in self.inputIVIVCs:
            experiment=PKPDExperiment()
            experiment.load(ptrExperiment.get().fnPKPD.get())
            x,y = experiment.getXYMeanValues("tvivo","tvitroReinterpolated")
            L1.append((x,y))
            x,y = experiment.getXYMeanValues("tvivo","Fabs")
            L2.append((x,y))
            x,y = experiment.getXYMeanValues("AdissolReinterpolated","FabsPredicted")
            L3.append((x,y))
        tvivo1, tvitroReinterpolated = computeXYmean(L1)
        tvivo2, Fabs = computeXYmean(L2)
        AdissolReinterpolated, FabsPredicted = computeXYmean(L3)

        tvivo1, tvitroReinterpolated = twoWayUniqueFloatValues(tvivo1, tvitroReinterpolated)
        tvivo2, Fabs = twoWayUniqueFloatValues(tvivo2, Fabs)
        AdissolReinterpolated, FabsPredicted = twoWayUniqueFloatValues(AdissolReinterpolated, FabsPredicted)

        Bt=InterpolatedUnivariateSpline(tvivo1, tvitroReinterpolated, k=1)
        BtA=InterpolatedUnivariateSpline(tvitroReinterpolated, AdissolReinterpolated, k=1)
        BtF=InterpolatedUnivariateSpline(tvivo2, Fabs, k=1)
        BAF=InterpolatedUnivariateSpline(AdissolReinterpolated, FabsPredicted, k=1)

        vtvitroReinterpolated = np.zeros(tvivo1.size)
        vAdissolReinterpolated = np.zeros(tvivo1.size)
        vFabs = np.zeros(tvivo1.size)
        vFabsPredicted = np.zeros(tvivo1.size)
        for i in range(tvivo1.size):
            tvivoi = tvivo1[i]
            vtvitroReinterpolated[i]=Bt(tvivoi)
            vAdissolReinterpolated[i]=BtA(vtvitroReinterpolated[i])
            vFabs[i]=BtF(tvivoi)
            vFabsPredicted[i]=BAF(vAdissolReinterpolated[i])

        tvitroReinterpolatedVar=experiment.variables["tvitroReinterpolated"]
        AdissolReinterpolatedVar=experiment.variables["AdissolReinterpolated"]
        tvivoVar=experiment.variables["tvivo"]
        FabsVar=experiment.variables["Fabs"]
        FabsPredictedVar=experiment.variables["FabsPredicted"]

        self.outputExperimentFabsSingle = PKPDExperiment()
        self.outputExperimentFabsSingle.variables[tvitroReinterpolatedVar.varName] = tvitroReinterpolatedVar
        self.outputExperimentFabsSingle.variables[AdissolReinterpolatedVar.varName] = AdissolReinterpolatedVar
        self.outputExperimentFabsSingle.variables[tvivoVar.varName] = tvivoVar
        self.outputExperimentFabsSingle.variables[FabsVar.varName] = FabsVar
        self.outputExperimentFabsSingle.variables[FabsPredictedVar.varName] = FabsPredictedVar
        self.outputExperimentFabsSingle.general["title"] = "In-vitro In-vivo correlation"
        self.outputExperimentFabsSingle.general["comment"] = "Fabs vs Predicted Fabs"

        sampleName="jointIVIVC"
        newSampleFabsSingle = PKPDSample()
        newSampleFabsSingle.sampleName = sampleName
        newSampleFabsSingle.variableDictPtr = self.outputExperimentFabsSingle.variables
        newSampleFabsSingle.descriptors = {}
        newSampleFabsSingle.addMeasurementColumn("tvitroReinterpolated", vtvitroReinterpolated)
        newSampleFabsSingle.addMeasurementColumn("AdissolReinterpolated", vAdissolReinterpolated)
        newSampleFabsSingle.addMeasurementColumn("tvivo", tvivo1)
        newSampleFabsSingle.addMeasurementColumn("FabsPredicted", vFabsPredicted)
        newSampleFabsSingle.addMeasurementColumn("Fabs",vFabs)
        self.outputExperimentFabsSingle.addLabelToSample(sampleName, "from", "individual---vesel", "AvgVivo---AvgVitro")

        self.outputExperimentFabsSingle.samples[sampleName] = newSampleFabsSingle
        self.outputExperimentFabsSingle.write(self._getPath("experimentFabsSingle.pkpd"))
    def createOutputStep(self):
        self._defineOutputs(outputExperimentFabsSingle=self.outputExperimentFabsSingle)
        for ptrExperiment in self.inputIVIVCs:
            self._defineSourceRelation(ptrExperiment.get(), self.outputExperimentFabsSingle)

    def _validate(self):
        retval = []
        for ptrExperiment in self.inputIVIVCs:
            if not "experimentFabs" in ptrExperiment.get().fnPKPD.get():
                retval.append("You can only take Fabs files")
        return retval

    def _summary(self):
        return []