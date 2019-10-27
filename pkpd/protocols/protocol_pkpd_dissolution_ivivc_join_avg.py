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


class ProtPKPDDissolutionIVIVCJoin(ProtPKPD):
    """ Join several IVIVCs into a single one. The strategy is to compute the average of all the plots involved in the
        IVIVC process: 1) tvivo -> tvitro; 2) tvitro -> Adissol; 3) Adissol->FabsPredicted. The plot tvivo-Fabs comes
        after the IVIVC process, while the plot tvivo-FabsOrig is the observed one in the input files. These two
        plots need not be exactly the same. """

    _label = 'dissol ivivc join avg'

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
        L4=[]
        L5=[]
        for ptrExperiment in self.inputIVIVCs:
            experiment=PKPDExperiment()
            experiment.load(ptrExperiment.get().fnPKPD.get())
            x,y = experiment.getXYMeanValues("tvivo","tvitroReinterpolated")
            L1.append((x,y))
            x,y = experiment.getXYMeanValues("tvivo","Fabs")
            L2.append((x,y))
            x,y = experiment.getXYMeanValues("AdissolReinterpolated","FabsPredicted")
            L3.append((x,y))
            x, y = experiment.getXYMeanValues("FabsPredicted", "Fabs")
            L4.append((x, y))
            x, y = experiment.getXYMeanValues("tvitroReinterpolated", "AdissolReinterpolated")
            L5.append((x, y))
        tvivo1, tvitroReinterpolatedY = computeXYmean(L1)
        tvivoOrig, FabsOrig = computeXYmean(L2)
        AdissolReinterpolatedX, FabsPredictedY = computeXYmean(L3)
        FabsPredictedX, Fabs = computeXYmean(L4)
        tvitroReinterpolatedX, AdissolReinterpolatedY = computeXYmean(L5)

        x, y = twoWayUniqueFloatValues(tvivo1, tvitroReinterpolatedY)
        Bt=InterpolatedUnivariateSpline(x, y, k=1)
        x, y = twoWayUniqueFloatValues(tvitroReinterpolatedX, AdissolReinterpolatedY)
        BtA=InterpolatedUnivariateSpline(x, y, k=1)
        x, y = twoWayUniqueFloatValues(tvivoOrig, FabsOrig)
        BtF=InterpolatedUnivariateSpline(x, y, k=1)
        x, y=twoWayUniqueFloatValues(AdissolReinterpolatedX, FabsPredictedY)
        BAF=InterpolatedUnivariateSpline(x, y, k=1)
        x, y= twoWayUniqueFloatValues(FabsPredictedX, Fabs)
        BFF=InterpolatedUnivariateSpline(x, y, k=1)

        vtvitroReinterpolated = np.zeros(len(tvivo1))
        vAdissolReinterpolated = np.zeros(len(tvivo1))
        vFabs = np.zeros(len(tvivo1))
        vFabsPredicted = np.zeros(len(tvivo1))
        vFabsOrig = np.zeros(len(tvivo1))
        for i in range(len(tvivo1)):
            tvivoi = tvivo1[i]
            vtvitroReinterpolated[i]=Bt(tvivoi)
            vAdissolReinterpolated[i]=BtA(vtvitroReinterpolated[i])
            vFabsPredicted[i]=BAF(vAdissolReinterpolated[i])
            vFabs[i]=BFF(vFabsPredicted[i])
            vFabsOrig[i]=BtF(tvivoi)

        tvitroReinterpolatedVar=experiment.variables["tvitroReinterpolated"]
        AdissolReinterpolatedVar=experiment.variables["AdissolReinterpolated"]
        tvivoVar=experiment.variables["tvivo"]
        FabsOrigVar=copy.copy(experiment.variables["Fabs"])
        FabsOrigVar.varName = "FabsOriginal"
        FabsVar=experiment.variables["Fabs"]
        FabsVar.comment += ". After IVIVC: tvivo->tvitro->Adissol->Fabs "
        FabsPredictedVar=experiment.variables["FabsPredicted"]

        self.outputExperimentFabsSingle = PKPDExperiment()
        self.outputExperimentFabsSingle.variables[tvitroReinterpolatedVar.varName] = tvitroReinterpolatedVar
        self.outputExperimentFabsSingle.variables[AdissolReinterpolatedVar.varName] = AdissolReinterpolatedVar
        self.outputExperimentFabsSingle.variables[tvivoVar.varName] = tvivoVar
        self.outputExperimentFabsSingle.variables[FabsVar.varName] = FabsVar
        self.outputExperimentFabsSingle.variables[FabsPredictedVar.varName] = FabsPredictedVar
        self.outputExperimentFabsSingle.variables[FabsOrigVar.varName] = FabsOrigVar
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
        newSampleFabsSingle.addMeasurementColumn("FabsOriginal",vFabsOrig)

        self.outputExperimentFabsSingle.samples[sampleName] = newSampleFabsSingle
        self.outputExperimentFabsSingle.addLabelToSample(sampleName, "from", "individual---vesel", "AvgVivo---AvgVitro")

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