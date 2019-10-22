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
from itertools import izip

from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.utils import uniqueFloatValues, twoWayUniqueFloatValues
from pkpd.pkpd_units import createUnit
from pkpd.utils import computeXYmean
from .protocol_pkpd import ProtPKPD


# tested in test_workflow_levyplot
# tested in test_workflow_deconvolution2

class ProtPKPDDissolutionLevyPlot(ProtPKPD):
    """ Calculate the Levy plot between two dissolution experiments. Each experiment may have
        several profiles and all vs all profiles are calculated"""

    _label = 'dissol Levy'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputInVitro', params.PointerParam, label="Dissolution profiles in vitro",
                      pointerClass='ProtPKPDDissolutionFit', help='Select an experiment with dissolution profiles')
        form.addParam('inputInVivo', params.PointerParam, label="Dissolution profiles in vivo",
                      pointerClass='ProtPKPDDeconvolve, ProtPKPDDeconvolutionWagnerNelson, ProtPKPDDeconvolutionLooRiegelman, ProtPKPDDeconvolveFourier, PKPDExperiment',
                      help='Select an experiment with dissolution profiles')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllLevy',self.inputInVitro.get().getObjId(),self.inputInVivo.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def getInVivoProfiles(self):
        if hasattr(self.inputInVivo.get(),"outputExperiment"):
            fnPKPD = self.inputInVivo.get().outputExperiment.fnPKPD
        elif hasattr(self.inputInVivo.get(),"fnPKPD"):
            fnPKPD = self.inputInVivo.get().fnPKPD
        else:
            raise Exception("Cannot find a suitable filename for reading the experiment")
        experiment = self.readExperiment(fnPKPD)
        allY = []
        sampleNames = []
        for sampleName, sample in experiment.samples.iteritems():
            t=sample.getValues("t")
            y=sample.getValues("A")
            allY.append((np.asarray(t,dtype=np.float64),np.asarray(y,dtype=np.float64)))
            sampleNames.append(sampleName)
        self.experimentInVivo = experiment
        return allY, sampleNames

    def getInVitroModels(self):
        allParameters = []
        self.protFit = self.inputInVitro.get()
        experiment = self.readExperiment(self.protFit.outputExperiment.fnPKPD)
        self.fitting = self.readFitting(self.protFit.outputFitting.fnFitting)
        self.varNameX = self.fitting.predictor.varName
        self.varNameY = self.fitting.predicted.varName

        self.protFit.model = self.protFit.createModel()
        self.protFit.model.setExperiment(experiment)
        self.protFit.model.setXVar(self.varNameX)
        self.protFit.model.setYVar(self.varNameY)
        self.protFit.setupFromFormParameters()
        self.protFit.experiment = experiment
        self.protFit.setupModel()

        parameterNames = self.protFit.model.getParameterNames()
        vesselNames=[]
        for sampleName, sample in experiment.samples.iteritems():
            vesselNames.append(sampleName)
            parameters0 = []
            for parameterName in parameterNames:
                parameters0.append(float(sample.descriptors[parameterName]))
            allParameters.append(parameters0)
        self.experimentInVitro = experiment
        return allParameters, vesselNames

    def produceLevyPlot(self,tvivo,parameterInVitro,Avivo):
        Avivounique, tvivoUnique=uniqueFloatValues(Avivo,tvivo)
        B = InterpolatedUnivariateSpline(Avivounique, tvivoUnique, k=1)
        tvitro = np.arange(0, 10*np.max(tvivo), 1)
        self.protFit.model.x=tvitro
        Avitro = self.protFit.model.forwardModel(parameterInVitro)[0]
        tvivo=[]
        for i in range(tvitro.shape[0]):
            tvivo.append(B(Avitro[i]))
        return (tvitro,np.asarray(tvivo,dtype=np.float64),Avitro)

    def addSample(self, outputExperiment, sampleName, tvitro, tvivo, individualFrom, vesselFrom):
        newSample = PKPDSample()
        newSample.sampleName = sampleName
        newSample.variableDictPtr = self.outputExperiment.variables
        newSample.descriptors = {}
        newSample.addMeasurementColumn("tvivo",tvivo)
        newSample.addMeasurementColumn("tvitro", tvitro)
        outputExperiment.samples[sampleName] = newSample
        outputExperiment.addLabelToSample(sampleName, "from", "individual---vesel", "%s---%s"%(individualFrom,vesselFrom))

    def makeSuggestions(self, levyName, tvitro, tvivo):
        print("Polynomial fitting suggestions for %s"%levyName)
        for degree in range(1,10):
            coeffs =np.polyfit(tvivo,tvitro,degree)
            p=np.poly1d(coeffs, variable='tvivo')
            residuals=tvitro-p(tvivo)
            R2=1-np.var(residuals)/np.var(tvivo)
            print("Degree %d, R2: %f, tvitro="%(degree,R2))
            print(p)
            print(" ")

        logtvivo=np.log10(tvivo)
        logtvitro=np.log10(tvitro)
        idx=np.logical_and(np.isfinite(logtvitro),np.isfinite(logtvivo))
        logtvivo=logtvivo[idx]
        logtvitro=logtvitro[idx]
        for degree in range(1,10):
            coeffs =np.polyfit(logtvivo,logtvitro,degree)
            p=np.poly1d(coeffs, variable='log10(tvivo)')
            residuals=logtvitro-p(logtvivo)
            R2=1-np.var(residuals)/np.var(logtvivo)
            print("Degree %d, R2: %f, log10(tvitro)="%(degree,R2))
            print(p)
            print(" ")

    def calculateAllLevy(self, objId1, objId2):
        parametersInVitro, vesselNames =self.getInVitroModels()
        profilesInVivo, sampleNames =self.getInVivoProfiles()

        self.outputExperiment = PKPDExperiment()
        tvitrovar = PKPDVariable()
        tvitrovar.varName = "tvitro"
        tvitrovar.varType = PKPDVariable.TYPE_NUMERIC
        tvitrovar.role = PKPDVariable.ROLE_MEASUREMENT
        tvitrovar.units = createUnit("min")

        tvivovar = PKPDVariable()
        tvivovar.varName = "tvivo"
        tvivovar.varType = PKPDVariable.TYPE_NUMERIC
        tvivovar.role = PKPDVariable.ROLE_MEASUREMENT
        tvivovar.units = createUnit("min")

        self.outputExperiment.variables[tvivovar.varName] = tvivovar
        self.outputExperiment.variables[tvitrovar.varName] = tvitrovar
        self.outputExperiment.general["title"] = "Levy plots"
        self.outputExperiment.general["comment"] = "Time in vivo vs time in vitro"

        i=1
        levyList = []
        for parameterInVitro, vesselFrom in izip(parametersInVitro,vesselNames):
            for aux, sampleFrom in izip(profilesInVivo,sampleNames):
                t, profileInVivo = aux
                tvitro, tvivo, _ = self.produceLevyPlot(t,parameterInVitro,profileInVivo)
                tvitroUnique, tvivoUnique = twoWayUniqueFloatValues(tvitro, tvivo)
                idx = np.logical_and(tvitroUnique>0,tvivoUnique>0)
                tvitroUnique=tvitroUnique[idx]
                tvivoUnique=tvivoUnique[idx]
                if tvitroUnique[0]>0 and tvivoUnique[0]>0:
                    tvitroUnique=np.insert(tvitroUnique,0,0.0)
                    tvivoUnique = np.insert(tvivoUnique, 0, 0.0)
                levyName = "levy_%04d"%i
                levyList.append((tvitroUnique, tvivoUnique))
                self.addSample(self.outputExperiment, levyName, tvitroUnique, tvivoUnique, sampleFrom, vesselFrom)

                self.makeSuggestions(levyName, tvitroUnique, tvivoUnique)
                i+=1

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

        # Single Levy plot
        self.outputExperimentSingle = PKPDExperiment()
        self.outputExperimentSingle.variables[tvivovar.varName] = tvivovar
        self.outputExperimentSingle.variables[tvitrovar.varName] = tvitrovar
        self.outputExperimentSingle.general["title"] = "Levy plots"
        self.outputExperimentSingle.general["comment"] = "Time in vivo vs time in vitro"

        tvitroSingle, tvivoSingle = computeXYmean(levyList)
        tvitroUnique, tvivoUnique = twoWayUniqueFloatValues(tvitroSingle, tvivoSingle)
        idx = np.logical_and(tvitroUnique > 0, tvivoUnique > 0)
        tvitroUnique = tvitroUnique[idx]
        tvivoUnique = tvivoUnique[idx]
        if tvitroUnique[0] > 0 and tvivoUnique[0] > 0:
            tvitroUnique = np.insert(tvitroUnique, 0, 0.0)
            tvivoUnique = np.insert(tvivoUnique, 0, 0.0)

        self.addSample(self.outputExperimentSingle, "levyAvg", tvitroUnique, tvivoUnique, "vivoAvg", "vitroAvg")

        self.outputExperimentSingle.write(self._getPath("experimentSingle.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        self._defineSourceRelation(self.inputInVitro.get(), self.outputExperiment)
        self._defineSourceRelation(self.inputInVivo.get(), self.outputExperiment)

        self._defineOutputs(outputExperimentSingle=self.outputExperimentSingle)
        self._defineSourceRelation(self.inputInVitro.get(), self.outputExperimentSingle)
        self._defineSourceRelation(self.inputInVivo.get(), self.outputExperimentSingle)

    def _validate(self):
        return []

    def _summary(self):
        retval = []
        return retval
