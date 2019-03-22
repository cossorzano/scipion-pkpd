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
import math
import numpy as np

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDExperiment, PKPDAllometricScale
from pkpd.pkpd_units import strUnit
from pyworkflow.protocol.constants import LEVEL_ADVANCED


class ProtPKPDAllometricScaling(ProtPKPD):
    """ Compute the allometric scaling between several spicies. This fits a model of the form Y=k*X^a where X can be
        any variable (although weight is normally used). The protocol can also leave out some parameters of the
        input model, and for these parameters, a simple average is calculated."""

    _label = 'allometric scaling'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExps', params.MultiPointerParam, label="Input experiments",
                      pointerClass='PKPDExperiment',
                      help='Select several results of the same ODE model. The different runs correspond to different '\
                           'spicies. Normally they are run over populations and only the first sample of each run is taken.')
        form.addParam('predictor', params.StringParam, label="Predictor variable (X)", default="weight",
                      help='Y is predicted as an exponential function of X, Y=f(X)')
        form.addParam('allometricParameters', params.StringParam, label="ODE parameters to allometric scale", default="",
                      help='List of variables separated by commas. The output for these parameters will be  fitted allometric scale model (Y=k*x^a).')
        form.addParam('avgParameters', params.StringParam, label="ODE parameters to average", default="",
                      help='List of variables separated by commas. The output for these parameters will simply be an average instead of an allometric model')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runFit', self.predictor.get(), self.avgParameters.get(), self.confidenceInterval.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runFit(self, predictor, avgParameters, confidenceInterval):
        self.scaleModel = PKPDAllometricScale()
        self.scaleModel.confidence = float(self.confidenceInterval.get())
        self.scaleModel.predictor = self.predictor.get()

        if self.allometricParameters.get().strip()!="":
            for token in self.allometricParameters.get().split(','):
                self.scaleModel.scaled_vars.append(token.strip())
        if self.avgParameters.get().strip()!="":
            for token in self.avgParameters.get().split(','):
                self.scaleModel.averaged_vars.append(token.strip())

        X=[]
        Y={}
        for varName in self.scaleModel.scaled_vars+self.scaleModel.averaged_vars:
            Y[varName]=[]

        for experimentPtr in self.inputExps:
            # Load experiment
            experiment = PKPDExperiment()
            experiment.load(experimentPtr.get().fnPKPD)

            # Get the first sample
            sample = experiment.samples.values()[0]

            # Extract X
            retval = sample.getDescriptorValue(self.predictor.get())
            if retval:
                X.append(float(retval))
            else:
                raise Exception("Cannot find %s in %s"%(self.predictor.get(),experiment.fnPKPD))

            # Extract Y
            for varName in self.scaleModel.scaled_vars+self.scaleModel.averaged_vars:
                retval = sample.getDescriptorValue(varName)
                if retval:
                    Y[varName].append(float(retval))
                else:
                    raise Exception("Cannot find %s in %s"%(varName,experiment.fnPKPD))

        self.scaleModel.predictorUnits = strUnit(experiment.variables[self.predictor.get()].units.unit)
        varList = []
        for varName in self.scaleModel.scaled_vars:
            varList.append((varName,strUnit(experiment.variables[varName].units.unit)))
        self.scaleModel.scaled_vars = copy.copy(varList)
        varList = []
        for varName in self.scaleModel.averaged_vars:
            varList.append((varName,strUnit(experiment.variables[varName].units.unit)))
        self.scaleModel.averaged_vars = copy.copy(varList)

        print("X: %s"%str(X))
        logx=np.log10(np.asarray(X))
        self.scaleModel.X = X
        self.scaleModel.Y = Y

        fhSummary=open(self._getPath("summary.txt"),"w")
        fhSummary.write("Predictor variable: %s\n"%self.predictor.get())
        fhSummary.write("Parameters to average: %s\n"%self.avgParameters.get())
        fhSummary.write("Confidence interval: %f\n"%self.confidenceInterval.get())
        for varName, varUnits in self.scaleModel.scaled_vars:
            print("%s: %s"%(varName,str(Y[varName])))
            y = np.asarray(Y[varName])
            logy = np.log10(y)

            A=np.ones((logx.size,2))
            A[:,0]=logx.T
            F=np.linalg.inv(np.dot(A.T,A))
            p = np.dot(F,np.dot(A.T,logy.T))
            logyPredicted=np.dot(A,p)
            e=logy-logyPredicted
            sigma2=np.var(e)
            V=sigma2*F

            # We manually implement the regression because np.polyfit has an unbiased estimator of V that
            # does not work well with 4 samples
            # p, V = np.polyfit(logx, logy, 1, cov=True)
            # e = logy - logyPredicted
            R2 = (1 - np.var(e) / np.var(logy))

            from scipy.stats import norm
            nstd = norm.ppf(1-(1-self.scaleModel.confidence/100)/2)
            perr = np.sqrt(np.diag(V))
            lowerBound = p - nstd * perr
            upperBound = p + nstd * perr

            # Back to natural units
            p[1]=math.pow(10.0,p[1])
            lowerBound[1]=math.pow(10.0,lowerBound[1])
            upperBound[1]=math.pow(10.0,upperBound[1])

            self.scaleModel.models[varName] = [p[1],p[0]]
            self.scaleModel.qualifiers[varName] = [R2,lowerBound[1],upperBound[1],lowerBound[0],upperBound[0]]
            self.doublePrint(fhSummary,
                "%s=(%f)*%s^(%f) R2=%f %f%% Confidence intervals (y=k*x^a): k=[%f,%f] a=[%f,%f]; %s" % \
                (varName, p[1], self.predictor, p[0], R2, self.confidenceInterval, lowerBound[1],upperBound[1],lowerBound[0],upperBound[0],
                 varUnits))

        for varName, varUnits in self.scaleModel.averaged_vars:
            print("%s: %s"%(varName,str(Y[varName])))
            y = np.asarray(Y[varName])
            self.scaleModel.models[varName] = [np.mean(y)]
            self.scaleModel.qualifiers[varName] = [np.std(y)]
            self.doublePrint(fhSummary,
                             "%s avg=%f std=%f; %s"%(varName,self.scaleModel.models[varName][0],self.scaleModel.qualifiers[varName][0],varUnits))

        fhSummary.close()
        self.scaleModel.write(self._getPath("allometricScale.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputScaleModel=self.scaleModel)
        for experimentPtr in self.inputExps:
            self._defineSourceRelation(experimentPtr.get(), self.scaleModel)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        self.addFileContentToMessage(msg,self._getPath("summary.txt"))
        return msg
