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

import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from .protocol_pkpd_fit_base import ProtPKPDFitBase
from pkpd.models.dissolution_models import *
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.pkpd_units import createUnit, strUnit

# Tested by test_workflow_dissolution
# Tested by test_workflow_levyplot
# Tested by test_workflow_deconvolution2


class ProtPKPDDissolutionFit(ProtPKPDFitBase):
    """ Fit a dissolution model. The observed measurement is modelled as Y=f(t).\n
Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'fit dissolution'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineParams1(form,"t","C")
        form.addParam('allowTlag', params.BooleanParam,
                      label="Allow lag", default=False,
                      help='Allow lag time before starting dissolution (t-tlag)')
        form.addParam('modelType', params.EnumParam, choices=["Zero order","First order","Fractional","Weibull",
                                                              "Double Weibull", "Higuchi",
                                                              "Korsmeyer-Peppas","Hixson-Crowell","Hopfenberg","Hill",
                                                              "Makoid-Banakar",
                                                              "Splines2", "Splines3", "Splines4", "Spline5", "Splines6",
                                                              "Splines7", "Splines8", "Splines9", "Splines10"],
                      label="Dissolution model", default=3,
                      help='Zero order: Y=K*(t-[tlag])\n'\
                           'First order: Y=Ymax*(1-exp(-beta*(t-[tlag])))\n'\
                           'Fractional order: Y=Ymax-pow(Amax^alpha-alpha*beta*t,1/alpha))\n'\
                           'Weibull: Y=Ymax*(1-exp(-lambda*t^b))\n'\
                           'Double Weibull: Y=Ymax*(F1*(1-exp(-lambda1*t^b1))+(1-F1)*(1-exp(-lambda2*(t-tlag2)^b2)))\n'\
                           'Higuchi: Y=Ymax*t^0.5\n'\
                           'Korsmeyer-Peppas: Y=Ymax*t^m\n'\
                           'Hixson-Crowell: Y=Ymax*(1-(1-K*t)^3)\n'
                           'Hopfenberg: Y=Ymax*(1-(1-K*t)^m)\n'
                           'Hill: Y = Ymax*t^d/(g^d+t^d)\n'
                           'Makoid-Banakar: Ymax*(t/tmax)^b*exp(b*(1-t/tmax))\n'
                           'SplinesN: Y= Ymax*Bspline(t;N,tmax)\n')
        form.addParam('fitType', params.EnumParam, choices=["Linear","Logarithmic","Relative"], label="Fit mode", default=0,
                      expertLevel=LEVEL_ADVANCED,
                      help='Linear: sum (Cobserved-Cpredicted)^2\nLogarithmic: sum(log10(Cobserved)-log10(Cpredicted))^2\n'\
                           "Relative: sum ((Cobserved-Cpredicted)/Cobserved)^2")
        form.addParam('bounds', params.StringParam, label="Bounds (optional)", default="",
                      help='Parameter values for the simulation.\nExample: (1,10);(0,0.05) is (1,10) for the first parameter, (0,0.05) for the second parameter\n'
                           'Zero order: [tlag];K\n'
                           'First order: [tlag];Ymax;beta\n'
                           'Fractional order: [tlag]; Ymax;beta;alpha\n'
                           'Weibull: [tlag]; Ymax;lambda;b\n'
                           'Double Weibull: [tlag]; Ymax; lambda1; b1; F1; tlag2; lambda2; b2\n'
                           'Higuchi: [tlag]; Ymax\n'
                           'Korsmeyer-Peppas: [tlag]; Ymax; m\n'
                           'Hixson-Crowell: [tlag]; Ymax; K\n'
                           'Hopfenberg: [tlag]; Ymax; K; m\n'
                           'Hill: [tlag]; Ymax; g; d\n'
                           'Makoid-Banakar: [tlag]; Ymax; Tmax; b\n'
                           'SplinesN: [tlag]; Ymax; tmax; c1; c2; ...; cN\n')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval=", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')
        form.addParam('resampleT',params.FloatParam, label='Simulation model time step=', default=-1,
                      help='If this value is greater than 0, then the fitted models will be sampled at this sampling period '
                           'and the created profiles will be collected in a new output experiment. '
                           'The time unit of the simulation is the same as the one of the predictor variable.')
        form.addParam('resampleT0',params.FloatParam, label='Initial time for simulation', default=0,
                      condition='resampleT>0',
                      help='The time unit of the simulation is the same as the one of the predictor variable.')
        form.addParam('resampleTF',params.FloatParam, label='Final time for simulation', default=0,
                      condition='resampleT>0', help='If set to 0, then the maximum time of the sample will be taken. '
                      'The time unit of the simulation is the same as the one of the predictor variable.')

    def getListOfFormDependencies(self):
        return [self.allowTlag.get(), self.modelType.get(), self.bounds.get(), self.confidenceInterval.get()]

    #--------------------------- STEPS functions --------------------------------------------
    def createModel(self):
        if self.modelType.get() == 0:
            return Dissolution0()
        elif self.modelType.get() == 1:
            return Dissolution1()
        elif self.modelType.get() == 2:
            return DissolutionAlpha()
        elif self.modelType.get() == 3:
            return DissolutionWeibull()
        elif self.modelType.get() == 4:
            return DissolutionDoubleWeibull()
        elif self.modelType.get() == 5:
            return DissolutionHiguchi()
        elif self.modelType.get() == 6:
            return DissolutionKorsmeyer()
        elif self.modelType.get() == 7:
            return DissolutionHixson()
        elif self.modelType.get() == 8:
            return DissolutionHopfenberg()
        elif self.modelType.get() == 9:
            return DissolutionHill()
        elif self.modelType.get() == 10:
            return DissolutionMakoidBanakar()
        elif self.modelType.get() == 11:
            return DissolutionSplines2()
        elif self.modelType.get() == 12:
            return DissolutionSplines3()
        elif self.modelType.get() == 13:
            return DissolutionSplines4()
        elif self.modelType.get() == 14:
            return DissolutionSplines5()
        elif self.modelType.get() == 15:
            return DissolutionSplines6()
        elif self.modelType.get() == 16:
            return DissolutionSplines7()
        elif self.modelType.get() == 17:
            return DissolutionSplines8()
        elif self.modelType.get() == 18:
            return DissolutionSplines9()
        elif self.modelType.get() == 19:
            return DissolutionSplines10()

    def setupFromFormParameters(self):
        self.model.allowTlag = self.allowTlag.get()

    def prepareForAnalysis(self):
        if self.resampleT.get()>0:
            self.experimentSimulated = PKPDExperiment()
            self.experimentSimulated.variables[self.fitting.predicted.varName]=self.fitting.predicted
            self.experimentSimulated.variables[self.fitting.predictor.varName]=self.fitting.predictor
            self.experimentSimulated.general["title"]="Simulated response from dissolution profiles"
            self.experimentSimulated.general["comment"]="Simulated response from dissolution profiles"

    def postSampleAnalysis(self, sampleName):
        self.experiment.addParameterToSample(sampleName, "tvitroMax",
                                             self.experiment.variables[self.varNameX].units.unit,
                                             "Maximum tvitro for which this fitting is valid", np.max(self.model.x))
        if self.resampleT.get()>0:
            newSample = PKPDSample()
            newSample.sampleName = sampleName+"_simulated"
            newSample.variableDictPtr = self.experimentSimulated.variables
            newSample.doseDictPtr = self.experimentSimulated.doses
            newSample.descriptors = {}
            newSample.addMeasurementPattern([self.fitting.predicted.varName])

            t0=self.resampleT0.get()
            tF=self.resampleTF.get()
            deltaT=self.resampleT.get()
            if tF==0:
                tF=np.max(self.model.x)
            t=np.arange(t0,tF+deltaT,deltaT)
            y=self.model.forwardModel(self.model.parameters,t)
            newSample.addMeasurementColumn(self.fitting.predictor.varName,t)
            newSample.addMeasurementColumn(self.fitting.predicted.varName,y[0])

            self.experimentSimulated.samples[sampleName] = newSample

    def postAnalysis(self):
        if self.resampleT.get()>0:
            self.experimentSimulated.write(self._getPath("experimentSimulated.pkpd"))

    def createOutputStep(self):
        ProtPKPDFitBase.createOutputStep(self)
        if self.resampleT.get()>0:
            self._defineOutputs(outputExperimentSimulated=self.experimentSimulated)
            self._defineSourceRelation(self.getInputExperiment(), self.experimentSimulated)
