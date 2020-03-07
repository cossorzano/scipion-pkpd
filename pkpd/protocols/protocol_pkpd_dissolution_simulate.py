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
from pkpd.models.dissolution_models import *
from pkpd.protocols.protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from pkpd.pkpd_units import createUnit, strUnit

# tested in test_workflow_dissolution.py

class ProtPKPDDissolutionSimulate(ProtPKPD):
    """ Simulate a dissolution profile\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'simulate dissolution'

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('allowTlag', params.BooleanParam,
                      label="Allow lag", default=False,
                      help='Allow lag time before starting dissolution (t-tlag)')
        form.addParam('modelType', params.EnumParam, choices=["Zero order", "First order", "Fractional", "Weibull",
                                                              "Double Weibull", "Higuchi",
                                                              "Korsmeyer-Peppas", "Hixson-Crowell", "Hopfenberg",
                                                              "Hill",
                                                              "Makoid-Banakar",
                                                              "Splines2", "Splines3", "Splines4", "Spline5", "Splines6",
                                                              "Splines7", "Splines8", "Splines9", "Splines10"],
                      label="Dissolution model", default=3,
                      help='Zero order: Y=K*(t-[tlag])\n' \
                           'First order: Y=Ymax*(1-exp(-beta*(t-[tlag])))\n' \
                           'Fractional order: Y=Ymax-pow(Amax^alpha-alpha*beta*t,1/alpha))\n' \
                           'Weibull: Y=Ymax*(1-exp(-lambda*t^b))\n' \
                           'Double Weibull: Y=Ymax*(F1*(1-exp(-lambda1*t^b1))+(1-F1)*(1-exp(-lambda2*(t-tlag2)^b2)))\n' \
                           'Higuchi: Y=Ymax*t^0.5\n' \
                           'Korsmeyer-Peppas: Y=Ymax*t^m\n' \
                           'Hixson-Crowell: Y=Ymax*(1-(1-K*t)^3)\n'
                           'Hopfenberg: Y=Ymax*(1-(1-K*t)^m)\n'
                           'Hill: Y = Ymax*t^d/(g^d+t^d)\n'
                           'Makoid-Banakar: Ymax*(t/tmax)^b*exp(b*(1-t/tmax))\n'
                           'SplinesN: Y= Ymax*Bspline(t;N,tmax)\n')
        form.addParam('parameters', params.StringParam, label="Parameters", default="",
                      help='Parameter values for the simulation.\nExample: 2;5 is 2 for the first parameter, 5 for the second parameter\n'
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
        form.addParam('timeUnits', params.EnumParam, choices=["min", "h"],
                      label="Time units", default=0)
        form.addParam('resampleT', params.FloatParam, label='Simulation model time step=', default=1)
        form.addParam('resampleT0', params.FloatParam, label='Initial time for simulation', default=0)
        form.addParam('resampleTF', params.FloatParam, label='Final time for simulation', default=100)
        form.addParam('AUnits', params.StringParam, label="Dissolution units", default="%",
                      help="%, mg/dL, ...")

        form.addParam('noiseType', params.EnumParam, label="Type of noise to add",
                      choices=["None", "Additive", "Multiplicative"],
                      default=0, expertLevel=LEVEL_ADVANCED,
                      help='Additive: noise is normally distributed (mean=0 and standard deviation=sigma)\n' \
                           'Multiplicative: noise is normally distributed (mean=0 and standard deviation=sigma*X)\n')
        form.addParam('noiseSigma', params.FloatParam, label="Noise sigma",
                      default=0.0, expertLevel=LEVEL_ADVANCED, condition="noiseType>0",
                      help='See help of Type of noise to add\n')

    # --------------------------- STEPS functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulate')
        self._insertFunctionStep('createOutputStep')

    def addNoise(self, y):
        if self.noiseType.get() == 0:
            return y
        elif self.noiseType.get() == 1:
            return y + np.random.normal(0.0, self.noiseSigma.get(), y.shape)
        elif self.noiseType.get() == 2:
            return y * (1 + np.random.normal(0.0, self.noiseSigma.get(), y.shape))

    def runSimulate(self):
        tvar = PKPDVariable()
        tvar.varName = "t"
        tvar.varType = PKPDVariable.TYPE_NUMERIC
        tvar.role = PKPDVariable.ROLE_TIME
        if self.timeUnits.get() == 0:
            tvar.units = createUnit("min")
        elif self.timeUnits.get() == 1:
            tvar.units = createUnit("h")

        Avar = PKPDVariable()
        Avar.varName = "A"
        Avar.varType = PKPDVariable.TYPE_NUMERIC
        Avar.role = PKPDVariable.ROLE_MEASUREMENT
        if self.AUnits.get() != "%":
            Avar.units = createUnit(self.Aunits.get())
        else:
            Avar.units = createUnit("none")

        self.experimentSimulated = PKPDExperiment()
        self.experimentSimulated.variables["t"] = tvar
        self.experimentSimulated.variables["A"] = Avar
        self.experimentSimulated.general["title"] = "Simulated dissolution profile"
        self.experimentSimulated.general["comment"] = ""

        if self.modelType.get() == 0:
            self.model = Dissolution0()
        elif self.modelType.get() == 1:
            self.model = Dissolution1()
        elif self.modelType.get() == 2:
            self.model = DissolutionAlpha()
        elif self.modelType.get() == 3:
            self.model = DissolutionWeibull()
        elif self.modelType.get() == 4:
            self.model = DissolutionDoubleWeibull()
        elif self.modelType.get() == 5:
            self.model = DissolutionHiguchi()
        elif self.modelType.get() == 6:
            self.model = DissolutionKorsmeyer()
        elif self.modelType.get() == 7:
            self.model = DissolutionHixson()
        elif self.modelType.get() == 8:
            self.model = DissolutionHopfenberg()
        elif self.modelType.get() == 9:
            self.model = DissolutionHill()
        elif self.modelType.get() == 10:
            self.model = DissolutionMakoidBanakar()
        elif self.modelType.get() == 11:
            self.model = DissolutionSplines2()
        elif self.modelType.get() == 12:
            self.model = DissolutionSplines3()
        elif self.modelType.get() == 13:
            self.model = DissolutionSplines4()
        elif self.modelType.get() == 14:
            self.model = DissolutionSplines5()
        elif self.modelType.get() == 15:
            self.model = DissolutionSplines6()
        elif self.modelType.get() == 16:
            self.model = DissolutionSplines7()
        elif self.modelType.get() == 17:
            self.model = DissolutionSplines8()
        elif self.modelType.get() == 18:
            self.model = DissolutionSplines9()
        elif self.modelType.get() == 19:
            self.model = DissolutionSplines10()
        self.model.allowTlag = self.allowTlag.get()
        self.model.parameters = [float(x) for x in self.parameters.get().split(';')]

        newSample = PKPDSample()
        newSample.sampleName = "simulatedProfile"
        newSample.variableDictPtr = self.experimentSimulated.variables
        newSample.descriptors = {}
        newSample.addMeasurementPattern(["A"])

        t0 = self.resampleT0.get()
        tF = self.resampleTF.get()
        deltaT = self.resampleT.get()
        t = np.arange(t0, tF + deltaT, deltaT)
        y = self.model.forwardModel(self.model.parameters, t)
        newSample.addMeasurementColumn("t", t)
        newSample.addMeasurementColumn("A", y[0])

        self.experimentSimulated.samples[newSample.sampleName] = newSample
        self.experimentSimulated.write(self._getPath("experimentSimulated.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experimentSimulated)
