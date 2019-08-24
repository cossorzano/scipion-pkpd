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

import sys

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDVariable, PKPDUnit
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pkpd.utils import parseRange
from pkpd.models.pd_models import *


class ProtPKPDSimulateGenericPD(ProtPKPD):
    """ Simulate a generic pharmacodynamic response Y=f(X).\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'simulate generic'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form, fullForm=True):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('predictor', params.StringParam, label="Predictor variable (X)", default="Cp",
                      help='The predictor variable is X in Y=f(X). X must be one of the variables of the experiment. '\
                           'It is normally a concentration.')
        form.addParam('predicted', params.StringParam, label="Predicted variable (Y)", default="E",
                      help='The predicted variable is Y in Y=f(X). Y must not exist in the experiment. It is normally the effect')
        form.addParam('predictedUnit', params.StringParam, label="Units for Y", default="none", expertLevel=LEVEL_ADVANCED,
                      help='Valid units are none, g, h, ...')
        form.addParam('predictedComment', params.StringParam, label="Comment for Y", default="Simulated pharmacodynamic feature",
                      help='This comment will appear in the output experiment')
        form.addParam('modelType', params.EnumParam, choices=["Linear (Polynomial 1)","Log-linear","Saturated","Sigmoid","Gompertz",
                                                              "Logistic1 ","Logistic 2","Logistic 3","Logistic 4",
                                                              "Richards","Morgan-Mercer-Flodin","Weibull"],
                      label="Generic model", default=0,
                      help='Linear: Y=e1*X+e0\nLog-linear: Y=m*log(X-X0)\n'\
                           'Saturated: Y=emin+emax*X/(ec50+X)\nSigmoid: Y=e0+(emax*(X**h)/((eC50**h)+(X**h)))\n'\
                           'Gompertz: Y=e0+a*exp(-exp(b-g*X))\nLogistic 1: Y=e0+(a/(1+exp(b-g*X)))\n'\
                           'Logistic 2: Y=e0+(1/(a+exp(b-g*X)))\nLogistic 3: Y=e0+(a/(1+b*exp(-g*X)))\n'\
                           'Logistic 4: Y=e0+(1/(a+b*exp(-g*X)))\nRichards: Y=e0+(a/((1+exp(b-g*X))^(1/d)))\n'\
                           'Morgan-Mercer-Flodin: Y=e0+((b*g+a*(X^d))/(g+(X^d)))\nWeibull: Y=a-b*exp(-g*(X^d))\n'\
                           'Polynomial n: en*X^n+...+e1*X+e0')

        form.addParam('paramValues', params.StringParam, label="Parameter values", default="",
                      help='Parameter values for the simulation.\nExample: 3.5;-1 is 3.5 for the first parameter, -1 for the second parameter\n'
                           'Linear: e1;e0\n'\
                           'Log-linear: m;X0\n'\
                           'Saturated: e0;emax;eC50\n'\
                           'Sigmoid: e0;emax;eC50;h\n'\
                           'Gompertz: e0;a;b;g\n'\
                           'Logistic 1: e0;a;b;g\n'\
                           'Logistic 2: e0;a;b;g\n'\
                           'Logistic 3: e0;a;b;g\n'\
                           'Logistic 4: e0;a;b;g\n'\
                           'Richards: e0;a;b;g;d\n'\
                           'Morgan-Mercer-Flodin: e0;b;g;a;d\n'\
                           'Weibull: a;b;g;d\n'\
                           'Polynomial n: en; ...; e1; e0\n')
        form.addParam('reportX', params.StringParam, label="Evaluate at X=", default="", expertLevel=LEVEL_ADVANCED,
                      help='Evaluate the model at these X values\nExample 1: [0,5,10,20,40,100]\nExample 2: 0:2:10, from 0 to 10 in steps of 2')
        form.addParam('noiseType', params.EnumParam, label="Type of noise to add", choices=["None","Additive","Multiplicative"],
                      default=0, expertLevel=LEVEL_ADVANCED,
                      help='Additive: noise is normally distributed (mean=0 and standard deviation=sigma)\n'\
                           'Multiplicative: noise is normally distributed (mean=0 and standard deviation=sigma*X)\n')
        form.addParam('noiseSigma', params.FloatParam, label="Noise sigma",
                      default=0.0, expertLevel=LEVEL_ADVANCED, condition="noiseType>0",
                      help='See help of Type of noise to add\n')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulate',self.inputExperiment.get().getObjId(), self.predictor.get(), \
                                 self.predicted.get(), self.modelType.get(), self.paramValues.get(), self.reportX.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addNoise(self,y):
        if self.noiseType.get()==0:
            return y
        elif self.noiseType.get()==1:
            return y+np.random.normal(0.0,self.noiseSigma.get(),y.shape)
        elif self.noiseType.get()==2:
            return y*(1+np.random.normal(0.0,self.noiseSigma.get(),y.shape))

    def runSimulate(self, objId, X, Y, modelType, paramValues, reportX):
        reportX = parseRange(self.reportX.get())
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        # Setup model
        self.printSection("Model setup")
        if self.modelType.get()==0:
            model = PDPolynomial1()
        elif self.modelType.get()==1:
            model = PDLogLinear()
        elif self.modelType.get()==2:
            model = PDSaturated()
        elif self.modelType.get()==3:
            model = PDSigmoid()
        elif self.modelType.get()==4:
            model = PDGompertz()
        elif self.modelType.get()==5:
            model = PDLogistic1()
        elif self.modelType.get()==6:
            model = PDLogistic2()
        elif self.modelType.get()==7:
            model = PDLogistic3()
        elif self.modelType.get()==8:
            model = PDLogistic4()
        elif self.modelType.get()==9:
            model = PDRichards()
        elif self.modelType.get()==10:
            model = PDMorgan()
        elif self.modelType.get()==11:
            model = PDWeibull()
        elif self.modelType.get()==12:
            model = PDPolynomial2()
        elif self.modelType.get()==13:
            model = PDPolynomial3()
        elif self.modelType.get()==14:
            model = PDPolynomial4()
        elif self.modelType.get()==15:
            model = PDPolynomial5()
        elif self.modelType.get()==16:
            model = PDPolynomial6()
        elif self.modelType.get()==17:
            model = PDPolynomial7()
        elif self.modelType.get()==18:
            model = PDPolynomial8()
        elif self.modelType.get()==19:
            model = PDPolynomial9()

        model.setExperiment(self.experiment)
        model.setXVar(self.predictor.get())
        model.printSetup()

        # Create list of parameters
        tokens=self.paramValues.get().split(';')
        if len(tokens)!=model.getNumberOfParameters():
                raise Exception("The list of parameter values has not the same number of parameters as the model")
        model.parameters=[]
        for token in tokens:
            try:
                model.parameters.append(float(token.strip()))
            except:
                raise Exception("Cannot convert %s to float"%token)
        print("Simulated model: %s"%model.getEquation())
        if self.noiseType.get()==1:
            print("Adding additive noise with sigma=%f"%self.noiseSigma.get())
        elif self.noiseType.get()==2:
            print("Adding multiplicative noise with sigma=%f*X"%self.noiseSigma.get())

        # Add new variable to experiment
        newVariable = PKPDVariable()
        newVariable.varName = self.predicted.get()
        newVariable.varType = PKPDVariable.TYPE_NUMERIC
        newVariable.displayString = "%f"
        newVariable.role = PKPDVariable.ROLE_MEASUREMENT
        newVariable.comment = self.predictedComment.get()
        newVariable.units = PKPDUnit(self.predictedUnit.get())
        self.experiment.variables[newVariable.varName] = newVariable

        for sampleName, sample in self.experiment.samples.iteritems():
            model.x = np.array(sample.getValues(self.predictor.get()),dtype=np.float)
            y = model.forwardModel(model.parameters,model.x)
            y = self.addNoise(y)
            sample.addMeasurementColumn(newVariable.varName,y)
            print("==========================================")
            sample._printToStream(sys.stdout)
            print("==========================================")
            sample._printMeasurements(sys.stdout)
            print(" ")
            if reportX!=None:
                print("Evaluation of the model at specified values")
                yReportX = model.forwardModel(model.parameters, reportX)
                yReportX = self.addNoise(yReportX)
                print("==========================================")
                print("X     Ypredicted     log10(Ypredicted)")
                print("==========================================")
                for n in range(0,reportX.shape[0]):
                    print("%f %f %f"%(reportX[n],yReportX[n],math.log10(yReportX[n])))
                print(' ')

    def createOutputStep(self):
        self.experiment.write(self._getPath("experiment.pkpd"))
        self._defineOutputs(outputFitting=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        modelTypeStr="unknown"
        if self.modelType.get()==0:
            modelTypeStr = "linear"
        elif self.modelType.get()==1:
           modelTypeStr = "Log-Linear"
        elif self.modelType.get()==2:
            modelTypeStr = "Saturated"
        elif self.modelType.get()==3:
            modelTypeStr = "Sigmoid"
        elif self.modelType.get()==4:
            modelTypeStr = "Gompertz"
        elif self.modelType.get()==5:
            modelTypeStr == "Logistic1"
        elif self.modelType.get() == 6:
            modelTypeStr == "Logistic 2"
        elif self.modelType.get() == 7:
            modelTypeStr == "Logistic 3"
        elif self.modelType.get() == 8:
            modelTypeStr == "Logistic 4"
        elif self.modelType.get() == 9:
            modelTypeStr == "Richards"
        elif self.modelType.get() == 10:
            modelTypeStr == "Morgan-Mercer-Flodin"
        elif self.modelType.get() == 11:
            modelTypeStr == "Weibull"
        elif self.modelType.get() == 12:
            modelTypeStr == "Polynomial 2"
        elif self.modelType.get() == 13:
            modelTypeStr == "Polynomial 3"
        elif self.modelType.get() == 14:
            modelTypeStr == "Polynomial 4"
        elif self.modelType.get() == 15:
            modelTypeStr == "Polynomial 5"
        elif self.modelType.get() == 16:
            modelTypeStr == "Polynomial 6"
        elif self.modelType.get() == 17:
            modelTypeStr == "Polynomial 7"
        elif self.modelType.get() == 18:
            modelTypeStr == "Polynomial 8"
        elif self.modelType.get() == 19:
            modelTypeStr == "Polynomial 9"

        msg.append("Variable %s added to experiment by simulation of a %s model"%(self.predicted.get(),modelTypeStr))
        msg.append("Parameter values: %s"%self.paramValues)
        return msg

    def _validate(self):
        msg=[]
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        if self.predicted.get() in experiment.variables:
            msg.append("The experiment already has a column called %s"%self.predicted.get())
        units = PKPDUnit()
        if units._fromString(self.predictedUnit.get()) is None:
            msg.append("Predicted unit is not valid")
        return msg

    def filterVarForWizard(self, v):
        """ Define the type of variables required (used in wizard). """
        return v.isMeasurement() and v.isNumeric()
