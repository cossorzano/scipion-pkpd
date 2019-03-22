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
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDDoseResponse
from pkpd.utils import parseRange
from pkpd.models.pd_models import *
from numpy.random import uniform


class ProtPKPDSimulateDoseEscalation(ProtPKPD):
    """ Simulate a dose escalation\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'simulate dose escalation'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form, fullForm=True):
        form.addSection('Input')
        form.addParam('modelType', params.EnumParam, choices=["OQuigley0","OQuigley1", "OQuigley2", "Sigmoid", "Gompertz", "Logistic", "Richards"],
                      label="Response model", default=0,
                      help='OQuigley0: Y=((tanh(X)+1)/2)^a. Order: a\n'\
                           'OQuigley1: Y=((tanh(X-X0)+1)/2)^a. Order: X0;a\n'\
                           'OQuigley2: Y=exp(g*(X-X0))/(1+exp(g*(X-X0)). Order: X0;g\n'\
                           'Sigmoid: Y=((X**h)/((X50**h)+(X**h))). Order X50;h\n'\
                           'Gompertz: Y=exp(-exp(g*(X-X0))). Order: X0;g\n'\
                           'Logistic: Y=1/(1+exp(g*(X-X0))). Order: X0;g\n'\
                           'Richards: Y=1/((1+exp(g*(X-X0)))^(1/d)). Order: X0;g;d\n')

        form.addParam('paramValues', params.StringParam, label="Parameter values", default="",
                      help='Parameter values for the simulation.\nExample: 3.5;-1 is 3.5 for the first parameter, -1 for the second parameter\n'
                           'OQuigley0: a\n'\
                           'OQuigley1: X0;a\n'\
                           'OQuigley2: X0;g\n'\
                           'Sigmoid: X50;h\n'\
                           'Gompertz: X0;g\n'\
                           'Logistic: X0;g\n'\
                           'Richards: X0;g;d\n')

        form.addParam('reportX', params.StringParam, label="Evaluate at X=", default="",
                      help='Evaluate the model at these X values\nExample 1: [0,5,10,20,40,100]\nExample 2: 0:2:10, from 0 to 10 in steps of 2')
        form.addParam('doLog', params.BooleanParam, label="Take log10 in the dose", default=False,
                      help='In the formulas, X is substituted by log10(X)')
        form.addParam('prot1ptr', params.PointerParam, label="Previous dose response",
                      pointerClass='ProtPKPDDoseEscalation,ProtPKPDSimulateDoseEscalation',allowsNull=True,
                      help='Optional. You may link this simulation to a previous dose escalation analysis, so that you may know the sequence')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulate', self.modelType.get(), self.paramValues.get(), self.reportX.get())

    #--------------------------- STEPS functions --------------------------------------------
    def runSimulate(self, modelType, paramValues, reportX):
        reportXorig = parseRange(self.reportX.get())
        reportX = parseRange(self.reportX.get())
        if self.doLog:
            reportX=np.log10(np.clip(reportX,0.0,None))
            print("Taking logarithms: X is really log10(X)")

        # Setup model
        self.printSection("Model setup")
        if self.modelType.get()==0:
            model = PDOQuigley0()
        elif self.modelType.get()==1:
            model = PDOQuigley1()
        elif self.modelType.get()==2:
            model = PDOQuigley2()
        elif self.modelType.get()==3:
            model = PDSigmoid()
        elif self.modelType.get()==4:
            model = PDGompertz()
        elif self.modelType.get()==5:
            model = PDLogistic1()
        elif self.modelType.get()==6:
            model = PDRichards()

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

        if reportX!=None:
            dr = PKPDDoseResponse()
            if self.prot1ptr.get() is not None:
                dr.read(self.prot1ptr.get()._getPath("dose_response.txt"), True)
            print("Evaluation of the model at specified values")
            yReportX = model.forwardModel(model.parameters, [reportX])
            yReportX = yReportX[0] # From [array(...)] to array(...)
            print("==========================================")
            print("X     Xused     Ypredicted     RandomResponse")
            print("==========================================")
            for n in range(0,reportX.shape[0]):
                response=uniform()<yReportX[n]
                dr.appendResponse(reportXorig[n],response)
                print("%f %f %f %s"%(reportXorig[n],reportX[n],yReportX[n],response))
            print(' ')
            dr.write(self._getPath("dose_response.txt"))


    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        modelTypeStr="unknown"
        if self.modelType.get()==0:
            modelTypeStr = "OQuigley 0"
        elif self.modelType.get()==1:
            modelTypeStr = "OQuigley 1"
        elif self.modelType.get()==2:
           modelTypeStr = "OQuigley 2"
        elif self.modelType.get()==3:
            modelTypeStr = "Sigmoid"
        elif self.modelType.get()==4:
            modelTypeStr = "Gompertz"
        elif self.modelType.get()==5:
            modelTypeStr == "Logistic"
        elif self.modelType.get()==6:
            modelTypeStr == "Richards"
        msg.append("Model type: %s"%modelTypeStr)
        msg.append("Parameter values: %s"%self.paramValues)
        return msg
