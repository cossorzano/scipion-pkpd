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
from pkpd.models.pd_models import *

# TESTED in test_workflow_gabrielsson_pk11.py
# TESTED in test_workflow_gabrielsson_pk15.py
# TESTED in test_workflow_gabrielsson_pd03.py
# TESTED in test_workflow_gabrielsson_pd11.py

class ProtPKPDGenericFit(ProtPKPDFitBase):
    """ Fit a generic model. The observed measurement is modelled as Y=f(X).\n
Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'fit pd generic'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        self._defineParams1(form,"Cp","E")
        form.addParam('modelType', params.EnumParam, choices=["Linear (polynomial 1)","Log-linear","Saturated","Sigmoid","Gompertz",
                                                              "Logistic1 ","Logistic 2","Logistic 3","Logistic 4",
                                                              "Richards","Morgan-Mercer-Flodin","Weibull","Hill",
                                                              "Polynomial 2", "Polynomial 3", "Polynomial 4",
                                                              "Polynomial 5", "Polynomial 6", "Polynomial 7",
                                                              "Polynomial 8", "Polynomial 9"],
                      label="Generic model", default=0,
                      help='Linear: Y=e1*X+e0\nLog-linear: Y=m*log(X-X0)\nSaturated: Y=e0+(emax*X/(eC50+X))\n'\
                            'Sigmoid: Y=e0+(emax*(X**h)/((eC50**h)+(X**h)))\nGompertz: Y=e0+a*exp(-exp(b-g*X))\n'\
                            'Logistic 1: Y=e0+(a/(1+exp(b-g*X)))\nLogistic 2: Y=e0+(1/(a+exp(b-g*X)))\nLogistic 3: Y=e0+(a/(1+b*exp(-g*X)))\n'\
                            'Logistic 4: Y=e0+(1/(a+b*exp(-g*X)))\nRichards: Y=e0+(a/((1+exp(b-g*X))^(1/d)))\n'\
                            'Morgan-Mercer-Flodin: Y=e0+((b*g+a*(X^d))/(g+(X^d)))\nWeibull: Y=a-b*exp(-g*(X^d))\n'
                            'Hill: Y=e0+b*X^d/(g^d+X^d)\nPolynomial 2: e2*X^2+e1*X+e0\nPolynomial N: eN*X^N+...+e1*X+e0')
        form.addParam('fitType', params.EnumParam, choices=["Linear","Logarithmic","Relative"], label="Fit mode", default=1,
                      help='Linear: sum (Cobserved-Cpredicted)^2\nLogarithmic: sum(log10(Cobserved)-log10(Cpredicted))^2\n'\
                           "Relative: sum ((Cobserved-Cpredicted)/Cobserved)^2")
        form.addParam('bounds', params.StringParam, label="Bounds (optional)", default="",
                      help='Parameter values for the simulation.\nExample: (1,10);(0,0.05) is (1,10) for the first parameter, (0,0.05) for the second parameter\n'
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
                           'Hill: e0;b;g;d\n'\
                           'Polynomial N: eN;...;e1;e0\n')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval=", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')
        form.addParam('reportX', params.StringParam, label="Evaluate at X=", default="", expertLevel=LEVEL_ADVANCED,
                      help='Evaluate the model at these X values\nExample 1: [0,5,10,20,40,100]\nExample 2: 0:2:10, from 0 to 10 in steps of 2')

    def getListOfFormDependencies(self):
        return [self.modelType.get(), self.fitType.get(), self.bounds.get(), self.confidenceInterval.get(),
                self.reportX.get()]

    #--------------------------- STEPS functions --------------------------------------------
    def createModel(self):
        if self.modelType.get() == 0:
            return PDPolynomial1()
        elif self.modelType.get() == 1:
            return PDLogLinear()
        elif self.modelType.get() == 2:
            return PDSaturated()
        elif self.modelType.get() == 3:
            return PDSigmoid()
        elif self.modelType.get() == 4:
            return PDGompertz()
        elif self.modelType.get() == 5:
            return PDLogistic1()
        elif self.modelType.get() == 6:
            return PDLogistic2()
        elif self.modelType.get() == 7:
            return PDLogistic3()
        elif self.modelType.get() == 8:
            return PDLogistic4()
        elif self.modelType.get() == 9:
            return PDRichards()
        elif self.modelType.get() == 10:
            return PDMorgan()
        elif self.modelType.get() == 11:
            return PDWeibull()
        elif self.modelType.get() == 12:
            return PDHill()
        elif self.modelType.get() == 13:
            return PDPolynomial2()
        elif self.modelType.get() == 14:
            return PDPolynomial3()
        elif self.modelType.get() == 15:
            return PDPolynomial4()
        elif self.modelType.get() == 16:
            return PDPolynomial5()
        elif self.modelType.get() == 17:
            return PDPolynomial6()
        elif self.modelType.get() == 18:
            return PDPolynomial7()
        elif self.modelType.get() == 19:
            return PDPolynomial8()
        elif self.modelType.get() == 20:
            return PDPolynomial9()
