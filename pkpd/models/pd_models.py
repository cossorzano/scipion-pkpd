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
"""
PD models
"""

import numpy as np

from pkpd.objects import PKPDModel
from pkpd.pkpd_units import inverseUnits, divideUnits, unitFromString, PKPDUnit

import math

class PDModel(PKPDModel):
    pass

class PDGenericModel(PDModel):
    pass


class PDPolynomial(PDGenericModel):
    def __init__(self):
        self.N = 0

    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        self.yPredicted = np.polyval(parameters,xToUse)
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Polynomial of degree %d (%s)"%(self.N, self.__class__.__name__)

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            p = np.polyfit(xToUse,yToUse,self.N)
            print("First estimate of polynomial: ")
            print(self.getEquation(p))

            self.bounds = []
            for i in range(self.N+1):
                self.bounds.append((-10*np.abs(p[i]),10*np.abs(p[i])))

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Bounds: "+str(getattr(self,"bounds",[])))

    def getModelEquation(self):
        toPrint=""
        for i in range(0,self.N+1):
            if i>0:
                toPrint+="+"
            toPrint+="e%d*X^%d"%(self.N-i,self.N-i)
        return toPrint

    def getEquation(self, p=None):
        if p is None:
            p=self.parameters
        toPrint=""
        for i in range(0,self.N+1):
            if i>0:
                toPrint+="+"
            toPrint+="(%f)*X^%d"%(p[i],self.N-i)
        return toPrint

    def getParameterNames(self):
        retval = []
        for i in range(0,self.N+1):
            retval+=['e%d'%(self.N-i)]
        return retval

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=%s'%self.getModelEquation()]*self.getNumberOfParameters()

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        sunits = divideUnits(yunits,xunits)
        self.parameterUnits=[]
        for i in range(0,self.N+1):
            deg=self.N-i
            if deg==0:
                self.parameterUnits.append(yunits)
            elif deg==1:
                self.parameterUnits.append(sunits)
            else:
                self.parameterUnits.append(PKPDUnit.UNIT_NONE)
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        for i in range(0,self.N+1):
            retval.append(lowerBound[i]>0 or upperBound[i]<0)
        return retval

    def areParametersValid(self, p):
        return True

class PDPolynomial1(PDPolynomial):
    def __init__(self):
        self.N=1

class PDPolynomial2(PDPolynomial):
    def __init__(self):
        self.N=2

class PDPolynomial3(PDPolynomial):
    def __init__(self):
        self.N=3

class PDPolynomial4(PDPolynomial):
    def __init__(self):
        self.N=4

class PDPolynomial5(PDPolynomial):
    def __init__(self):
        self.N=5

class PDPolynomial6(PDPolynomial):
    def __init__(self):
        self.N=6

class PDPolynomial7(PDPolynomial):
    def __init__(self):
        self.N=7

class PDPolynomial8(PDPolynomial):
    def __init__(self):
        self.N=8

class PDPolynomial9(PDPolynomial):
    def __init__(self):
        self.N=9


class PDLogLinear(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])

        m = parameters[0]
        C0 = parameters[1]
        self.yPredicted = [m*math.log(xi - C0) if np.isfinite(xi) and xi>C0 else float("inf") for xi in xToUse]
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Log-Linear (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Cmin=np.min(xToUse)
            xprime = xToUse - 0.9*Cmin
            idx = np.where(xprime>0)[0]
            p = np.polyfit(np.log(xprime[idx]), yToUse[idx], 1)
            C0=0.9*Cmin
            m=p[0]
            print("First estimate of log-linear term: ")
            print("Y=(%f)*log(X - (%f))" % (m, C0))

            self.bounds = []
            self.bounds.append((0.1 * m, 10 * m))
            self.bounds.append((-9 * C0, 11 * C0))

    def printSetup(self):
        print ("Model: %s " %self.getModelEquation())
        print ("Bounds:  " + str(self.bounds))

    def getModelEquation(self):
        return "Y = m*log(X - C0)"

    def getEquation(self):
        toPrint = "Y=(%f)*log(X - (%f))" % (self.parameters[0], self.parameters[1])
        return toPrint

    def getParameterNames(self):
        return ['m','C0']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y = m*log(X - C0)'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, xunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        return retval

    def areParametersValid(self, p):
        return True


class PDSaturated(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])

        e0 = parameters[0]
        emax = parameters[1]
        eC50 = parameters[2]

        self.yPredicted = e0 + (emax*xToUse / (eC50 + xToUse))
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Saturated (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            e0 = np.min(yToUse)
            emax = np.max(yToUse-e0)
            eC50 = 0.5*(np.max(xToUse)+np.min(xToUse))
            print("First estimate of saturated term: ")
            print("Y=(%f) + ( (%f)*X / ((%f) + X) )" % (e0, emax, eC50))

            self.bounds = []
            self.bounds.append((0.1 * e0, 10 * e0))
            self.bounds.append((min(0.1*emax, 10*emax),max(0.1*emax,10*emax)))
            self.bounds.append((0.1 * eC50, 10 * eC50))

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0 + (emax*X / (eC50 + X))"

    def getEquation(self):
        toPrint = "Y = (%f) + ( (%f)*X / ((%f) + X) )"%(self.parameters[0], self.parameters[1], self.parameters[2])
        return toPrint

    def getParameterNames(self):
        return ['e0', 'emax', 'eC50']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y = e0 + (emax*X / (eC50 + X))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, yunits, xunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(True)
        retval.append(True)
        return retval

    def areParametersValid(self, p):
        return True


class PDSigmoid(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        e0 = parameters[0]
        emax = parameters[1]
        eC50 = parameters[2]
        h = parameters[3]
        eC50prime = eC50**h
        xprime = xToUse**h
        self.yPredicted = e0  + ( (emax*(xprime)) / ( (eC50prime) + (xprime)))
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Sigmoid (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            e0 = np.min(yToUse)
            emax = np.max(yToUse - e0)
            eC50 = 0.5 * (np.max(xToUse) + np.min(xToUse))
            h = 1 #por defecto
            print("First estimate of saturated term: ")
            print("Y = (%f) + ( (%f)*(X**(%f)) / ( ((%f)**(%f)) + (X**(%f))))" % (e0, emax, h, eC50, h, h))

            self.bounds = []
            self.bounds.append((0.1 * e0, 10 * e0))
            self.bounds.append((min(0.1*emax, 10*emax),max(0.1*emax,10*emax)))
            self.bounds.append((0.1 * eC50, 10 * eC50))
            self.bounds.append((1 * h, 10 * h))

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0 + ( emax*(X**h) / ( (eC50**h) + (X**h)))"

    def getEquation(self):
        toPrint = "Y = (%f) + ( (%f)*(X**(%f)) / ( ((%f)**(%f)) + (X**(%f))))"%(self.parameters[0], self.parameters[1],
                                                                                self.parameters[3], self.parameters[2],
                                                                                self.parameters[3], self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['e0', 'emax', 'eC50', 'h']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form  Y = e0 + ( emax*(X**h) / ( (eC50**h) + (X**h)) )'] * self.getNumberOfParameters()

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, yunits, xunits, PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True


class PDGompertz(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])

        e0 = parameters[0]
        a = parameters[1]
        b = parameters[2]
        g = parameters[3]

        d = np.exp(b - (g*xToUse))
        self.yPredicted = e0 + a * np.exp(-d)
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted


    def getDescription(self):
        return "Gompertz (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds==None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            e0 = np.min(yToUse)

            emax = np.max(yToUse - e0 )
            a=emax

            emin = 0.9*e0+0.1*(emax+e0)
            b = math.log(-math.log((emin-e0)/emax))
            g = (b + 5) / np.max(xToUse)


            print("First estimate of Gompertz term: ")
            print("Y = (%f) + (%f) *exp(-exp((%f) - (%f) * X))" %(e0,a,b,g))

            self.bounds = []
            self.bounds.append((min(0.1*e0, 10*e0),max(0.1*e0,10*e0)))
            self.bounds.append((min(0.1*a, 10*a),max(0.1*a,10*a)))
            self.bounds.append((min(0.1*b, 10*b),max(0.1*b,10*b)))
            self.bounds.append((min(0.1*g, 10*g),max(0.1*g,10*g)))

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0 + a*exp(-exp(b-g*X))"

    def getEquation(self):
        toPrint = "Y = (%f) + (%f)*exp(-exp((%f)-(%f)*X))"%(self.parameters[0], self.parameters[1], self.parameters[2], self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['e0','a', 'b', 'g']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=e0+a*exp(-exp(b-g*X))']*self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, yunits, PKPDUnit.UNIT_NONE, inverseUnits(xunits)]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True


class PDLogistic1(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        e0 = parameters[0]
        a = parameters[1]
        b = parameters[2]
        g = parameters[3]

        d = np.exp(b - (g * xToUse))

        self.yPredicted = e0 + (a / (1 + d))
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Logistic1 (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            e0 = np.min(yToUse)
            emax = np.max(yToUse - e0)
            a = emax

            emin = 0.9 * e0 + 0.1 * (emax + e0)

            b = math.log(emax/(emin-e0) - 1)
            g = (b + 5)/np.max(xToUse)

            print("First estimate of Logistic 1 term: ")
            print("Y = (%f) + ( (%f) / (1+exp((%f) - (%f) * X)) )" % (e0,a,b,g))

            self.bounds = []

            self.bounds.append((min(0.1 * e0, 10 * e0), max(0.1 * e0, 10 * e0)))
            self.bounds.append((min(0.1 * a, 10 * a), max(0.1 * a, 10 * a)))
            self.bounds.append((min(0.1 * b, 10 * b), max(0.1 * b, 10 * b)))
            self.bounds.append((min(0.1 * g, 10 * g), max(0.1 * g, 10 * g)))


    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0 + a/(1+exp(b-g*X))"

    def getEquation(self):
        toPrint = "Y = (%f) + ( (%f)/(1+exp((%f)-(%f)*X)) )" % (self.parameters[0], self.parameters[1],
                                                            self.parameters[2], self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['e0', 'a', 'b', 'g']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=e0 + a/(1+exp(b-g*X))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, yunits, PKPDUnit.UNIT_NONE, inverseUnits(xunits)]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True


class PDLogistic2(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        e0 = parameters[0]
        a = parameters[1]
        b = parameters[2]
        g = parameters[3]

        d = np.exp(b - (g * xToUse))

        self.yPredicted =  e0 + (1 / (a + d))
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Logistic2 (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            e0 = np.min(yToUse)
            emax = np.max(yToUse - e0)

            a = (1 / emax)

            emin = 0.9*e0 + 0.1*(emax + e0)

            b = math.log((1/(emin-e0)) - a)
            g = (b + 5) / np.max(xToUse)

            print("First estimate of Logistic 2 term: ")
            print("Y = (%f) + ( 1 / ((%f) + exp((%f) - (%f) * X)) )" % (e0, a, b, g))


            self.bounds = []
            self.bounds.append((min(0.1 * e0, 10 * e0), max(0.1 * e0, 10 * e0)))
            self.bounds.append((min(0.01 * a, 10 * a), max(0.01 * a, 10 * a)))
            self.bounds.append((min(0.1 * b, 10 * b), max(0.1 * b, 10 * b)))
            self.bounds.append((min(0.1 * g, 10 * g), max(0.1 * g, 10 * g)))

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0+(1/(a+exp(b-g*X)))"

    def getEquation(self):
        toPrint = "Y = (%f)+(1/((%f)+exp((%f)-(%f)*X)))" % (self.parameters[0], self.parameters[1],
                                                          self.parameters[2], self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['e0', 'a', 'b', 'g']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=e0+(1/(a+exp(b-g*X)))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, yunits, PKPDUnit.UNIT_NONE, inverseUnits(xunits)]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True


class PDLogistic3(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])

        e0 = parameters[0]
        a = parameters[1]
        b = parameters[2]
        g = parameters[3]

        d = np.exp(-(g * xToUse))

        self.yPredicted = e0 + ( a / (1 + b*d) )
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Logistic3 (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            e0 = np.min(yToUse)
            emax = np.max(yToUse - e0)
            a = emax

            emin = 0.9*e0 + 0.1*(emax + e0)

            b = (emax/(emin - e0)) - 1
            g = 5 / np.max(xToUse)

            print("First estimate of Logistic 3 term: ")
            print("Y = (%f) + ( (%f) / (1 + (%f) * exp(-(%f) * X)) )" % (e0, a, b, g))

            self.bounds = []
            self.bounds.append((min(0.1 * e0, 10 * e0), max(0.1 * e0, 10 * e0)))
            self.bounds.append((min(0.1 * a, 10 * a), max(0.1 * a, 10 * a)))
            self.bounds.append((min(0.1 * b, 10 * b), max(0.1 * b, 10 * b)))
            self.bounds.append((min(0.1 * g, 10 * g), max(0.1 * g, 10 * g)))


    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0 + ( a/(1+b*exp(-g*X)) )"

    def getEquation(self):
        toPrint = "Y = (%f) + ( (%f)/(1+(%f)*exp(-(%f)*X)) )" % (self.parameters[0], self.parameters[1],
                                                                 self.parameters[2], self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['e0', 'a', 'b', 'g']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=e0+(a/(1+b*exp(-g*X)))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, yunits, PKPDUnit.UNIT_NONE, inverseUnits(xunits)]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True


class PDLogistic4(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        e0 = parameters[0]
        a = parameters[1]
        b = parameters[2]
        g = parameters[3]

        d = np.exp(-(g * xToUse))

        self.yPredicted = e0 + ( 1 / (a + b*d) )
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Logistic4 (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            e0 = np.min(yToUse)
            emax = np.max(yToUse - e0)
            a = 1 / emax

            emin = 0.9 * e0 + 0.1 * (emax + e0)
            b = (1 / (emin-e0)) - emax
            g = 5 / np.max(xToUse)

            print("First estimate of Logistic 4 term: ")
            print("Y = (%f) + ( 1 / ((%f)+ (%f)*exp(-(%f)*X)) )" % (e0, a, b, g))

            self.bounds = []
            self.bounds.append((min(0.1 * e0, 10 * e0), max(0.1 * e0, 10 * e0)))
            self.bounds.append((min(0.1 * a, 10 * a), max(0.1 * a, 10 * a)))
            self.bounds.append((min(0.1 * b, 10 * b), max(0.1 * b, 10 * b)))
            self.bounds.append((min(0.1 * g, 10 * g), max(0.1 * g, 10 * g)))


    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0 + ( 1/(a+b*exp(-g*X)) )"

    def getEquation(self):
        toPrint = "Y = (%f) + ( 1/((%f)+(%f)*exp(-(%f)*X)) )" % (self.parameters[0], self.parameters[1],
                                                                 self.parameters[2], self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['e0', 'a', 'b', 'g']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=e0+(1/(a+b*exp(-g*X)))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, yunits, PKPDUnit.UNIT_NONE, inverseUnits(xunits)]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True


class PDRichards(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        e0 = parameters[0]
        a = parameters[1]
        b = parameters[2]
        g = parameters[3]
        d = parameters[4]

        p = np.exp(b-(g * xToUse))

        self.yPredicted = e0 + ( a / ((1+p)**(1/d)) )
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Richards (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            d = 1

            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            e0 = np.min(yToUse)
            emax = np.max(yToUse - e0)
            a = emax

            emin = 0.9 * e0 + 0.1 * (emax + e0)

            b = math.log(emax/(emin-e0) -1)
            g = (b + 5) / np.max(xToUse)

            print("First estimate of Richards term: ")
            print("Y = (%f) + ( (%f)/ ((1 + exp((%f) - (%f)*X))^(1/(%f))) ) " % (e0, a, b, g, d))

            self.bounds = []
            self.bounds.append((min(0.1 * e0, 10 * e0), max(0.1 * e0, 10 * e0)))
            self.bounds.append((min(0.1 * a, 10 * a), max(0.1 * a, 10 * a)))
            self.bounds.append((min(0.1 * b, 10 * b), max(0.1 * b, 10 * b)))
            self.bounds.append((min(0.1 * g, 10 * g), max(0.1 * g, 10 * g)))
            self.bounds.append((min(0.1 * d, 10 * d), max(0.1 * d, 10 * d)))


    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0 + ( a/ ((1 + exp(b - g*X))^(1/d)) )"

    def getEquation(self):
        toPrint = "Y = (%f) + ( (%f)/ ((1 + exp((%f) - (%f)*X))^(1/(%f))) )" % (self.parameters[0], self.parameters[1],
                                                                                self.parameters[2], self.parameters[3],
                                                                                self.parameters[4])
        return toPrint

    def getParameterNames(self):
        return ['e0', 'a', 'b', 'g', 'd']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y = e0 + ( a/ ((1 + exp(b - g*X))^(1/d)) )'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, yunits, PKPDUnit.UNIT_NONE, inverseUnits(xunits), PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        retval.append(lowerBound[4] > 0 or upperBound[4] < 0)
        return retval

    def areParametersValid(self, p):
        return True


class PDMorgan(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        e0 = parameters[0]
        b = parameters[1]
        g = parameters[2]
        a = parameters[3]
        d = parameters[4]

        xprime = xToUse**d

        self.yPredicted = e0 + ( ((b*g) + (a*xprime)) / (g + xprime) )
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Morgan-Mercer-Flodin (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            d = 1
            g = 1

            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            e0 = np.min(yToUse)
            emax = np.max(yToUse - e0)
            a = emax

            emin = 0.9 * e0 + 0.1 * (emax + e0)

            b = emin - e0

            print("First estimate of Morgan term: ")
            print("Y = Y = (%f) + ( ((%f)*(%f) + (%f)*(X^(%f))) / ((%f) + (X^(%f))) ) " % (e0, b, g, a, d, g, d))

            self.bounds = []
            self.bounds.append((min(0.1 * e0, 10 * e0), max(0.1 * e0, 10 * e0)))
            self.bounds.append((min(0.1 * b, 10 * b), max(0.1 * b, 10 * b)))
            self.bounds.append((min(0.1 * g, 10 * g), max(0.1 * g, 10 * g)))
            self.bounds.append((min(0.1 * a, 10 * a), max(0.1 * a, 10 * a)))
            self.bounds.append((min(0.1 * d, 10 * d), max(0.1 * d, 10 * d)))

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0 + ( (b*g + a*(X^d)) / (g + (X^d)) )"

    def getEquation(self):
        toPrint = "Y = (%f) + ( ((%f)*(%f) + (%f)*(X^(%f))) / ((%f) + (X^(%f))) )" % (self.parameters[0],self.parameters[1],
                                                                                      self.parameters[2],self.parameters[3],
                                                                                      self.parameters[4], self.parameters[2],
                                                                                      self.parameters[4])
        return toPrint

    def getParameterNames(self):
        return ['e0', 'b', 'g', 'a', 'd']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=e0+((b*g+a*(X^d))/(g+(X^d)))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, yunits, PKPDUnit.UNIT_NONE, inverseUnits(xunits), PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        retval.append(lowerBound[4] > 0 or upperBound[4] < 0)
        return retval

    def areParametersValid(self, p):
        return True


class PDWeibull(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        a = parameters[0]
        b = parameters[1]
        g = parameters[2]
        d = parameters[3]

        xprime = xToUse**d

        self.yPredicted = a - (b*(np.exp(- g * xprime)))
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Weibull (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            d = 1

            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            emin = np.min(yToUse)
            emax = np.max(yToUse)
            a = emax

            b = emax - emin
            g = 5 / np.max(xToUse)

            print("First estimate of Weibull term: ")
            print("Y  = (%f) - (%f)*exp(-(%f)*(X^(%f))) " % (a, b, g, d))

            self.bounds = []
            self.bounds.append((min(0.1 * a, 10 * a), max(0.1 * a, 10 * a)))
            self.bounds.append((min(0.1 * b, 10 * b), max(0.1 * b, 10 * b)))
            self.bounds.append((min(0.1 * g, 10 * g), max(0.1 * g, 10 * g)))
            self.bounds.append((min(0.1 * d, 10 * d), max(0.1 * d, 10 * d)))

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = a - b*exp(-g*(X^d))"

    def getEquation(self):
        toPrint = "Y  = (%f) - (%f)*exp(-(%f)*(X^(%f)))" % (self.parameters[0], self.parameters[1], self.parameters[2],
                                                            self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['a', 'b', 'g', 'd']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=a-b*exp(-g*(X^d))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, PKPDUnit.UNIT_NONE, inverseUnits(xunits), PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True

class PDHill(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        e0 = parameters[0]
        b = parameters[1]
        g = parameters[2]
        d = parameters[3]

        xprime = xToUse**d

        self.yPredicted = e0 + (b*xprime)/(g**d+xprime)
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Hill (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            d = 1

            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            emin = np.min(yToUse)
            emax = np.max(yToUse)
            e0 = emin

            b = emax - emin
            g = 5 / np.max(xToUse)

            print("First estimate of Hill term: ")
            print("Y  = (%f) + (%f)*exp(-(%f)*(X^(%f))) " % (e0, b, g, d))

            self.bounds = []
            self.bounds.append((min(0.1 * e0, 10 * e0), max(0.1 * e0, 10 * e0)))
            self.bounds.append((min(0.1 * b, 10 * b), max(0.1 * b, 10 * b)))
            self.bounds.append((min(0.1 * g, 10 * g), max(0.1 * g, 10 * g)))
            self.bounds.append((min(0.1 * d, 10 * d), max(0.1 * d, 10 * d)))

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y = e0 + b*X^d/(g^d+X^d)"

    def getEquation(self):
        toPrint = "Y  = (%f) + (%f)*(X^(%f))/((%f)^(%f)+X^(%f))" % (self.parameters[0], self.parameters[1],
                                                                    self.parameters[3], self.parameters[2],
                                                                    self.parameters[3], self.parameters[3])
        return toPrint

    def getParameterNames(self):
        return ['e0', 'b', 'g', 'd']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=e0 + b*X^d/(g^d+X^d)'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = [yunits, PKPDUnit.UNIT_NONE, xunits, PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        retval.append(lowerBound[2] > 0 or upperBound[2] < 0)
        retval.append(lowerBound[3] > 0 or upperBound[3] < 0)
        return retval

    def areParametersValid(self, p):
        return True

class PDOQuigley0(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        a = parameters[0]

        self.yPredicted = ((np.tanh(xToUse)+1)/2)**a
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "OQuigley0 (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            a = 1

            print("First estimate of OQuigley0 term: ")
            print("Y  = ((tanh(X)+1)/2)^(%f) " % (a))

            self.bounds = []
            self.bounds.append((-5, 5))

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y=((tanh(X)+1)/2)^a"

    def getEquation(self):
        toPrint = "Y  = ((tanh(X)+1)/2)^(%f)" % self.parameters[0]
        return toPrint

    def getParameterNames(self):
        return ['a']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=((tanh(X)+1)/2)^a'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        self.parameterUnits = [PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        return retval

    def areParametersValid(self, p):
        return True

class PDOQuigley1(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        x0 = parameters[0]
        a = parameters[1]

        self.yPredicted = ((np.tanh(xToUse-x0)+1)/2)**a
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "OQuigley1 (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            a = 1
            x0=min(xToUse)

            print("First estimate of OQuigley1 term: ")
            print("Y  = ((tanh(X-(%f))+1)/2)^(%f) " % (x0,a))

            self.bounds = []
            self.bounds.append((min(0.1*x0), max(max(xToUse),1e3*x0)))
            self.bounds.append((-5, 5))

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y=((tanh(X-X0)+1)/2)^a"

    def getEquation(self):
        toPrint = "Y  = ((tanh(X-(%f))+1)/2)^(%f)" % (self.parameters[0], self.parameters[1])
        return toPrint

    def getParameterNames(self):
        return ['x0', 'a']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=((tanh(X-X0)+1)/2)^a'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        self.parameterUnits = [PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        return retval

    def areParametersValid(self, p):
        return True

class PDOQuigley2(PDGenericModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        x0 = parameters[0]
        g = parameters[1]
        expx = np.exp(g * (xToUse - x0))

        self.yPredicted = expx/(1+expx)
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "OQuigley2 (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            g = 1
            x0=min(xToUse)

            print("First estimate of OQuigley2 term: ")
            print("Y  = exp(g*(X-(%f)))/(1+exp((%f)*(X-(%f))) " % (x0,g,x0))

            self.bounds = []
            self.bounds.append((min(0.1*x0), max(max(xToUse),1e3*x0)))
            self.bounds.append((-5, 5))

    def printSetup(self):
        print("Model: %s" % self.getModelEquation())
        print("Bounds: " + str(self.bounds))

    def getModelEquation(self):
        return "Y=exp(g*(X-X0))/(1+exp(g*(X-X0))"

    def getEquation(self):
        toPrint = "Y=exp(g*(X-(%f)))/(1+exp((%f)*(X-(%f)))" % (self.parameters[0], self.parameters[1], self.parameters[0])
        return toPrint

    def getParameterNames(self):
        return ['x0', 'g']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=exp(g*(X-X0))/(1+exp(g*(X-X0))'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        self.parameterUnits = [PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        return retval

    def areParametersValid(self, p):
        return True
