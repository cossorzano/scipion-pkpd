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
from pkpd.pkpd_units import inverseUnits, divideUnits, PKPDUnit

import math

class DissolutionModel(PKPDModel):
    def __init__(self):
        self.allowTlag = False

    def setAllowTLag(self, _allowTlag):
        self.allowTlag = _allowTlag

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form %s'%self.getModelEquation()]*self.getNumberOfParameters()

    def areParametersValid(self, p):
        return np.sum(p>0)==p.size # All positive

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        for i in range(len(upperBound)):
            retval.append(lowerBound[i]>0 or upperBound[i]<0)
        return retval


class Dissolution0(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        if self.allowTlag:
            tlag=parameters[0]
            K=parameters[1]
        else:
            tlag=0.0
            K=parameters[0]
        xToUse=np.clip(xToUse-tlag,0.0,None) # u(t-tlag)

        self.yPredicted = K*xToUse
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Zero order dissolution (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            K = np.max(yToUse)/np.max(xToUse)
            print("First estimate of zero order dissolution: ")
            print("Y=%f*t"%K)

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((0.1*K,10*K))

    def getModelEquation(self):
        if self.allowTlag:
            return "Y=K*(t-tlag)"
        else:
            return "Y=K*t"

    def getEquation(self):
        if self.allowTlag:
            toPrint="Y=(%f)*(t-(%f))"%(self.parameters[1],self.parameters[0])
        else:
            toPrint="Y=(%f)*t"%self.parameters[0]
        return toPrint

    def getParameterNames(self):
        if self.allowTlag:
            return ['tlag','K']
        else:
            return ['K']

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        sunits = divideUnits(yunits,xunits)
        if self.allowTlag:
            self.parameterUnits=[xunits, sunits]
        else:
            self.parameterUnits=[sunits]
        return self.parameterUnits


class Dissolution1(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        if self.allowTlag:
            tlag=parameters[0]
            Vmax=parameters[1]
            beta=parameters[2]
        else:
            tlag=0.0
            Vmax=parameters[0]
            beta=parameters[1]
        xToUse=np.clip(xToUse-tlag,0.0,None) # u(t-tlag)

        self.yPredicted = Vmax*(1-np.exp(-beta*xToUse))
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "1st order dissolution (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Vmax = np.max(yToUse)
            beta=5/np.max(xToUse)
            print("First estimate of first order dissolution: ")
            print("Y=%f*(1-exp(-(%f)*t)"%(Vmax,beta))

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((0.1*Vmax,10*Vmax))
            self.bounds.append((0.1*beta,10*beta))

    def getModelEquation(self):
        if self.allowTlag:
            return "Y=Vmax*(1-exp(-beta*(t-tlag)))"
        else:
            return "Y=Vmax*(1-exp(-beta*t))"

    def getEquation(self):
        if self.allowTlag:
            toPrint="Y=(%f)*(1-exp(-(%f)*(t-(%f)))"%(self.parameters[1],self.parameters[2],self.parameters[0])
        else:
            toPrint="Y=(%f)*(1-exp(-(%f)*t))"%(self.parameters[0],self.parameters[1])
        return toPrint

    def getParameterNames(self):
        if self.allowTlag:
            return ['tlag','Vmax','beta']
        else:
            return ['Vmax','beta']

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        x1units = inverseUnits(xunits)
        if self.allowTlag:
            self.parameterUnits=[xunits, yunits, x1units]
        else:
            self.parameterUnits=[yunits, x1units]
        return self.parameterUnits


class DissolutionAlpha(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        if self.allowTlag:
            tlag=parameters[0]
            Vmax=parameters[1]
            beta=parameters[2]
            alpha=parameters[3]
        else:
            tlag=0.0
            Vmax=parameters[0]
            beta=parameters[1]
            alpha=parameters[2]
        xToUse=np.clip(xToUse-tlag,0.0,None) # u(t-tlag)
        alpha=np.clip(0.00001,alpha,None)
        argument = np.clip(pow(Vmax,alpha)-alpha*beta*xToUse,0.0,None)

        self.yPredicted = Vmax - np.power(argument,1/alpha)
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Fractional order dissolution (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Vmax = np.max(yToUse)
            beta=Vmax/np.max(xToUse)
            print("First estimate of fractional order dissolution: ")
            print("Y=%f-pow(%f^1-1*%f*t,1)"%(Vmax,Vmax,beta))

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((0.1*Vmax,10*Vmax))
            self.bounds.append((0.1*beta,10*beta))
            self.bounds.append((0.0001,3.0))

    def getModelEquation(self):
        if self.allowTlag:
            return "Y=Vmax-pow(Vmax^alpha-alpha*beta*(t-tlag),1/alpha)"
        else:
            return "Y=Vmax-pow(Vmax^alpha-alpha*beta*t,1/alpha)"

    def getEquation(self):
        if self.allowTlag:
            toPrint="Y=(%f)-pow((%f)^(%f)-(%f)*(%f)*(t-(%f)),1/(%f))"%(self.parameters[1],self.parameters[1],self.parameters[3],
                                                                       self.parameters[3],self.parameters[2],self.parameters[0],
                                                                       self.parameters[3])
        else:
            toPrint="Y=(%f)-pow((%f)^(%f)-(%f)*(%f)*t),1/(%f))"%(self.parameters[0],self.parameters[0],self.parameters[2],
                                                                 self.parameters[2],self.parameters[1],self.parameters[2])
        return toPrint

    def getParameterNames(self):
        if self.allowTlag:
            return ['tlag','Vmax','beta','alpha']
        else:
            return ['Vmax','beta','alpha']

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        x1units = inverseUnits(xunits)
        if self.allowTlag:
            self.parameterUnits=[xunits, yunits, x1units, PKPDUnit.UNIT_NONE]
        else:
            self.parameterUnits=[yunits, x1units, PKPDUnit.UNIT_NONE]
        return self.parameterUnits


class DissolutionWeibull(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        if self.allowTlag:
            tlag=parameters[0]
            Vmax=parameters[1]
            lambdda=parameters[2]
            b=parameters[3]
        else:
            tlag=0.0
            Vmax=parameters[0]
            lambdda=parameters[1]
            b=parameters[2]
        xToUse=np.clip(xToUse-tlag,0.0,None) # u(t-tlag)
        argument = lambdda*np.power(xToUse,b)

        self.yPredicted = Vmax*(1-np.exp(-argument))
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Weibull dissolution (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Vmax = np.max(yToUse)
            lambdda=5/np.max(xToUse)
            print("First estimate of Weibull order dissolution: ")
            print("Y=(%f)*(1-exp(-(%f)*t)"%(Vmax,lambdda))

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((0.1*Vmax,10*Vmax))
            self.bounds.append((0.1*lambdda,10*lambdda))
            self.bounds.append((0.0001,5.0))

    def getModelEquation(self):
        if self.allowTlag:
            return "Y=Vmax*(1-exp(-lambda*(t-tlag)^b)"
        else:
            return "Y=Vmax*(1-exp(-lambda*t^b)"

    def getEquation(self):
        if self.allowTlag:
            toPrint="Y=(%f)*(1-exp(-(%f)*(t-(%f))^(%f)))"%(self.parameters[1],self.parameters[2],self.parameters[0],
                                                    self.parameters[3])
        else:
            toPrint="Y=(%f)*(1-exp(-(%f)*t^(%f)))"%(self.parameters[0],self.parameters[1],self.parameters[2])
        return toPrint

    def getParameterNames(self):
        if self.allowTlag:
            return ['tlag','Vmax','lambda','b']
        else:
            return ['Vmax','lambda','b']

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        x1units = inverseUnits(xunits)
        if self.allowTlag:
            self.parameterUnits=[xunits, yunits, x1units, PKPDUnit.UNIT_NONE]
        else:
            self.parameterUnits=[yunits, x1units, PKPDUnit.UNIT_NONE]
        return self.parameterUnits


class DissolutionHiguchi(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        if self.allowTlag:
            tlag=parameters[0]
            Vmax=parameters[1]
        else:
            tlag=0.0
            Vmax=parameters[0]
        xToUse=np.clip(xToUse-tlag,0.0,None) # u(t-tlag)

        self.yPredicted = Vmax*np.sqrt(xToUse)
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Higuchi dissolution (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Vmax = np.max(yToUse)/np.sqrt(np.max(xToUse))
            print("First estimate of Higuchi order dissolution: ")
            print("Y=(%f)*t^0.5"%Vmax)

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((0.1*Vmax,10*Vmax))

    def getModelEquation(self):
        if self.allowTlag:
            return "Y=Vmax*(t-tlag)^0.5"
        else:
            return "Y=Vmax*t^0.5"

    def getEquation(self):
        if self.allowTlag:
            toPrint="Y=(%f)*(t-(%f))^0.5"%(self.parameters[1],self.parameters[0])
        else:
            toPrint="Y=(%f)*t^0.5"%(self.parameters[0])
        return toPrint

    def getParameterNames(self):
        if self.allowTlag:
            return ['tlag','Vmax']
        else:
            return ['Vmax']

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        if self.allowTlag:
            self.parameterUnits=[xunits, yunits]
        else:
            self.parameterUnits=[yunits]
        return self.parameterUnits


class DissolutionKorsmeyer(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        if self.allowTlag:
            tlag=parameters[0]
            Vmax=parameters[1]
            m=parameters[2]
        else:
            tlag=0.0
            Vmax=parameters[0]
            m = parameters[1]
        xToUse=np.clip(xToUse-tlag,0.0,None) # u(t-tlag)

        self.yPredicted = Vmax*np.power(xToUse,m)
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Korsmeyer-Peppas dissolution (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Vmax = np.max(yToUse)/np.sqrt(np.max(xToUse))
            print("First estimate of Korsmeyer-Peppas order dissolution: ")
            print("Y=(%f)*t^0.5"%Vmax)

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((0.1*Vmax,10*Vmax))
            self.bounds.append((0.0001,5.0))

    def getModelEquation(self):
        if self.allowTlag:
            return "Y=Vmax*(t-tlag)^m"
        else:
            return "Y=Vmax*t^m"

    def getEquation(self):
        if self.allowTlag:
            toPrint="Y=(%f)*(t-(%f))^(%f)"%(self.parameters[1],self.parameters[0],self.parameters[2])
        else:
            toPrint="Y=(%f)*t^(%f)"%(self.parameters[0],self.parameters[1])
        return toPrint

    def getParameterNames(self):
        if self.allowTlag:
            return ['tlag','Vmax','m']
        else:
            return ['Vmax','m']

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        if self.allowTlag:
            self.parameterUnits=[xunits, yunits, PKPDUnit.UNIT_NONE]
        else:
            self.parameterUnits=[yunits, PKPDUnit.UNIT_NONE]
        return self.parameterUnits


class DissolutionHixson(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        if self.allowTlag:
            tlag=parameters[0]
            Vmax=parameters[1]
            K=parameters[2]
        else:
            tlag=0.0
            Vmax=parameters[0]
            K = parameters[1]
        xToUse=np.clip(xToUse-tlag,0.0,None) # u(t-tlag)

        self.yPredicted = Vmax*(1-np.power(1-K*xToUse,3))
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Hixson-Crowell dissolution (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Vmax = np.max(yToUse)/np.sqrt(np.max(xToUse))
            K=0.5/np.max(xToUse)
            print("First estimate of Hixson-Crowell order dissolution: ")
            print("Y=(%f)*(1-(1-(%f)*t)^3)"%(Vmax,K))

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((0.1*Vmax,10*Vmax))
            self.bounds.append((0.1*K,10*K))

    def getModelEquation(self):
        if self.allowTlag:
            return "Y=Vmax*(1-(1-K*(t-tlag))^3)"
        else:
            return "Y=Vmax*(1-(1-K*t)^3)"

    def getEquation(self):
        if self.allowTlag:
            toPrint="Y=(%f)*(1-(1-(%f)*(t-(%f)))^3)"%(self.parameters[1],self.parameters[2],self.parameters[0])
        else:
            toPrint="Y=(%f)*(1-(1-(%f)*t)^3)"%(self.parameters[0],self.parameters[1])
        return toPrint

    def getParameterNames(self):
        if self.allowTlag:
            return ['tlag','Vmax','K']
        else:
            return ['Vmax','K']

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        x1units = inverseUnits(xunits)
        if self.allowTlag:
            self.parameterUnits=[xunits, yunits, x1units]
        else:
            self.parameterUnits=[yunits, x1units]
        return self.parameterUnits


class DissolutionHopfenberg(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        if self.allowTlag:
            tlag=parameters[0]
            Vmax=parameters[1]
            K=parameters[2]
            m=parameters[3]
        else:
            tlag=0.0
            Vmax=parameters[0]
            K = parameters[1]
            m = parameters[2]
        xToUse=np.clip(xToUse-tlag,0.0,None) # u(t-tlag)

        self.yPredicted = Vmax*(1-np.power(1-K*xToUse,m))
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Hopfenberg dissolution (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Vmax = np.max(yToUse)/np.sqrt(np.max(xToUse))
            K=0.5/np.max(xToUse)
            print("First estimate of Hopfenberg order dissolution: ")
            print("Y=(%f)*(1-(1-(%f)*t)^3)"%(Vmax,K))

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((0.1*Vmax,10*Vmax))
            self.bounds.append((0.1*K,10*K))
            self.bounds.append((0.0001,20.0))

    def getModelEquation(self):
        if self.allowTlag:
            return "Y=Vmax*(1-(1-K*(t-tlag))^m)"
        else:
            return "Y=Vmax*(1-(1-K*t)^m)"

    def getEquation(self):
        if self.allowTlag:
            toPrint="Y=(%f)*(1-(1-(%f)*(t-(%f)))^(%f))"%(self.parameters[1],self.parameters[2],self.parameters[0],
                                                         self.parameters[3])
        else:
            toPrint="Y=(%f)*(1-(1-(%f)*t)^(%f))"%(self.parameters[0],self.parameters[1],self.parameters[2])
        return toPrint

    def getParameterNames(self):
        if self.allowTlag:
            return ['tlag','Vmax','K','m']
        else:
            return ['Vmax','K','m']

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        x1units = inverseUnits(xunits)
        if self.allowTlag:
            self.parameterUnits=[xunits, yunits, x1units, PKPDUnit.UNIT_NONE]
        else:
            self.parameterUnits=[yunits, x1units, PKPDUnit.UNIT_NONE]
        return self.parameterUnits
