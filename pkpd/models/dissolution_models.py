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

import copy
import numpy as np

from pkpd.objects import PKPDModel
from pkpd.pkpd_units import inverseUnits, divideUnits, PKPDUnit
from pkpd.utils import uniqueFloatValues
from scipy.interpolate import InterpolatedUnivariateSpline, PchipInterpolator

# Tested by test_workflow_dissolution

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
        if x is None:
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
        if x is None:
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
        if x is None:
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
        if x is None:
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


class DissolutionDoubleWeibull(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x is None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        if self.allowTlag:
            tlag=parameters[0]
            idx=1
        else:
            tlag=0.0
            idx=0
        Vmax=parameters[idx]
        lambdda1=parameters[idx+1]
        b1=parameters[idx+2]
        F1=parameters[idx+3]
        tlag2=parameters[idx+4]
        lambdda2=parameters[idx+5]
        b2=parameters[idx+6]
        xToUse=np.clip(xToUse-tlag,0.0,None) # u(t-tlag)
        argument1 = lambdda1*np.power(xToUse,b1)
        xToUse2=np.clip(xToUse-tlag-tlag2,0.0,None)
        argument2 = lambdda2*np.power(xToUse2,b2)

        self.yPredicted = Vmax*(F1*(1-np.exp(-argument1))+(1-F1)*(1-np.exp(-argument2)))
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Double Weibull dissolution (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Vmax = np.max(yToUse)
            tmax=np.max(xToUse)
            lambdda=5/tmax
            print("First estimate of Double Weibull order dissolution: ")
            print("Y=(%f)*(1-exp(-(%f)*t)"%(Vmax,lambdda))

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((0.1*Vmax,10*Vmax)) # Vmax
            self.bounds.append((0.1*lambdda,10*lambdda)) # lambda1
            self.bounds.append((0.0001,5.0)) # b1
            self.bounds.append((0,1)) # F1
            self.bounds.append((0.0,tmax)) # tlag2
            self.bounds.append((0.1*lambdda,10*lambdda)) # lambda2
            self.bounds.append((0.0001,5.0)) # b2

    def getModelEquation(self):
        if self.allowTlag:
            tStr="(t-tlag)"
        else:
            tStr="t"
        return "Y=Vmax*(F1*(1-exp(-lambda1*%s^b1))+(1-F1)*(1-exp(-lambda2*(%s-tlag2)^b2)))"%(tStr,tStr)

    def getEquation(self):
        if self.allowTlag:
            tStr="(t-(%f))"%self.parameters[0]
            idx=1
        else:
            tStr="t"
            idx=0
        Vmax=self.parameters[idx]
        lambdda1=self.parameters[idx+1]
        b1=self.parameters[idx+2]
        F1=self.parameters[idx+3]
        tlag2=self.parameters[idx+4]
        lambdda2=self.parameters[idx+5]
        b2=self.parameters[idx+6]
        toPrint="Y=(%f)*((%f)*(1-exp(-(%f)*%s^(%f)))+(%f)*(1-exp(-(%f)*(%s-(%f))^(%f))))"%\
                (Vmax,F1,lambdda1,tStr,b1,1-F1,lambdda2,tStr,tlag2,b2)
        return toPrint

    def getParameterNames(self):
        if self.allowTlag:
            return ['tlag','Vmax','lambda1','b1','F1','tlag2','lambda2','b2']
        else:
            return ['Vmax','lambda1','b1','F1','tlag2','lambda2','b2']

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        x1units = inverseUnits(xunits)
        self.parameterUnits = []
        if self.allowTlag:
            self.parameterUnits.append(xunits)
        self.parameterUnits+=[yunits, x1units, PKPDUnit.UNIT_NONE, PKPDUnit.UNIT_NONE, xunits, x1units, PKPDUnit.UNIT_NONE]
        return self.parameterUnits


class DissolutionHiguchi(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x is None:
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
        if x is None:
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
        if x is None:
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

        self.yPredicted = Vmax*(1-np.power(np.clip(1-K*xToUse,0.0,None),3))
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
        if x is None:
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

        self.yPredicted = Vmax*(1-np.power(np.clip(1-K*xToUse,0.0,None),m))
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


class DissolutionHill(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x == None:
            x = self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        if self.allowTlag:
            tlag=parameters[0]
            Vmax = parameters[1]
            g = parameters[2]
            d = parameters[3]
        else:
            tlag=0.0
            Vmax = parameters[0]
            g = parameters[1]
            d = parameters[2]
        xToUse=np.clip(xToUse-tlag,0.0,None) # u(t-tlag)

        self.yPredicted = np.zeros(xToUse.shape[0])

        xprime = xToUse**d

        self.yPredicted = (Vmax*xprime)/(g**d+xprime)
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Hill (%s)" % self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            d = 1

            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Vmax = np.max(yToUse)
            g = 5 / np.max(xToUse)

            print("First estimate of Hill term: ")
            print("Y  = (%f)*(t^(%f))/((%f)^(%f)+t^(%f))" % (Vmax,d,g,d,d))

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((min(0.1 * Vmax, 10 * Vmax), max(0.1 * Vmax, 10 * Vmax)))
            self.bounds.append((min(0.1 * g, 10 * g), max(0.1 * g, 10 * g)))
            self.bounds.append((min(0.1 * d, 10 * d), max(0.1 * d, 10 * d)))

    def getModelEquation(self):
        if self.allowTlag:
            return "Y = Vmax*(t-tlag)^d/(g^d+(t-tlag)^d)"
        else:
            return "Y = Vmax*t^d/(g^d+t^d)"

    def getEquation(self):
        if self.allowTlag:
            toPrint = "Y  = (%f)*((t-(%f))^(%f))/((%f)^(%f)+(t-(%f))^(%f))" % (self.parameters[1],self.parameters[0],
                                                                 self.parameters[3], self.parameters[2],
                                                                 self.parameters[3], self.parameters[0], self.parameters[3])
        else:
            toPrint = "Y  = (%f)*(t^(%f))/((%f)^(%f)+t^(%f))" % (self.parameters[0],
                                                                 self.parameters[2], self.parameters[1],
                                                                 self.parameters[2], self.parameters[2])
        return toPrint

    def getParameterNames(self):
        if self.allowTlag:
            return ['tlag', 'Vmax', 'g', 'd']
        else:
            return ['Vmax', 'g', 'd']

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form Y=Vmax*t^d/(g^d+t^d)'] * self.getNumberOfParameters()

    def calculateParameterUnits(self, sample):
        xunits = self.experiment.getVarUnits(self.xName)
        if self.allowTlag:
            self.parameterUnits = [xunits, PKPDUnit.UNIT_NONE, xunits, PKPDUnit.UNIT_NONE]
        else:
            self.parameterUnits = [PKPDUnit.UNIT_NONE, xunits, PKPDUnit.UNIT_NONE]
        return self.parameterUnits


class DissolutionMakoidBanakar(DissolutionModel):
    def forwardModel(self, parameters, x=None):
        if x is None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        if self.allowTlag:
            tlag=parameters[0]
            idx=1
        else:
            tlag=0.0
            idx=0
        Vmax=parameters[idx]
        Tmax=parameters[idx+1]
        b=parameters[idx+2]
        xToUse=np.clip(xToUse-tlag,0.0,None) # u(t-tlag)
        xToUseTmax=xToUse/Tmax

        self.yPredicted = np.clip(Vmax*np.power(xToUseTmax,b)*np.exp(b*(1-xToUseTmax)),0.0,Vmax)
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Makoid-Banakar dissolution (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Vmax = np.max(yToUse)
            tmax=np.max(xToUse)
            b=1.0
            print("First estimate of Makoid-Banakar order dissolution: ")
            print("Y=(%f)*(t/%f)*exp(-t/%f)"%(Vmax,tmax,tmax))

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((0.1*Vmax,10*Vmax)) # Vmax
            self.bounds.append((0.0,tmax)) # tmax
            self.bounds.append((0.0001,5.0)) # b

    def getModelEquation(self):
        if self.allowTlag:
            tStr="(t-tlag)/tmax"
        else:
            tStr="t/tmax"
        return "Y=Vmax*((%s)^b*exp(b*(1-%s)) if %s<1 and Vmax if %s>=1"%(tStr,tStr,tStr,tStr)

    def getEquation(self):
        if self.allowTlag:
            tStr="(t-(%f))/(%f)"%(self.parameters[0],self.parameters[2])
            idx=1
        else:
            tStr="t/(%f)"%self.parameters[1]
            idx=0
        Vmax=self.parameters[idx]
        b=self.parameters[idx+2]
        toPrint="Y=(%f)*(%s)^(%s)*exp((%f)*(1-%s)) if %s<1 and %f if %s>=1"%\
                (Vmax,tStr,b,b,tStr,tStr,Vmax,tStr)
        return toPrint

    def getParameterNames(self):
        if self.allowTlag:
            return ['tlag','Vmax','tmax','b']
        else:
            return ['Vmax','tmax','b']

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits = []
        if self.allowTlag:
            self.parameterUnits.append(xunits)
        self.parameterUnits+=[yunits, xunits, PKPDUnit.UNIT_NONE]
        return self.parameterUnits


class DissolutionSplinesGeneric(DissolutionModel):
    def __init__(self):
        self.nknots=0
        self.parametersPrepared=None

    def rearrange(self,parameters):
        retval = parameters
        if self.allowTlag:
            retval[3:]=np.sort(retval[3:])
        else:
            retval[2:]=np.sort(retval[2:])
        return retval

    def forwardModel(self, parameters, x=None):
        if x is None:
            x=self.x
        xToUse = x[0] if type(x)==list else x # From [array(...)] to array(...)
        xToUse = np.copy(xToUse) # Just in case it is modified
        self.yPredicted = np.zeros(xToUse.shape[0])

        if self.allowTlag:
            tlag=parameters[0]
            Vmax=parameters[1]
            tmax=parameters[2]
            coefs=parameters[3:]
            xToUse-=tlag
        else:
            Vmax=parameters[0]
            tmax=parameters[1]
            coefs=parameters[2:]

        xToUse = np.clip(xToUse, 0.0, tmax)

        if self.parametersPrepared is None or not np.array_equal(self.parametersPrepared,parameters):
            self.knots = np.linspace(0, tmax, self.nknots+2)
            self.knotsY = np.append(np.insert(coefs,0,0),1)
            self.knotsY=np.sort(self.knotsY)
            knotsUnique, knotsYUnique=uniqueFloatValues(self.knots, self.knotsY)
            # self.B=InterpolatedUnivariateSpline(knotsUnique, knotsYUnique, k=1)
            self.B=PchipInterpolator(knotsUnique, knotsYUnique)
            self.parametersPrepared=copy.copy(parameters)

        fraction=self.B(xToUse)
        fraction=np.clip(fraction,0.0,1.0)
        self.yPredicted = Vmax*fraction
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Splines dissolution (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse=self.x[0] # From [array(...)] to array(...)
            yToUse=self.y[0] # From [array(...)] to array(...)
            Vmax = np.max(yToUse)

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)/2))
            self.bounds.append((0.1*Vmax,10*Vmax))
            self.bounds.append((0.0, np.max(xToUse)))
            for i in range(self.nknots):
                self.bounds.append((0.0,1.0))

    def getModelEquation(self):
        if self.allowTlag:
            return "Y=Vmax*Bspline(t-tlag,%d,tmax)"%self.nknots
        else:
            return "Y=Vmax*Bspline(t,%d,tmax)"%self.nknots

    def getEquation(self):
        if self.allowTlag:
            tlag=self.parameters[0]
            Vmax=self.parameters[1]
            tmax=self.parameters[2]
            coefs=self.parameters[3:]
            toPrint="Y=(%f)*Bspline(t-(%f),%d,%f,coefs=%s)"%(Vmax,tlag,self.nknots,tmax,np.array2string(coefs,max_line_width=1000))
        else:
            Vmax=self.parameters[0]
            tmax=self.parameters[1]
            coefs=self.parameters[2:]
            toPrint="Y=(%f)*Bspline(t,%d,%f,coefs=%s)"%(Vmax,self.nknots,tmax,np.array2string(coefs,max_line_width=1000))
        return toPrint

    def getParameterNames(self):
        retval=['Vmax','tmax']
        if self.allowTlag:
            retval=['tlag']+retval
        retval+=['dissol_spline%d_A%d'%(self.nknots,i) for i in range(self.nknots)]
        return retval

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        self.parameterUnits=[]
        if self.allowTlag:
            self.parameterUnits+=[xunits]
        self.parameterUnits+=[yunits,xunits]+[PKPDUnit.UNIT_NONE]*(self.nknots)
        return self.parameterUnits

class DissolutionSplines2(DissolutionSplinesGeneric):
    def __init__(self):
        DissolutionSplinesGeneric.__init__(self)
        self.nknots = 2

class DissolutionSplines3(DissolutionSplinesGeneric):
    def __init__(self):
        DissolutionSplinesGeneric.__init__(self)
        self.nknots = 3

class DissolutionSplines4(DissolutionSplinesGeneric):
    def __init__(self):
        DissolutionSplinesGeneric.__init__(self)
        self.nknots = 4

class DissolutionSplines5(DissolutionSplinesGeneric):
    def __init__(self):
        DissolutionSplinesGeneric.__init__(self)
        self.nknots = 5

class DissolutionSplines6(DissolutionSplinesGeneric):
    def __init__(self):
        DissolutionSplinesGeneric.__init__(self)
        self.nknots = 6

class DissolutionSplines7(DissolutionSplinesGeneric):
    def __init__(self):
        DissolutionSplinesGeneric.__init__(self)
        self.nknots = 7

class DissolutionSplines8(DissolutionSplinesGeneric):
    def __init__(self):
        DissolutionSplinesGeneric.__init__(self)
        self.nknots = 8

class DissolutionSplines9(DissolutionSplinesGeneric):
    def __init__(self):
        DissolutionSplinesGeneric.__init__(self)
        self.nknots = 9

class DissolutionSplines10(DissolutionSplinesGeneric):
    def __init__(self):
        DissolutionSplinesGeneric.__init__(self)
        self.nknots = 10
