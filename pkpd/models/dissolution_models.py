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

class DissolutionModel(PKPDModel):
    def __init__(self):
        self.allowTlag = False

    def setAllowTLag(self, _allowTlag):
        self.allowTlag = _allowTlag

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
        xToUse=np.clip(xToUse,0.0,None) # u(t-tlag)

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

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())

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

    def getParameterDescriptions(self):
        if self.allowTlag:
            eqStr='Y=K*(t-tlag)'
        else:
            eqStr='Y=K*t'
        return ['Automatically fitted model of the form %s'%eqStr]*self.getNumberOfParameters()

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        sunits = divideUnits(yunits,xunits)
        if self.allowTlag:
            self.parameterUnits=[xunits, sunits]
        else:
            self.parameterUnits=[sunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        if self.allowTlag:
            retval.append(lowerBound[0]>0 or upperBound[0]<0)
            retval.append(lowerBound[1]>0 or upperBound[1]<0)
        else:
            retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
        return retval

    def areParametersValid(self, p):
        return np.sum(p>0)==p.size # All positive


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
        xToUse=np.clip(xToUse,0.0,None) # u(t-tlag)

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
            print("First estimate of zero order dissolution: ")
            print("Y=%f*(1-exp(-(%f)*t)"%(Vmax,beta))

            self.bounds = []
            if self.allowTlag:
                self.bounds.append((0.0,np.max(xToUse)))
            self.bounds.append((0.1*Vmax,10*Vmax))
            self.bounds.append((0.1*beta,10*beta))

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())

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

    def getParameterDescriptions(self):
        if self.allowTlag:
            eqStr='Y=Vmax*(1-exp(-beta*(t-tlag)))'
        else:
            eqStr='Y=Vmax*(1-exp(-beta*t))'
        return ['Automatically fitted model of the form %s'%eqStr]*self.getNumberOfParameters()

    def calculateParameterUnits(self,sample):
        yunits = self.experiment.getVarUnits(self.yName)
        xunits = self.experiment.getVarUnits(self.xName)
        x1units = inverseUnits(xunits)
        if self.allowTlag:
            self.parameterUnits=[xunits, yunits, x1units]
        else:
            self.parameterUnits=[yunits, x1units]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        if self.allowTlag:
            retval.append(lowerBound[0]>0 or upperBound[0]<0)
            retval.append(lowerBound[1]>0 or upperBound[1]<0)
            retval.append(lowerBound[2]>0 or upperBound[2]<0)
        else:
            retval.append(lowerBound[0] > 0 or upperBound[0] < 0)
            retval.append(lowerBound[1] > 0 or upperBound[1] < 0)
        return retval

    def areParametersValid(self, p):
        return np.sum(p>0)==p.size # All positive

