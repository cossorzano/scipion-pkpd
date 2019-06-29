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
PK models
"""
import math

import numpy as np

from pkpd.objects import PKPDModel, PKPDODEModel
from pkpd.pkpd_units import (inverseUnits, divideUnits, multiplyUnits,
                             unitFromString, createUnit, PKPDUnit)


class PKModel(PKPDModel):
    pass

class PKGenericModel(PKModel):
    pass

def signifcantMessage1(varName,lower,upper):
    if lower < 0 and upper > 0:
        return "Suspicious, %s looks like a constant"%varName
    elif lower < 0:
        return "Suspicious, %s may be unstable"%varName
    elif lower > 0:
        return "True"
    else:
        return "NA"

def signifcantMessage2(varName,lower,upper):
    if lower < 0 and upper > 0:
        return "Suspicious, %s looks like 0"%varName
    elif upper < 0:
        return "Suspicious, %s seems to be negative"%varName
    else:
        return "True"

def signifcantMessage3(varName,lower,upper):
    if lower < 0 and upper > 0:
        return "Suspicious, %s looks like 0"%varName
    else:
        return "True"

class PKPDExponentialModel(PKGenericModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        xToUse = x[0] # From [array(...)] to array(...)
        self.yPredicted = np.zeros(xToUse.shape[0])
        proceed=True
        for k in range(1,self.Nexp):
            ck = parameters[2*k]
            ck_1 = parameters[2*(k-1)]
            if ck_1<ck:
                proceed=False
                self.yPredicted = -1000*np.ones(xToUse.shape[0])
                break
        if proceed:
            for k in range(0,self.Nexp):
                ck = parameters[2*k]
                lk = parameters[2*k+1]
                self.yPredicted += ck*np.exp(-lk*xToUse)
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Sum of exponentials (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse = self.x[0] # From [array(...)] to array(...)
            ylogToUse = self.ylog[0]
            p = np.polyfit(xToUse,ylogToUse,1)
            print("First estimate of 1 exponential term: ")
            print("Y=%f*exp(-%f*X)"%(math.exp(p[1]),-p[0]))

            cBound = (math.exp(p[1])*0.01,math.exp(p[1])*100.0)
            lambdaBound = (-p[0]*0.01,-p[0]*100.0)
            self.bounds = []
            for i in range(self.Nexp):
                self.bounds.append(cBound)
                self.bounds.append(lambdaBound)

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Number of exponentials: "+str(self.Nexp))
        print("Bounds: "+str(self.bounds))

    def getModelEquation(self):
        return "Y=sum_i c_i*exp(-lambda_i * X)"

    def getEquation(self):
        toPrint="Y="
        for i in range(self.Nexp):
            toPrint+= "+[%f*exp(-%f*X)]"%(self.parameters[2*i],self.parameters[2*i+1])
        return toPrint

    def getParameterNames(self):
        parameterList = []
        for i in range(self.Nexp):
            parameterList.append('c%d'%(i+1))
            parameterList.append('lambda%d'%(i+1))
        return parameterList

    def calculateParameterUnits(self,sample):
        xunits = self.experiment.getVarUnits(self.xName)
        yunits = self.experiment.getVarUnits(self.yName)
        self.parameterUnits = [yunits,inverseUnits(xunits)]*self.Nexp
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        for i in range(self.Nexp):
            cLower = lowerBound[2*i]
            cUpper = upperBound[2*i]
            if cLower<0 and cUpper>0:
                retval.append("False")
            elif cLower>0:
                retval.append("True")
            elif cUpper<0:
                retval.append("Suspicious, this term may be negative")
            else:
                retval.append("NA")

            decayLower = lowerBound[2*i+1]
            decayUpper = upperBound[2*i+1]
            if decayLower<0 and decayUpper>0:
                retval.append("Suspicious, looks like a constant")
            elif decayLower<0:
                retval.append("Suspicious, this term may be unstable")
            elif decayLower>0:
                retval.append("True")
            else:
                retval.append("NA")
        return retval

    def areParametersValid(self, p):
        return np.sum(p<0)==0

class PKPDSimpleEVModel(PKModel):
    def __init__(self, includeTlag=True):
        PKModel.__init__(self)
        self.includeTlag = includeTlag

    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        xToUse = x[0] # From [array(...)] to array(...)

        Ka=parameters[0]
        Vd=parameters[1]
        if self.includeTlag:
            tlag=parameters[2]
        else:
            tlag = 0.0

        self.yPredicted = np.zeros(xToUse.shape)
        for i in range(xToUse.shape[0]):
            if xToUse[i]>=tlag:
                td=xToUse[i]-tlag
                self.yPredicted[i] = Ka*self.F*self.D/(Vd*(Ka-self.Ke))*(np.exp(-self.Ke*td)-np.exp(-Ka*td))
        self.yPredicted = [self.yPredicted] # From array(...) to [array(...)]
        return self.yPredicted

    def getDescription(self):
        return "Simple non-iv model (%s)"%self.__class__.__name__

    def prepare(self):
        if self.bounds == None:
            xToUse = self.x[0] # From [array(...)] to array(...)
            yToUse = self.y[0]
            # Keep only the ascending part of the curve
            idx = []
            for i in range(1,yToUse.shape[0]):
                idx.append(i-1)
                if yToUse[i-1]>yToUse[i]:
                    break
            if len(idx)<=1:
                print("The first estimate of Ka and Vd cannot be determined")
                self.bounds = None
                return
            xAscending = xToUse[idx]
            yAscending = yToUse[idx]

            ylogE = np.polyval(np.asarray([-self.Ke,math.log(self.C0)],np.double),xAscending)
            ylogToFit = np.log(yAscending)-ylogE
            p = np.polyfit(xAscending,-ylogToFit,1)
            Ka = -p[0]
            Vd = self.D / np.max(self.y) # Incorrect: (Ka*self.F*self.D)/(math.exp(p[1])*(Ka-self.Ke))
            self.bounds=[(Ka*0.2,Ka*5),(0.2*Vd,5*Vd)]
            print("First estimate of Ka: %f"%Ka)
            print("First estimate of Vd: %f"%Vd)
            if self.includeTlag:
                tlag = (p[1]-math.log(self.C0))/(Ka-self.Ke)
                print("First estimate of tlag: %f"%tlag)
                self.bounds.append((-5*abs(tlag),5*abs(tlag)))

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Bounds: "+str(self.bounds))

    def getModelEquation(self):
        if self.includeTlag:
            return "Y=Ka*F*D/(Vd*(Ka-Ke))*(exp(-Ke*(t-tlag))-exp(-Ka*(t-tlag))"
        else:
            return "Y=Ka*F*D/(Vd*(Ka-Ke))*(exp(-Ke*t)-exp(-Ka*t)"

    def getEquation(self):
        Ka=self.parameters[0]
        Vd=self.parameters[1]
        if self.includeTlag:
            tlag=self.parameters[2]
            return "Y=%f*%f*D/(%f*(%f-%f))*(exp(-%f*(t-(%f)))-exp(-%f*(t-(%f))))"%(Ka,self.F,Vd,Ka,self.Ke,self.Ke,tlag,Ka,tlag)
        else:
            return "Y=%f*%f*D/(%f*(%f-%f))*(exp(-%f*t)-exp(-%f*t))"%(Ka,self.F,Vd,Ka,self.Ke,self.Ke,Ka)

    def getParameterNames(self):
        if self.includeTlag:
            return ['Ka','Vd','tlag']
        else:
            return ['Ka','Vd']

    def calculateParameterUnits(self,sample):
        xunits = self.experiment.getVarUnits(self.xName)
        Vunits = divideUnits(self.Dunits,self.C0units)
        self.parameterUnits = [inverseUnits(xunits),Vunits]
        if self.includeTlag:
            self.parameterUnits.append(xunits)
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Ka
        decayLower = lowerBound[0]
        decayUpper = upperBound[0]
        if decayLower<0 and decayUpper>0:
            retval.append("Suspicious, Ka looks like a constant")
        elif decayLower<0:
            retval.append("Suspicious, Ka may be unstable")
        elif decayLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Vd
        VLower = lowerBound[1]
        VUpper = upperBound[1]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, Vd looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, Vd seems to be negative")
        else:
            retval.append("True")

        # tlag
        if self.includeTlag:
            tLower = lowerBound[2]
            tUpper = upperBound[2]
            if tLower<0 and tUpper>0:
                retval.append("Suspicious, tlag looks like 0")
            elif tUpper<0:
                retval.append("Suspicious, tlag seems to be negative")
            else:
                retval.append("True")

        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:1]<0)==0

class PK_Monocompartment(PKPDODEModel):
    def F(self, t, y):
        Cl=self.parameters[0]
        V=self.parameters[1]
        # print("F t=",t," C=",y, " incC=",-Cl/V*y)
        return -Cl/V*y

    def G(self, t, dD):
        V=self.parameters[1]
        # print("G t=",t," dD=",dD," incC=",dD/V)
        return dD/V

    def getResponseDimension(self):
        return 1

    def getStateDimension(self):
        return 1

    def getDescription(self):
        return "Monocompartmental model (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Cl/V * C + 1/V * dD/dt"

    def getEquation(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        return "dC/dt = -(%f)/(%f) * C + 1/(%f) dD/dt"%(Cl,V,V)

    def getParameterNames(self):
        return ['Cl','V']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName)
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [Clunits,Vunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Cl
        ClLower = lowerBound[0]
        ClUpper = upperBound[0]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Cl looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Cl may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # V
        VLower = lowerBound[1]
        VUpper = upperBound[1]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, V looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, V seems to be negative")
        else:
            retval.append("True")
        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:1]<0)==0

    def printOtherParameterization(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        print(" ")
        print("Model parameterization ===========")
        print("The one compartment model (excluding the source) has been estimated according to the homogeneous equation:")
        print("This is called the physiological parameterization")
        print("dC/dt = -Cl/V * C")
        print("where:")
        print("C is the concentration in the central compartment")
        print("with parameters:")
        print("Cl=%f (%s) (Clearance from central compartment)"%(Cl,PKPDUnit.codeToString(self.parameterUnits[0])))
        print("V=%f (%s) (Apparent volume of distribution in the central compartment)"%(V,PKPDUnit.codeToString(self.parameterUnits[1])))
        print(" ")
        print("An alternative formulation is given by its microconstants formulation:")
        print("dA/dt = -ke * A")
        print("where:")
        print("A is the amount of drug in the central compartment (A=C*V)")
        print("with parameters:")
        ke=Cl/V
        kunits = divideUnits(self.parameterUnits[0],self.parameterUnits[1])
        print("ke=%f (%s) (Elimination from central compartment, ke=Cl/V)"%(ke,PKPDUnit.codeToString(kunits)))
        print(" ")
        print("An alternative formulation is given by its impulse response (response to an intravenous bolus of amount A0):")
        print("A=A0*exp(-ke*t)")
        print("where:")
        print("A is the amount of drug in the central compartment")
        print("A0 is the amount of drug given in the bolus")


class PK_MonocompartmentClint(PKPDODEModel):
    # https://www.nps.org.au/australian-prescriber/articles/pharmacokinetics-made-easy-9-non-linear-pharmacokinetics
    def F(self, t, y):
        Vmax=self.parameters[0]
        Km=self.parameters[1]
        V=self.parameters[2]
        Clint=Vmax/(Km+y)
        return -Clint/V*y

    def G(self, t, dD):
        V=self.parameters[2]
        return dD/V

    def getResponseDimension(self):
        return 1

    def getStateDimension(self):
        return 1

    def getDescription(self):
        return "Monocompartmental intrinsic clearance model (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Clint/V * C + 1/V * dD/dt and Clint=Vmax/(Km+C)"

    def getEquation(self):
        Vmax=self.parameters[0]
        Km=self.parameters[1]
        V=self.parameters[2]
        return "dC/dt = -Clint/(%f) * C + 1/(%f) dD/dt and Clint=(%f)/((%f)+C)"%(V,V,Vmax,Km)

    def getParameterNames(self):
        return ['Vmax','Km','V']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName)
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        Vmaxunits = multiplyUnits(Clunits,yunits)
        self.parameterUnits = [Vmaxunits,yunits,Vunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Vmax
        VmaxLower = lowerBound[0]
        VmaxlUpper = upperBound[0]
        if VmaxLower<0 and VmaxlUpper>0:
            retval.append("Suspicious, Vmax looks like a constant")
        elif VmaxLower<0:
            retval.append("Suspicious, Vmax may be unstable")
        elif VmaxLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Km
        KmLower = lowerBound[1]
        KmUpper = upperBound[1]
        if KmLower<0 and KmUpper>0:
            retval.append("Suspicious, Km looks like a constant")
        elif KmLower<0:
            retval.append("Suspicious, Km may be unstable")
        elif KmLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # V
        VLower = lowerBound[2]
        VUpper = upperBound[2]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, V looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, V seems to be negative")
        else:
            retval.append("True")
        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:1]<0)==0

class PK_Twocompartments(PKPDODEModel):
    def F(self, t, y):
        Cl=self.parameters[0]
        V=self.parameters[1]
        Clp=self.parameters[2]
        Vp=self.parameters[3]
        C=y[0]
        Cp=y[1]

        Q12 = Clp * (C-Cp)
        return np.array([-(Cl*C + Q12)/V, Q12/Vp],np.double)

    def G(self, t, dD):
        V=self.parameters[1]
        return np.array([dD/V,0.0],np.double)

    def getResponseDimension(self):
        return 1

    def getStateDimension(self):
        return 2

    def getDescription(self):
        return "Two-compartments model (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Cl/V * C - Clp/V * (C-Cp) + 1/V * dD/dt and dCp/dt = Clp/Vp * (C-Cp)"

    def getEquation(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        Clp=self.parameters[2]
        Vp=self.parameters[3]
        return "dC/dt = -(%f)/(%f) * C - (%f)/(%f) * (C-Cp) + 1/(%f) * dD/dt; dCp/dt = (%f)/(%f) * (C-Cp)"%(Cl,V,Clp,V,V,Clp,Vp)

    def getParameterNames(self):
        return ['Cl','V','Clp','Vp']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName)
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [Clunits,Vunits,Clunits,Vunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Cl
        ClLower = lowerBound[0]
        ClUpper = upperBound[0]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Cl looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Cl may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # V
        VLower = lowerBound[1]
        VUpper = upperBound[1]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, V looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, V seems to be negative")
        else:
            retval.append("True")

        # Clp
        ClLower = lowerBound[2]
        ClUpper = upperBound[2]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Clp looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Clp may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Vp
        VLower = lowerBound[3]
        VUpper = upperBound[3]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, Vp looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, Vp seems to be negative")
        else:
            retval.append("True")

        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:1]<0)==0

    def printOtherParameterization(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        Clp=self.parameters[2]
        Vp=self.parameters[3]
        print(" ")
        print("Model parameterization ===========")
        print("The two compartment model (excluding the source) has been estimated according to the homogeneous equations:")
        print("This is called the physiological parameterization")
        print("dC/dt = -Cl/V * C - Clp/V * (C-Cp) and dCp/dt = Clp/Vp * (C-Cp)")
        print("where:")
        print("C is the concentration in the central compartment")
        print("Cp is the concentration in the peripheral compartment")
        print("with parameters:")
        print("Cl=%f (%s) (Clearance from central compartment)"%(Cl,PKPDUnit.codeToString(self.parameterUnits[0])))
        print("V=%f (%s) (Apparent volume of distribution in the central compartment)"%(V,PKPDUnit.codeToString(self.parameterUnits[1])))
        print("Clp=%f (%s) (Clearance from the central to the peripheral compartment and viceversa)"%(Clp,PKPDUnit.codeToString(self.parameterUnits[2])))
        print("Vp=%f (%s) (Apparent volume of distribution in the peripheral compartment)"%(Vp,PKPDUnit.codeToString(self.parameterUnits[3])))
        print(" ")
        print("An alternative formulation is given by its microconstants formulation:")
        print("dA/dt = -k10 * A - k12 * A + k21 * Ap and dAp/dt = k12 * A - k21 * Ap")
        print("where:")
        print("A is the amount of drug in the central compartment (A=C*V)")
        print("Ap is the amount of drug in the peripheral compartment (Ap=Cp*Vp)")
        print("with parameters:")
        k12=Clp/V
        k21=Clp/Vp
        k10=Cl/V
        kunits = divideUnits(self.parameterUnits[0],self.parameterUnits[1])
        print("k10=%f (%s) (Elimination from central compartment, k10=Cl/V)"%(k10,PKPDUnit.codeToString(kunits)))
        print("k12=%f (%s) (Transfer from central to peripheral compartment, k12=Clp/V)"%(k12,PKPDUnit.codeToString(kunits)))
        print("k21=%f (%s) (Transfer from peripheral to central compartment, k21=Clp/Vp)"%(k21,PKPDUnit.codeToString(kunits)))
        print(" ")
        print("An alternative formulation is given by its impulse response (response to an intravenous bolus of amount A0):")
        print("A=A0*a*exp(-alpha*t)+A0*b*exp(-beta*t)")
        print("where:")
        print("A is the amount of drug in the central compartment")
        print("A0 is the amount of drug given in the bolus")
        print("The first exponential (alpha) should extinguish faster than the second (beta)")
        print("with parameters:")
        u=k21+k12+k10
        v=k21*k10
        k2units = multiplyUnits(kunits,kunits)
        print("u=%f (%s) (auxiliary variable, u=k12+k21+k10)"%(u,PKPDUnit.codeToString(kunits)))
        print("v=%f (%s) (auxiliary variable, u=k21*k10)"%(v,PKPDUnit.codeToString(k2units)))
        print("alpha and beta are solved from:")
        print("alpha+beta=u")
        print("alpha*beta=v")
        alpha=0.5*(u+math.sqrt(u*u-4*v))
        beta= 0.5*(u-math.sqrt(u*u-4*v))
        a=(alpha-k21)/(alpha-beta)
        b=(k21-beta) /(alpha-beta)
        print("a=%f (%s)     (Amplitude of first exponential,  a=(alpha-k21)/(alpha-beta))"%(a,PKPDUnit.codeToString(PKPDUnit.UNIT_NONE)))
        print("alpha=%f (%s) (Decay of first exponential)"%(alpha,PKPDUnit.codeToString(PKPDUnit.UNIT_INVTIME_MIN)))
        print("b=%f (%s)     (Amplitude of second exponential, b=(k21-beta)/(alpha-beta))"%(b,PKPDUnit.codeToString(PKPDUnit.UNIT_NONE)))
        print("beta=%f (%s)  (Decay of second exponential)"%(beta,PKPDUnit.codeToString(PKPDUnit.UNIT_INVTIME_MIN)))

class PK_TwocompartmentsClint(PKPDODEModel):
    def F(self, t, y):
        Vmax=self.parameters[0]
        Km=self.parameters[1]
        V=self.parameters[2]
        Clp=self.parameters[3]
        Vp=self.parameters[4]
        C=y[0]
        Cp=y[1]

        Clint=Vmax/(Km+C)
        Q12 = Clp * (C-Cp)
        return np.array([-(Clint*C + Q12)/V, Q12/Vp],np.double)

    def G(self, t, dD):
        V=self.parameters[2]
        return np.array([dD/V,0.0],np.double)

    def getResponseDimension(self):
        return 1

    def getStateDimension(self):
        return 2

    def getDescription(self):
        return "Two-compartments intrinsic clearance model (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Clint/V * C - Clp/V * (C-Cp) + 1/V * dD/dt and dCp/dt = Clp/Vp * (C-Cp) and Clint=Vmax/(Km+C)"

    def getEquation(self):
        Vmax=self.parameters[0]
        Km=self.parameters[1]
        V=self.parameters[2]
        Clp=self.parameters[3]
        Vp=self.parameters[4]
        return "dC/dt = -Clint/(%f) * C - (%f)/(%f) * (C-Cp) + 1/(%f) * dD/dt; dCp/dt = (%f)/(%f) * (C-Cp); Clint=(%f)/((%f)+C)"%\
               (V,Clp,V,V,Clp,Vp,Vmax,Km)

    def getParameterNames(self):
        return ['Vmax','Km','V','Clp','Vp']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName)
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        Vmaxunits = multiplyUnits(Clunits,yunits)
        self.parameterUnits = [Vmaxunits,yunits,Vunits,Clunits,Vunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Vmax
        VmaxLower = lowerBound[0]
        VmaxlUpper = upperBound[0]
        if VmaxLower<0 and VmaxlUpper>0:
            retval.append("Suspicious, Vmax looks like a constant")
        elif VmaxLower<0:
            retval.append("Suspicious, Vmax may be unstable")
        elif VmaxLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Km
        KmLower = lowerBound[1]
        KmUpper = upperBound[1]
        if KmLower<0 and KmUpper>0:
            retval.append("Suspicious, Km looks like a constant")
        elif KmLower<0:
            retval.append("Suspicious, Km may be unstable")
        elif KmLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # V
        VLower = lowerBound[2]
        VUpper = upperBound[2]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, V looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, V seems to be negative")
        else:
            retval.append("True")

        # Clp
        ClLower = lowerBound[3]
        ClUpper = upperBound[3]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Clp looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Clp may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Vp
        VLower = lowerBound[4]
        VUpper = upperBound[4]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, Vp looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, Vp seems to be negative")
        else:
            retval.append("True")

        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:1]<0)==0

class PK_TwocompartmentsClintMetabolite(PKPDODEModel):
    def F(self, t, y):
        Vmax=self.parameters[0]
        Km=self.parameters[1]
        V=self.parameters[2]
        Clp=self.parameters[3]
        Vp=self.parameters[4]
        Clm=self.parameters[5]
        Vm=self.parameters[6]
        C=y[0]
        Cm=y[1]
        Cp=y[2]

        Clint=Vmax/(Km+C)
        Q12 = Clp * (C-Cp)
        return np.array([-(Clint*C + Q12)/V, (Clint*C-Clm*Cm)/Vm, Q12/Vp],np.double)

    def G(self, t, dD):
        V=self.parameters[2]
        return np.array([dD/V,0.0,0.0],np.double)

    def getResponseDimension(self):
        return 2

    def getStateDimension(self):
        return 3

    def getDescription(self):
        return "Two-compartments intrinsic clearance model with metabolites (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Clint/V * C - Clp/V * (C-Cp) + 1/V * dD/dt and dCp/dt = Clp/Vp * (C-Cp) and Clint=Vmax/(Km+C) and dCm/dt=Clint*C/Vm-Clm*Cm/Vm"

    def getEquation(self):
        Vmax=self.parameters[0]
        Km=self.parameters[1]
        V=self.parameters[2]
        Clp=self.parameters[3]
        Vp=self.parameters[4]
        Clm=self.parameters[5]
        Vm=self.parameters[6]
        return "dC/dt = -Clint/(%f) * C - (%f)/(%f) * (C-Cp) + 1/(%f) * dD/dt; dCp/dt = (%f)/(%f) * (C-Cp); Clint=(%f)/((%f)+C); dCm/dt=Clint*C/(%f)-(%f)*Cm/(%f)"%\
               (V,Clp,V,V,Clp,Vp,Vmax,Km,Vm,Clm,Vm)

    def getParameterNames(self):
        return ['Vmax','Km','V','Clp','Vp','Clm','Vm']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName[0])
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        Vmaxunits = multiplyUnits(Clunits,yunits)
        self.parameterUnits = [Vmaxunits,yunits,Vunits,Clunits,Vunits,Clunits,Vunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Vmax
        VmaxLower = lowerBound[0]
        VmaxlUpper = upperBound[0]
        if VmaxLower<0 and VmaxlUpper>0:
            retval.append("Suspicious, Vmax looks like a constant")
        elif VmaxLower<0:
            retval.append("Suspicious, Vmax may be unstable")
        elif VmaxLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Km
        KmLower = lowerBound[1]
        KmUpper = upperBound[1]
        if KmLower<0 and KmUpper>0:
            retval.append("Suspicious, Km looks like a constant")
        elif KmLower<0:
            retval.append("Suspicious, Km may be unstable")
        elif KmLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # V
        VLower = lowerBound[2]
        VUpper = upperBound[2]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, V looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, V seems to be negative")
        else:
            retval.append("True")

        # Clp
        ClLower = lowerBound[3]
        ClUpper = upperBound[3]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Clp looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Clp may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Vp
        VLower = lowerBound[4]
        VUpper = upperBound[4]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, Vp looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, Vp seems to be negative")
        else:
            retval.append("True")

        # Clm
        ClLower = lowerBound[5]
        ClUpper = upperBound[5]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Clm looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Clp may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Vm
        VLower = lowerBound[6]
        VUpper = upperBound[6]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, Vm looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, Vm seems to be negative")
        else:
            retval.append("True")

        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:1]<0)==0

class PK_TwocompartmentsAutoinduction(PKPDODEModel):
    def F(self, t, y):
        E0=self.parameters[0]
        a=self.parameters[1]
        kout=self.parameters[2]
        V=self.parameters[3]
        Clp=self.parameters[4]
        Vp=self.parameters[5]
        C=y[0]
        Cp=y[1]
        E=y[2]

        Cl=a*(1+E)
        Q12 = Clp * (C-Cp)
        return np.array([-(Cl*C + Q12)/V, Q12/Vp, kout*(E0+C-E)],np.double)

    def G(self, t, dD):
        V=self.parameters[3]
        return np.array([dD/V,0.0,0.0],np.double)

    def getResponseDimension(self):
        return 1

    def getStateDimension(self):
        return 3

    def getDescription(self):
        return "Two-compartments model with autoinduction (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Cl/V * C - Clp/V * (C-Cp) + 1/V * dD/dt, dCp/dt = Clp/Vp * (C-Cp), dE/dt=kout*(E0+C-E) and Cl=a*(1+E)"

    def getEquation(self):
        E0=self.parameters[0]
        a=self.parameters[1]
        kout=self.parameters[2]
        V=self.parameters[3]
        Clp=self.parameters[4]
        Vp=self.parameters[5]
        return "dC/dt = -Cl/(%f) * C - (%f)/(%f) * (C-Cp) + 1/(%f) * dD/dt; dCp/dt = (%f)/(%f) * (C-Cp); dE/dt=(%f)*(%f+C-E) and Cl=(%f)*(1+E)"%(V,Clp,V,V,Clp,Vp,kout,E0,a)

    def getParameterNames(self):
        return ['E0', 'a', 'kout', 'V','Clp','Vp']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName)
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [PKPDUnit.UNIT_NONE,Clunits,PKPDUnit.UNIT_NONE,Vunits,Clunits,Vunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # E0
        E0Lower = lowerBound[0]
        E0Upper = upperBound[0]
        if E0Lower<0 and E0Upper>0:
            retval.append("Suspicious, Cl looks like a constant")
        elif E0Lower<0:
            retval.append("Suspicious, Cl may be unstable")
        elif E0Lower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # A
        ALower = lowerBound[1]
        AUpper = upperBound[1]
        if ALower<0 and AUpper>0:
            retval.append("Suspicious, A looks like a constant")
        elif ALower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # kout
        koutLower = lowerBound[2]
        koutUpper = upperBound[2]
        if koutLower<0 and koutUpper>0:
            retval.append("Suspicious, A looks like a constant")
        elif koutLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # V
        VLower = lowerBound[3]
        VUpper = upperBound[3]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, V looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, V seems to be negative")
        else:
            retval.append("True")

        # Clp
        ClLower = lowerBound[4]
        ClUpper = upperBound[4]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Clp looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Clp may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Vp
        VLower = lowerBound[5]
        VUpper = upperBound[5]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, Vp looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, Vp seems to be negative")
        else:
            retval.append("True")

        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:1]<0)==0


class PK_MonocompartmentUrine(PKPDODEModel):
    def F(self, t, y):
        C = y[0]
        Cl=self.parameters[0]
        V=self.parameters[1]
        fe=self.parameters[2]
        return np.array([-Cl/V*C, fe*Cl*C],np.double)

    def G(self, t, dD):
        V=self.parameters[1]
        return np.array([dD/V,0.0],np.double)

    def getResponseDimension(self):
        return 2

    def getStateDimension(self):
        return 2

    def getDescription(self):
        return "Monocompartmental model urine (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Cl/V * C + 1/V * dD/dt and dAu/dt=fe*Cl*C"

    def getEquation(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        fe=self.parameters[2]
        return "dC/dt = -(%f)/(%f) * C + 1/(%f) dD/dt and dAu/dt = (%f)*(%f)*C"%(Cl,V,V,fe,Cl)

    def getParameterNames(self):
        return ['Cl','V','fe']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName[0])
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [Clunits,Vunits,PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Cl
        ClLower = lowerBound[0]
        ClUpper = upperBound[0]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Ka looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Ka may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # V
        VLower = lowerBound[1]
        VUpper = upperBound[1]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, V looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, V seems to be negative")
        else:
            retval.append("True")

        feLower = lowerBound[2]
        feUpper = upperBound[2]
        if feLower<0 and feUpper>0:
            retval.append("Suspicious, fe looks like 0")
        elif feUpper<0:
            retval.append("Suspicious, fe seems to be negative")
        elif feLower>1:
            retval.append("Suspicious, fe seems to be larger than 1")
        else:
            retval.append("True")
        return retval

    def areParametersValid(self, p):
        return np.sum(p<0)==0 and p[2]<1

class PK_MonocompartmentLinkPD(PKPDODEModel):
    def F(self, t, y):
        C = y[0]
        Cb = y[2]

        Cl=self.parameters[0]
        V=self.parameters[1]
        Clb=self.parameters[2]
        Vb=self.parameters[3]
        Q12 = Clb * (C-Cb)

        return np.array([-(Cl*C + Q12)/V, 0.0, Q12/Vb],np.double)

    def G(self, t, dD):
        V=self.parameters[1]
        return np.array([dD/V,0.0,0.0],np.double)

    def H(self, y):
        Cb = y[2]

        E0=self.parameters[4]
        a=self.parameters[5]
        b=self.parameters[6]
        Cbm=self.parameters[7]
        Cbb=math.pow(Cb,b)

        y[1]=E0*(1+a*Cbb/(math.pow(Cbm,b)+Cbb))

    def getResponseDimension(self):
        return 2

    def getStateDimension(self):
        return 3

    def getDescription(self):
        return "One-compartmental model link pd (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Cl/V * C - Clb/V * (C-Cb) + 1/V * dD/dt, dCb/dt = Clb/Vb * (C-Cb) and E=E0+a*Cb^b/(Cbm^b+Cb^b)"

    def getEquation(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        Clp=self.parameters[2]
        Vp=self.parameters[3]
        E0=self.parameters[4]
        a=self.parameters[5]
        b=self.parameters[6]
        Cbm=self.parameters[7]
        return "dC/dt = -(%f)/(%f) * C - (%f)/(%f) * (C-Cb) + 1/(%f) * dD/dt, dCb/dt = (%f)/(%f) * (C-Cb) and E = (%f)+(%f)*Cb^(%f)/((%f)^(%f)+Cb^(%f))"%\
               (Cl,V,Clp,V,V,Clp,Vp,E0,a,b,Cbm,b,b)

    def getParameterNames(self):
        return ['Cl','V','Clp','Vp','E0','a','b','Cbm']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName[0])
        Eunits = self.experiment.getVarUnits(self.yName[1])
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [Clunits,Vunits,Clunits,Vunits,Eunits,PKPDUnit.UNIT_NONE,PKPDUnit.UNIT_NONE,yunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        retval.append(signifcantMessage1("Cl", lowerBound[0],upperBound[0]))
        retval.append(signifcantMessage2("V",  lowerBound[1],upperBound[1]))
        retval.append(signifcantMessage1("Clp",lowerBound[2],upperBound[2]))
        retval.append(signifcantMessage2("Vp", lowerBound[3],upperBound[3]))
        retval.append(signifcantMessage2("E0", lowerBound[4],upperBound[4]))
        retval.append(signifcantMessage2("a",  lowerBound[5],upperBound[5]))

        if lowerBound[6] < 1 and upperBound[6] > 1:
            retval.append("Suspicious, %s looks like 1" % "b")
        else:
            retval.append("True")

        retval.append(signifcantMessage2("Cbm",lowerBound[7],upperBound[7]))
        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:4]<0)==0

class PK_MonocompartmentPD(PKPDODEModel):
    def F(self, t, y):
        C = y[0]

        Cl=self.parameters[0]
        V=self.parameters[1]

        return np.array([-Cl*C/V, 0.0],np.double)

    def G(self, t, dD):
        V=self.parameters[1]
        return np.array([dD/V,0.0],np.double)

    def H(self, y):
        C = y[0]

        E0=self.parameters[2]
        a=self.parameters[3]
        b=self.parameters[4]
        Cm=self.parameters[5]
        try:
            Cb=math.pow(C,b)
            Cmb=math.pow(Cm,b)
            y[1] = E0*(1+a*Cb/(Cmb+Cb))
        except:
            y[1] = E0

    def getResponseDimension(self):
        return 2

    def getStateDimension(self):
        return 2

    def getDescription(self):
        return "One-compartmental model pd (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Cl/V * C + 1/V * dD/dt and E=E0+a*C^b/(Cm^b+C^b)"

    def getEquation(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        E0=self.parameters[2]
        a=self.parameters[3]
        b=self.parameters[4]
        Cm=self.parameters[5]
        return "dC/dt = -(%f)/(%f) * C + 1/(%f) * dD/dt, and E=(%f)+(%f)*Cb^(%f)/((%f)^(%f)+Cb^(%f))"%\
               (Cl,V,V,E0,a,b,Cm,b,b)

    def getParameterNames(self):
        return ['Cl','V','E0','a','b','Cbm']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName[0])
        Eunits = self.experiment.getVarUnits(self.yName[1])
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [Clunits,Vunits,Eunits,PKPDUnit.UNIT_NONE,PKPDUnit.UNIT_NONE,yunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        retval.append(signifcantMessage1("Cl", lowerBound[0],upperBound[0]))
        retval.append(signifcantMessage2("V",  lowerBound[1],upperBound[1]))
        retval.append(signifcantMessage2("E0", lowerBound[2],upperBound[2]))
        retval.append(signifcantMessage2("a",  lowerBound[3],upperBound[3]))

        if lowerBound[4] < 1 and upperBound[4] > 1:
            retval.append("Suspicious, %s looks like 1" % "b")
        else:
            retval.append("True")

        retval.append(signifcantMessage2("Cm",lowerBound[5],upperBound[5]))
        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:2]<0)==0

class PK_TwocompartmentsUrine(PKPDODEModel):
    def F(self, t, y):
        C = y[0]
        Cp = y[2]

        Cl=self.parameters[0]
        V=self.parameters[1]
        Clp=self.parameters[2]
        Vp=self.parameters[3]
        fe=self.parameters[4]
        Q12 = Clp * (C-Cp)

        return np.array([-(Cl*C + Q12)/V, fe*Cl*C, Q12/Vp],np.double)

    def G(self, t, dD):
        V=self.parameters[1]
        return np.array([dD/V,0.0,0.0],np.double)

    def getResponseDimension(self):
        return 2

    def getStateDimension(self):
        return 3

    def getDescription(self):
        return "Two-compartmental model urine (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Cl/V * C - Clp/V * (C-Cp) + 1/V * dD/dt, dCp/dt = Clp/Vp * (C-Cp) and dAu/dt=fe*Cl*C"

    def getEquation(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        Clp=self.parameters[2]
        Vp=self.parameters[3]
        fe=self.parameters[4]
        return "dC/dt = -(%f)/(%f) * C - (%f)/(%f) * (C-Cp) + 1/(%f) * dD/dt, dCp/dt = (%f)/(%f) * (C-Cp) and dAu/dt = (%f)*(%f)*C"%\
               (Cl,V,Clp,V,V,Clp,Vp,fe,Cl)

    def getParameterNames(self):
        return ['Cl','V','Clp','Vp','fe']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName[0])
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [Clunits,Vunits,Clunits,Vunits,PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        # Cl
        ClLower = lowerBound[0]
        ClUpper = upperBound[0]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Ka looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Ka may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # V
        VLower = lowerBound[1]
        VUpper = upperBound[1]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, V looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, V seems to be negative")
        else:
            retval.append("True")

        # Clp
        ClLower = lowerBound[2]
        ClUpper = upperBound[2]
        if ClLower<0 and ClUpper>0:
            retval.append("Suspicious, Clp looks like a constant")
        elif ClLower<0:
            retval.append("Suspicious, Clp may be unstable")
        elif ClLower>0:
            retval.append("True")
        else:
            retval.append("NA")

        # Vp
        VLower = lowerBound[3]
        VUpper = upperBound[3]
        if VLower<0 and VUpper>0:
            retval.append("Suspicious, Vp looks like 0")
        elif VUpper<0:
            retval.append("Suspicious, Vp seems to be negative")
        else:
            retval.append("True")

        # fe
        feLower = lowerBound[4]
        feUpper = upperBound[4]
        if feLower<0 and feUpper>0:
            retval.append("Suspicious, fe looks like 0")
        elif feUpper<0:
            retval.append("Suspicious, fe seems to be negative")
        elif feLower>1:
            retval.append("Suspicious, fe seems to be larger than 1")
        else:
            retval.append("True")
        return retval

    def areParametersValid(self, p):
        return np.sum(p<0)==0 and p[4]<1

class PK_TwocompartmentsBoth(PK_Twocompartments):
    def getResponseDimension(self):
        return 2

    def getDescription(self):
        return "Two-compartmental model both (%s)"%self.__class__.__name__

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName[0])
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [Clunits,Vunits,Clunits,Vunits]
        return self.parameterUnits

class PK_TwocompartmentsBothPD(PKPDODEModel):
    def F(self, t, y):
        Cl=self.parameters[0]
        V=self.parameters[1]
        Clp=self.parameters[2]
        Vp=self.parameters[3]
        C=y[0]
        Cp=y[1]

        Q12 = Clp * (C-Cp)
        return np.array([-(Cl*C + Q12)/V, Q12/Vp, 0.0],np.double)

    def G(self, t, dD):
        V=self.parameters[1]
        return np.array([dD/V,0.0,0.0],np.double)

    def H(self, y):
        Cp = y[1]

        E0=self.parameters[4]
        a=self.parameters[5]
        b=self.parameters[6]
        Cpm=self.parameters[7]
        try:
            Cpb=math.pow(Cp,b)
            Cpmb=math.pow(Cpm,b)
            y[2] = E0*(1+a*Cpb/(math.pow(Cpm,b)+Cpb))
        except:
            y[2] = E0

    def getResponseDimension(self):
        return 3

    def getStateDimension(self):
        return 3

    def getDescription(self):
        return "Two-compartments both PD model (%s)"%self.__class__.__name__

    def getModelEquation(self):
        return "dC/dt = -Cl/V * C - Clp/V * (C-Cp) + 1/V * dD/dt, dCp/dt = Clp/Vp * (C-Cp) and E=E0+a*Cp^b/(Cpm^b+Cp^b)"

    def getEquation(self):
        Cl=self.parameters[0]
        V=self.parameters[1]
        Clp=self.parameters[2]
        Vp=self.parameters[3]
        return "dC/dt = -(%f)/(%f) * C - (%f)/(%f) * (C-Cp) + 1/(%f) * dD/dt; dCp/dt = (%f)/(%f) * (C-Cp)"%(Cl,V,Clp,V,V,Clp,Vp)

    def getParameterNames(self):
        return ['Cl','V','Clp','Vp','E0','a','b','Cpm']

    def calculateParameterUnits(self,sample):
        xunits = unitFromString("min")
        yunits = self.experiment.getVarUnits(self.yName[0])
        Eunits = self.experiment.getVarUnits(self.yName[1])
        Vunits = divideUnits(self.Dunits,yunits)
        Clunits = divideUnits(Vunits,xunits)
        self.parameterUnits = [Clunits,Vunits,Clunits,Vunits,Eunits,PKPDUnit.UNIT_NONE,PKPDUnit.UNIT_NONE,yunits]
        return self.parameterUnits

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        retval.append(signifcantMessage1("Cl", lowerBound[0],upperBound[0]))
        retval.append(signifcantMessage2("V",  lowerBound[1],upperBound[1]))
        retval.append(signifcantMessage1("Clp",lowerBound[2],upperBound[2]))
        retval.append(signifcantMessage2("Vp", lowerBound[3],upperBound[3]))
        retval.append(signifcantMessage2("E0", lowerBound[4],upperBound[4]))
        retval.append(signifcantMessage2("a",  lowerBound[5],upperBound[5]))

        if lowerBound[6] < 1 and upperBound[6] > 1:
            retval.append("Suspicious, %s looks like 1" % "b")
        else:
            retval.append("True")

        retval.append(signifcantMessage2("Cpm",lowerBound[7],upperBound[7]))
        return retval

    def areParametersValid(self, p):
        return np.sum(p[0:4]<0)==0

