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
Biopharmaceutics: Drug sources and how they dissolve
"""
import copy
from itertools import izip
import math
import numpy as np
from pyworkflow.em.pkpd_units import PKPDUnit, changeRateToWeight

class BiopharmaceuticsModel:
    def __init__(self):
        self.parameters = []

    def getNumberOfParameters(self):
        return len(self.getParameterNames())

    def getDescription(self):
        pass

    def getParameterNames(self):
        pass

    def calculateParameterUnits(self,sample):
        pass

    def setParameters(self, parameters):
        self.parameters = parameters

    def getAg(self,t):
        # Total amount of drug that is available at time t
        return 0.0

    def getEquation(self):
        return ""

    def getModelEquation(self):
        return ""

    def getDescription(self):
        return ""

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        for i in range(len(self.parameters)):
            lower = lowerBound[i]
            upper = upperBound[i]
            if lower<0 and upper>0:
                retval.append("False")
            elif lower>0 or upper<0:
                retval.append("True")
            else:
                retval.append("NA")
        return retval

    def areParametersValid(self, p):
        return np.sum(p<0)==0

class BiopharmaceuticsModelOrder0(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Constant absorption rate']

    def getParameterNames(self):
        return ['Rin']

    def calculateParameterUnits(self,sample):
        self.parameterUnits = [PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN]
        return self.parameterUnits

    def getAg(self,t):
        if t<0:
            return 0.0
        Rin = self.parameters[0]
        return max(self.Amax-Rin*t,0.0)

    def getEquation(self):
        Rin = self.parameters[0]
        return "D(t)=(%f)*t"%(Rin)

    def getModelEquation(self):
        return "D(t)=Rin*t"

    def getDescription(self):
        return "Zero order absorption (%s)"%self.__class__.__name__

class BiopharmaceuticsModelOrder01(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Constant absorption rate','Constant absorption time','Absorption rate']

    def getParameterNames(self):
        return ['Rin','t0','Ka']

    def calculateParameterUnits(self,sample):
        self.parameterUnits = [PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN, PKPDUnit.UNIT_TIME_MIN, PKPDUnit.UNIT_INVTIME_MIN]
        return self.parameterUnits

    def getAg(self,t):
        if t<0:
            return 0.0
        Rin = self.parameters[0]
        t0 = self.parameters[1]
        Ka = self.parameters[2]
        if t<t0:
            return max(self.Amax-Rin*t,0.0)
        else:
            A0=max(self.Amax-Rin*t0,0.0)
            return A0*math.exp(-Ka*(t-t0))

    def getEquation(self):
        Rin = self.parameters[0]
        t0 = self.parameters[1]
        Ka = self.parameters[2]
        return "D(t)=(%f)*t if t<%f and (%f)*t+((%f)-(%f)*(%f))*(1-exp(-(%f)*t) if t>%f"%(Rin, t0, Rin, self.Amax, Rin, t0, Ka, t0)

    def getModelEquation(self):
        return "D(t)=Rin*t if t<t0 and Rin*t0+(Amax-Rin*t0)*(1-exp(-Ka*(t-t0)) if t>t0"

    def getDescription(self):
        return "Zero-First Mixed order absorption (%s)"%self.__class__.__name__

class BiopharmaceuticsModelOrder1(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Absorption rate']

    def getParameterNames(self):
        return ['Ka']

    def calculateParameterUnits(self,sample):
        self.parameterUnits = [PKPDUnit.UNIT_INVTIME_MIN]
        return self.parameterUnits

    def getAg(self,t):
        if t<0:
            return 0.0
        Ka = self.parameters[0]
        return self.Amax*math.exp(-Ka*t)

    def getEquation(self):
        Ka = self.parameters[0]
        return "D(t)=(%f)*(1-exp(-(%f)*t)"%(self.Amax,Ka)

    def getModelEquation(self):
        return "D(t)=Amax*(1-exp(-Ka*t))"

    def getDescription(self):
        return "First order absorption (%s)"%self.__class__.__name__

class BiopharmaceuticsModelOrderFractional(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Initial amount','Constant absorption rate', 'alpha']

    def getParameterNames(self):
        return ['Amax','K','alpha']

    def calculateParameterUnits(self,sample):
        self.parameterUnits = [PKPDUnit.UNIT_WEIGHT_mg,PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN,PKPDUnit.UNIT_NONE]
        return self.paramterUnits

    def getAg(self,t):
        # COSS Hay que pensar si esto es correcto
        if t<=0:
            return 0.0
        Amax = self.parameters[0]
        K = self.parameters[1]
        alpha = self.parameters[2]
        aux = alpha*K*t
        if aux>Amax:
            return Amax
        return Amax-math.pow(math.pow(Amax,alpha)-aux,1.0/alpha)

    def getEquation(self):
        Amax = self.parameters[0]
        K = self.parameters[1]
        alpha = self.parameters[2]
        return "D(t)=(%f)-((%f)^(%f)-(%f)*(%f)*t)^(1/(%f))"%(Amax,Amax,alpha,K,alpha,alpha)

    def getModelEquation(self):
        return "D(t)=Amax-(Amax^alpha-alpha*K*t)^(1/alpha)"

    def getDescription(self):
        return "Fractional order absorption (%s)"%self.__class__.__name__

    def areParametersValid(self, p):
        return np.sum(p<0)==0 and p[2]>0 and p[2]<1

class BiopharmaceuticsModelImmediateAndOrder1(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Absorption rate','Immediate fraction']

    def getParameterNames(self):
        return ['Ka','F']

    def calculateParameterUnits(self,sample):
        self.parameterUnits = [PKPDUnit.UNIT_INVTIME_MIN,PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def getAg(self,t):
        if t<0:
            return 0.0
        Ka = self.parameters[0]
        F = self.parameters[1]
        return self.Amax*(1-F)*math.exp(-Ka*t)

    def getEquation(self):
        Ka = self.parameters[0]
        F = self.parameters[1]
        return "D(t)=(%f)*delta(t)+(%f)*(1-exp(-(%f)*t)"%(self.Amax*F,self.Amax*(1-F),Ka)

    def getModelEquation(self):
        return "D(t)=Amax*F*delta(t)+(1-F)*Amax*(1-exp(-Ka*t))"

    def getDescription(self):
        return "Immediate and First order absorption (%s)"%self.__class__.__name__

class BiopharmaceuticsModelOrder1AndOrder1(BiopharmaceuticsModel):
    def getDescription(self):
        return ['Absorption rate1','Absorption rate2','tlag2','Fraction 1']

    def getParameterNames(self):
        return ['Ka1','Ka2','tlag12','F1']

    def calculateParameterUnits(self,sample):
        self.parameterUnits = [PKPDUnit.UNIT_INVTIME_MIN,PKPDUnit.UNIT_INVTIME_MIN,PKPDUnit.UNIT_TIME_MIN,PKPDUnit.UNIT_NONE]
        return self.parameterUnits

    def getAg(self,t):
        if t<0:
            return 0.0
        Ka1 = self.parameters[0]
        Ka2 = self.parameters[1]
        tlag12 = self.parameters[2]
        F1 = self.parameters[3]
        A1=F1*math.exp(-Ka1*t)
        A2=1-F1
        if t>tlag12:
            A2=(1-F1)*math.exp(-Ka2*(t-tlag12))
        return self.Amax*(A1+A2)

    def getEquation(self):
        Ka1 = self.parameters[0]
        Ka2 = self.parameters[1]
        tlag12 = self.parameters[2]
        F1 = self.parameters[3]
        return "D(t)=(%f)*(1-exp(-(%f)*t)+(%f)*(1-exp(-(%f)*(t-%f))"%(self.Amax*F1,Ka1,self.Amax*(1-F1),Ka2,tlag12)

    def getModelEquation(self):
        return "D(t)=Amax*F1*(1-exp(-Ka1*t))+(1-F1)*Amax*(1-exp(-Ka2*(t-tlag12)))"

    def getDescription(self):
        return "First and First order absorption (%s)"%self.__class__.__name__

class PKPDVia:
    def __init__(self):
        self.viaName = None
        self.via = None
        self.viaProfile = None
        self.tlag = 0
        self.tunits = None
        self.bioavailability = 1
        self.paramsToOptimize = []
        self.paramsUnitsToOptimize=[]

    def parseTokens(self,tokens):
        # Intravenous; iv; [tlag=0 min]; [bioavailability=1]
        # Oral; ev1; tlag=0 min; bioavailability=1
        # Default values
        self.tlag = 0
        self.bioavailability = 1
        self.tunits = PKPDUnit("min")

        # Get name
        currentToken = 0
        self.viaName = tokens[currentToken].strip()
        currentToken+=1

        # Get via
        self.via = tokens[currentToken].strip()
        currentToken+=1

        while currentToken<len(tokens):
            optionalTokens=tokens[currentToken].strip()
            if '=' in optionalTokens:
                optionalTokens=optionalTokens.split('=')
                optionalVar=optionalTokens[0].strip()
                if optionalVar=="tlag":
                    optionalTokens=optionalTokens[1].split()
                    self.tlag=float(optionalTokens[0].strip())
                    unitString = optionalTokens[1].strip()
                    self.tunits = PKPDUnit(unitString)
                    if not self.tunits.unit:
                        raise Exception("Unrecognized unit: %s"%unitString)
                    if not self.tunits.isTime():
                        raise Exception("Time unit is not valid")
                elif optionalVar=="bioavailability":
                    self.bioavailability=float(optionalTokens[1].strip())
            else:
                optionalTokens=optionalTokens.split()[0]
                self.paramsToOptimize.append(optionalTokens)
                if optionalTokens=="tlag":
                    self.paramsUnitsToOptimize.append(PKPDUnit.UNIT_TIME_MIN)
                elif optionalTokens=="bioavailability":
                    self.paramsUnitsToOptimize.append(PKPDUnit.UNIT_NONE)
            currentToken+=1

    def _printToStream(self,fh):
        outStr="%s; %s; "%(self.viaName,self.via)
        if "tlag" in self.paramsToOptimize:
            outStr+=" tlag min; "
        else:
            outStr+="tlag=%f %s; "%(self.tlag,self.tunits._toString())
        if "bioavailability" in self.paramsToOptimize:
            outStr+=" bioavailability"
        else:
            outStr+="bioavailability=%f"%self.bioavailability
        fh.write("%s\n"%outStr)

    def prepare(self):
        if self.via=="iv":
            self.viaProfile=None
        else:
            if self.viaProfile==None:
                if self.via=="iv-ev1":
                    self.viaProfile=BiopharmaceuticsModelImmediateAndOrder1()
                elif self.via=="ev0":
                    self.viaProfile=BiopharmaceuticsModelOrder0()
                elif self.via=="ev01":
                    self.viaProfile=BiopharmaceuticsModelOrder01()
                elif self.via=="ev1":
                    self.viaProfile=BiopharmaceuticsModelOrder1()
                elif self.via=="evFractional":
                    self.viaProfile=BiopharmaceuticsModelOrderFractional()
                elif self.via=="ev1-ev1":
                    self.viaProfile=BiopharmaceuticsModelOrder1AndOrder1()

    def changeTimeUnitsToMinutes(self):
        if self.tunits.unit==PKPDUnit.UNIT_TIME_MIN:
            pass
        elif self.tunits.unit==PKPDUnit.UNIT_TIME_H:
            self.tlag *= 60
        else:
            raise Exception("Time for doses must be hours or minutes")

    def getEquation(self):
        if self.via == "iv":
            retval = "D=Div(t)"
        else:
            retval = self.viaProfile.getEquation()
        return "%s (tlag=%f, bioavailability=%f)"%(retval,self.tlag,self.bioavailability)

    def getModelEquation(self):
        if self.via == "iv":
            return "D=Div(t)"
        else:
            return self.viaProfile.getModelEquation()

    def getDescription(self):
        if self.via == "iv":
            return "Intravenous dose"
        else:
            return self.viaProfile.getDescription()

    def getParameterNames(self):
        names=copy.copy(self.paramsToOptimize)
        if self.via != "iv":
            names+=self.viaProfile.getParameterNames()
        return [self.viaName+"_"+x for x in names]

    def getNumberOfParameters(self):
        return len(self.getParameterNames())

    def calculateParameterUnits(self,sample):
        if self.via == "iv":
            return self.paramsUnitsToOptimize
        else:
            return self.paramsUnitsToOptimize+self.viaProfile.calculateParameterUnits(sample)

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        currentIdx=0
        if self.paramsToOptimize:
            for paramName in self.paramsToOptimize:
                if paramName=="tlag":
                    retval+=[str(lowerBound[currentIdx]>0)]
                elif paramName=="bioavailability":
                    retval+=[str(upperBound[currentIdx]<1)]
                currentIdx+=1
        if self.via == "iv":
            return retval
        else:
            return np.concatenate((np.asarray(retval),self.viaProfile.areParametersSignificant(lowerBound[currentIdx:], upperBound[currentIdx:])))

    def areParametersValid(self, p):
        retval=True
        currentIdx=0
        if self.paramsToOptimize:
            for paramName in self.paramsToOptimize:
                if paramName=="tlag":
                    retval=retval and p[currentIdx]>=0
                elif paramName=="bioavailability":
                    retval=retval and p[currentIdx]>=0 and p[currentIdx]<=1
                currentIdx+=1
        if self.via == "iv":
            return retval
        else:
            return retval and self.viaProfile.areParametersValid(p[currentIdx:])

    def setParameters(self, p):
        currentIdx=0
        if self.paramsToOptimize:
            for paramName in self.paramsToOptimize:
                if paramName=="tlag":
                    self.tlag=p[currentIdx]
                elif paramName=="bioavailability":
                    self.bioavailability=p[currentIdx]
                currentIdx+=1
        if self.via != "iv":
            self.viaProfile.setParameters(p[currentIdx:])

def createVia(line):
    via = PKPDVia()
    via.parseTokens(line.split(';'))
    via.prepare()
    return via

class PKPDDose:
    TYPE_BOLUS = 1
    TYPE_REPEATED_BOLUS = 2
    TYPE_INFUSION = 3

    def __init__(self):
        self.doseName = None
        self.via = None
        self.doseType = None
        self.doseAmount = None
        self.t0 = None
        self.tF = None
        self.every = None
        self.tunits = None
        self.dunits = None
        self.paramsToOptimize = []
        self.paramsUnitsToOptimize=[]

    def parseTokens(self,tokens,vias):
        # Dose1; via=Intravenous; bolus; t=0 min; d=60*$(weight)/1000 mg
        # Dose1; via=Oral; repeated_bolus; t=0:8:48 h; d=60*$(weight)/1000 mg
        # Dose1; via=Intravenous; infusion; t=0:59 min; d=1 mg/min

        # Get name
        currentToken = 0
        self.doseName = tokens[currentToken].strip()
        currentToken+=1

        # Get via
        viaName = tokens[currentToken].strip().split('=')[1]
        if viaName in vias:
            self.via=vias[viaName]
        else:
            raise Exception("Unrecognized via %s"%viaName)
        currentToken+=1

        # Get type
        doseTypeString = tokens[currentToken].strip()
        if doseTypeString=="bolus":
            self.doseType = PKPDDose.TYPE_BOLUS
        elif doseTypeString=="repeated_bolus":
            self.doseType = PKPDDose.TYPE_REPEATED_BOLUS
        elif doseTypeString=="infusion":
            self.doseType = PKPDDose.TYPE_INFUSION
        else:
            raise Exception("Unrecognized dose type %s"%doseTypeString)
        currentToken+=1

        # Get time description
        timeUnitsString = tokens[currentToken].strip().lower()
        timeTokens=timeUnitsString.split()
        if len(timeTokens)!=2:
            raise Exception("Time description is badly formed %s"%timeUnitsString)
        timeString = timeTokens[0].strip().split('=')[1]
        unitString = timeTokens[1].strip()
        if doseTypeString=="bolus":
            self.t0 = float(timeString)
        elif doseTypeString=="repeated_bolus":
            timeTokens = timeString.split(":")
            self.t0 = float(timeTokens[0].strip())
            self.every = float(timeTokens[1].strip())
            self.tF = float(timeTokens[2].strip())
        elif doseTypeString=="infusion":
            timeTokens = timeString.split(":")
            self.t0 = float(timeTokens[0].strip())
            self.tF = float(timeTokens[1].strip())
        else:
            raise Exception("Unrecognized dose type %s"%doseTypeString)

        # Get time units
        self.tunits = PKPDUnit(unitString)
        if not self.tunits.unit:
            raise Exception("Unrecognized unit: %s"%unitString)
        if not self.tunits.isTime():
            raise Exception("Time unit is not valid")
        currentToken+=1

        # Get dose units
        doseUnitsString = tokens[currentToken].strip().lower()
        doseTokens=doseUnitsString.split()
        if len(doseTokens)!=2:
            raise Exception("Dose description is badly formed %s"%doseUnitsString)

        self.doseAmount = doseTokens[0].strip().lower().split("=")[1]
        unitString = doseTokens[1].strip()
        self.dunits = PKPDUnit(unitString)
        if not self.dunits.unit:
            raise Exception("Unrecognized unit: %s"%unitString)
        if doseTypeString=="infusion":
            if not self.dunits.isWeightInvTime():
                raise Exception("After normalization, the dose must be a weight/time (=rate)")
        else:
            if not self.dunits.isWeight():
                raise Exception("After normalization, the dose must be a weight")
        currentToken+=1

    def _printToStream(self,fh):
        outStr=self.doseName+"; "+self.getDoseString2()
        fh.write("%s\n"%outStr)

    def getDoseString(self):
        if self.doseType == PKPDDose.TYPE_BOLUS:
            doseString = "bolus; t=%f" % self.t0
        elif self.doseType == PKPDDose.TYPE_REPEATED_BOLUS:
            doseString = "repeated_bolus; t=%f:%f:%f" % (self.t0, self.every, self.tF)
        elif self.doseType == PKPDDose.TYPE_INFUSION:
            doseString = "infusion; t=%f:%f" % (self.t0, self.tF)
        else:
            doseString = ""
        return doseString+" "+self.tunits._toString()

    def getDoseString2(self):
        outStr="via=%s; %s; d=%s %s" % (self.via.viaName,
                                    self.getDoseString(),
                                    self.doseAmount,
                                    self.dunits._toString())
        return outStr

    def changeTimeUnitsToMinutes(self):
        if self.tunits.unit==PKPDUnit.UNIT_TIME_MIN:
            pass
        elif self.tunits.unit==PKPDUnit.UNIT_TIME_H:
            if self.doseType == PKPDDose.TYPE_BOLUS:
                self.t0 *= 60
            elif self.doseType == PKPDDose.TYPE_REPEATED_BOLUS:
                self.t0 *= 60
                self.tF *= 60
                self.every *= 60
            elif self.doseType == PKPDDose.TYPE_INFUSION:
                self.t0 *= 60
                self.tF *= 60
        else:
            raise Exception("Time for doses must be hours or minutes")

    def prepare(self):
        self.via.prepare()

    def getDoseAt(self,t0,dt=0.5):
        """Dose between t0<=t<t0+dt, t0 is in minutes"""
        t0-=self.via.tlag
        t1=t0+dt-self.via.tlag
        if self.doseType == PKPDDose.TYPE_BOLUS:
            if t0<=self.t0 and self.t0<t1:
                return self.doseAmount
            else:
                return 0.0
        elif self.doseType == PKPDDose.TYPE_REPEATED_BOLUS:
            doseAmount=0
            for t in np.arange(self.t0,self.tF,self.every):
                if t0<=t and t<t1:
                    doseAmount+=self.doseAmount
            return doseAmount
        elif self.doseType == PKPDDose.TYPE_INFUSION:
            if t0>self.tF or t1<self.t0:
                return 0.0
            else:
                tLeft=max(t0,self.t0)
                tRight=min(t1,self.tF)
                return self.doseAmount*(tRight-tLeft)

    def getAmountReleasedAt(self,t0,dt=0.5):
        doseAmount = 0.0
        if self.via.viaProfile == None:
            doseAmount += self.getDoseAt(t0,dt)
        else:
            if self.doseType!=PKPDDose.TYPE_INFUSION:
                self.via.viaProfile.Amax = self.via.bioavailability*self.doseAmount
                doseAmount += self.via.viaProfile.getAg(t0-self.t0-self.via.tlag)-\
                              self.via.viaProfile.getAg(t0-self.t0-self.via.tlag+dt)
            else:
                raise Exception("getAmountReleasedAt not implemented for non-iv infusion")
        if doseAmount<0:
            doseAmount=0
        return doseAmount

    def isDoseABolus(self):
        if self.doseType != PKPDDose.TYPE_BOLUS:
            return False
        if self.t0 != 0:
            return False
        return True

    def getTUnitsString(self):
        return self.tunits._toString()

    def getDUnitsString(self):
        return self.dunits._toString()


def createDeltaDose(doseAmount,via,t=0,dunits="mg"):
    dose = PKPDDose()
    dose.doseName = "Bolus"
    dose.doseType = PKPDDose.TYPE_BOLUS
    dose.doseAmount = doseAmount
    dose.via = via
    dose.t0 = t
    dose.tunits = PKPDUnit("min")
    dose.dunits = PKPDUnit(dunits)
    return dose


class DrugSource:
    def __init__(self):
        self.parsedDoseList = []
        self.parameterNames = []
        self.parameterUnits = []
        self.vias = []

    def setDoses(self, parsedDoseList, t0, tF):
        self.originalDoseList = parsedDoseList
        self.parsedDoseList = []
        collectedVias = []
        self.parameterNames = []
        self.vias = []
        for dose in parsedDoseList:
            if dose.doseType == PKPDDose.TYPE_BOLUS:
                self.parsedDoseList.append(dose)
            elif dose.doseType == PKPDDose.TYPE_INFUSION:
                newDose = copy.copy(dose)
                newDose.dunits = copy.copy(dose.dunits)
                newDose.dunits.unit = changeRateToWeight(dose.dunits.unit)
                self.parsedDoseList.append(newDose)
            else:
                for t in np.arange(dose.t0,dose.tF,dose.every):
                    if t0<=t and t<=tF:
                        newDose = copy.copy(dose)
                        newDose.doseType = PKPDDose.TYPE_BOLUS
                        newDose.t0 = t
                        self.parsedDoseList.append(newDose)

            if not dose.via.viaName in collectedVias:
                collectedVias.append(dose.via.viaName)
                viaParameterNames = dose.via.getParameterNames()
                self.parameterNames += viaParameterNames
                self.vias.append((dose.via,len(viaParameterNames)))

    def getAmountReleasedAt(self,t0,dt=0.5):
        doseAmount = 0.0
        for dose in self.parsedDoseList:
            doseAmount+=dose.getAmountReleasedAt(t0,dt)
        return doseAmount

    def getEquation(self):
        retval = ""
        for via,_ in self.vias:
            retval+=via.getEquation()+" "
        return retval.strip()

    def getModelEquation(self):
        retval = ""
        for via,_ in self.vias:
            retval+=via.getModelEquation()+" "
        return retval.strip()

    def getDescription(self):
        retval = ""
        for via,_ in self.vias:
            retval+=via.getDescription()+" "
        return retval.strip()

    def getParameterNames(self):
        return self.parameterNames

    def getNumberOfParameters(self):
        return len(self.getParameterNames())

    def calculateParameterUnits(self,sample):
        if len(self.parameterUnits)>0:
            return self.parameterUnits
        else:
            retval = []
            for via,_ in self.vias:
                retval+=via.calculateParameterUnits(sample)
            self.parameterUnits=retval
            return retval

    def areParametersSignificant(self, lowerBound, upperBound):
        retval = []
        currentToken=0
        for via, viaPrmNo in self.vias:
            retval+=via.areParametersSignificant(lowerBound[currentToken:currentToken+viaPrmNo],
                                                 upperBound[currentToken:currentToken+viaPrmNo])
            currentToken+=viaPrmNo
        if retval:
            return retval
        else:
            return True

    def areParametersValid(self, p):
        retval = True
        currentToken=0
        for via,viaPrmNo in self.vias:
            retval=retval and via.areParametersValid(p[currentToken:currentToken+viaPrmNo])
            currentToken+=viaPrmNo
        return retval

    def setParameters(self, p):
        currentToken=0
        for via,viaPrmNo in self.vias:
            via.setParameters(p[currentToken:currentToken+viaPrmNo])
            currentToken+=viaPrmNo
