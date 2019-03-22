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
Signal Analysis models
"""
import numpy as np
import math
from scipy.optimize import fsolve
from pkpd.objects import PKPDModelBase
from pkpd.pkpd_units import multiplyUnits, divideUnits

class SAModel(PKPDModelBase):
    def calculateParameters(self, show=True):
        pass

class NCAObsIVModel(SAModel):
    def getDescription(self):
        return "Non-compartmental Analysis of intravascular bolus based on observations (%s)"%self.__class__.__name__

    def getParameterDescriptions(self):
        return ['Automatically estimated Area Under the Curve from 0 to t',
                'Automatically estimated Area Under the Curve from 0 to infinity',
                'Automatically estimated Area Under the 1st Moment Curve from 0 to t',
                'Automatically estimated Area Under the 1st Moment Curve from 0 to infinity',
                'Automatically estimated Mean Residence Time',
                'Automatically estimated Apparent Volume of Distribution using measurements from 0 to t',
                'Automatically estimated Apparent Volume of Distribution using measurements from 0 to infinity',
                'Automatically estimated Apparent Volume of Distribution at steady state',
                'Automatically estimated Clearance using measurements from 0 to t',
                'Automatically estimated Clearance using measurements from 0 to infinity',
                'Automatically estimated Half time']

    def getParameterNames(self):
        return ['AUC_0t','AUC_0inf','AUMC_0t','AUMC_0inf','MRT','Vd_0t','Vd_0inf','Vss','CL_0t','CL_0inf', 'thalf']

    def calculateParameterUnits(self, sample):
        tunits = self.experiment.getVarUnits(self.xName)
        Cunits = self.experiment.getVarUnits(self.yName)
        Dunits = sample.getDoseUnits()
        AUCunits = multiplyUnits(tunits,Cunits)
        AUMCunits = multiplyUnits(tunits,AUCunits)
        Vunits = divideUnits(Dunits,Cunits)
        CLunits = divideUnits(Vunits,tunits)
        self.parameterUnits = [AUCunits, AUCunits, AUMCunits, AUMCunits, tunits, Vunits, Vunits, Vunits, CLunits, CLunits, tunits]
        return self.parameterUnits

    def calculateParameters(self, show=True):
        t = self.x[0] # From [array(...)] to array(...)
        C = self.y[0]

        # AUC0t, AUMC0t
        AUC0t = 0
        AUMC0t = 0
        if self.areaCalc == "Trapezoidal":
            for i in range(len(C)-1):
                dt = (t[i+1]-t[i])
                AUC0t  += dt*(C[i]+C[i+1])
                AUMC0t += dt*(C[i]*t[i]+C[i+1]*t[i+1])
            AUC0t*=0.5
            AUMC0t*=0.5
        elif self.areaCalc == "Log-Trapezoidal":
            for i in range(len(C)-1):
                dt = (t[i+1]-t[i])
                decrement = C[i]/C[i+1]
                K = math.log(decrement)
                B = K/dt
                AUC0t  += dt*(C[i]-C[i+1])/K
                # AUMC0t += dt*(C[i+1]*t[i+1]-C[i]*t[i])/(-K)-1.0/decrement*dt*dt/(K*K)
                   # Eq. 7.44 of Domenech Berrozpe, ... Tratado general de biofarmacia y farmacocinetica Vol. 1 (2013)
                # AUMC0t += (C[i]*t[i]-C[i+1]*t[i+1])/B-(C[i]-C[i+1])/(B*B)
                   # http://www.agah.eu/fileadmin/_migrated/content_uploads/PK-glossary_PK_working_group_2004.pdf
                # AUMC0t += 1/K * dt*(C[i]*t[i]+C[i+1]*t[i+1])
                   # Eq. 8.32 Atkinson, Huang, ... Principles of Clinical Pharmacology (2012)
                AUMC0t += (C[i]*t[i]-C[i+1]*t[i+1])/B-(C[i+1]-C[i])/(B*B)
                   # Eq. 2.315 Gabrielsson and Weiner. Pharmacokinetic and Pharmacodynamic data analysis

        # AUC0inf, AUMC0inf
        AUC0inf = AUC0t+C[-1]/self.lambdaz
        AUMC0inf = AUMC0t+C[-1]*(t[-1]+1/self.lambdaz)/self.lambdaz

        # MRT
        MRT = AUMC0inf/AUC0inf

        # Volumes
        Vd0t = self.F*self.D/(AUC0t*self.lambdaz)
        Vd0inf = self.F*self.D/(AUC0inf*self.lambdaz)
        Vss = self.D*AUMC0inf/(AUC0inf*AUC0inf)

        # Clearances
        CL0t = self.F*self.D/AUC0t
        CL0inf = self.F*self.D/AUC0inf

        # Thalf
        thalf = math.log(2.0)*Vd0inf/CL0inf

        # Finish
        self.parameters = []
        self.parameters.append(AUC0t)
        self.parameters.append(AUC0inf)
        self.parameters.append(AUMC0t)
        self.parameters.append(AUMC0inf)
        self.parameters.append(MRT)
        self.parameters.append(Vd0t)
        self.parameters.append(Vd0inf)
        self.parameters.append(Vss)
        self.parameters.append(CL0t)
        self.parameters.append(CL0inf)
        self.parameters.append(thalf)

class NCAExpIVModel(SAModel):
    def getDescription(self):
        return "Non-compartmental Analysis of intravascular bolus based on exponential fitting (%s)"%self.__class__.__name__

    def getParameterDescriptions(self):
        return ['Automatically estimated Slope of decay at the end of the curve',
                'Automatically estimated Last decay rate (lambda_n)',
                'Automatically estimated Initial concentration at time 0',
                'Automatically estimated Area Under the Curve from 0 to infinity',
                'Automatically estimated Area Under the 1st Moment Curve from 0 to infinity',
                'Automatically estimated Mean Residence Time',
                'Automatically estimated Apparent Volume of the Central Compartment',
                'Automatically estimated Apparent Volume of Distribution',
                'Automatically estimated Apparent Volume of Distribution at steady state',
                'Automatically estimated Clearance',
                'Automatically estimated Half time'
                ]

    def getParameterNames(self):
        return ['slope','rate_constant','C0','AUC_0inf','AUMC_0inf','MRT','Vc','Vd','Vss','CL','thalf']

    def calculateParameterUnits(self, sample):
        tunits = self.experiment.getVarUnits(self.xName)
        Cunits = self.experiment.getVarUnits(self.yName)
        Dunits = sample.getDoseUnits()
        AUCunits = multiplyUnits(tunits,Cunits)
        AUMCunits = multiplyUnits(tunits,AUCunits)
        Vunits = divideUnits(Dunits,Cunits)
        CLunits = divideUnits(Vunits,tunits)
        lambdaUnits = self.lambdanUnits.unit
        self.parameterUnits = [lambdaUnits, lambdaUnits, Cunits, AUCunits, AUMCunits, tunits, Vunits, Vunits, Vunits, CLunits, tunits]
        return self.parameterUnits

    def calculateParameters(self, show=True):
        Cn = self.Cn
        lambdan = self.lambdan

        # Rate constant
        rateConstant = lambdan[-1]

        # Slope
        slope = -rateConstant/math.log(10.0)

        # C0
        C0 = sum(Cn)

        # AUC0inf, AUMC0inf
        AUC0inf = sum([ci/lambdai for (ci, lambdai) in zip(Cn, lambdan)])
        AUMC0inf = sum([ci/(lambdai*lambdai) for (ci, lambdai) in zip(Cn, lambdan)])

        # MRT
        # COSS: I think it is incorrect MRT = sum([1.0/lambdai for lambdai in lambdan])
        MRT = AUMC0inf/AUC0inf

        # Volumes
        Vc = self.D/C0
        Vd = self.F*self.D/(AUC0inf*self.lambdaz)
        Vss = self.D*AUMC0inf/(AUC0inf*AUC0inf)

        # Clearances
        CL0inf = self.F*self.D/AUC0inf

        # thalf
        # COSS: I think it is incorrect thalf = math.log(2.0)/rateConstant
        thalf = math.log(2.0)*Vd/CL0inf
        npLambdan = -np.asarray(lambdan,np.double)
        func = lambda thalf : 0.5*C0 - np.dot(Cn,np.exp(npLambdan*thalf))
        thalf = fsolve(func, thalf)
        thalf=thalf[0]

        # Finish
        self.parameters = []
        self.parameters.append(slope)
        self.parameters.append(rateConstant)
        self.parameters.append(C0)
        self.parameters.append(AUC0inf)
        self.parameters.append(AUMC0inf)
        self.parameters.append(MRT)
        self.parameters.append(Vc)
        self.parameters.append(Vd)
        self.parameters.append(Vss)
        self.parameters.append(CL0inf)
        self.parameters.append(thalf)

class NCAEVModel(SAModel):
    def getDescription(self):
        return "Non-compartmental Analysis of extravascular bolus (%s)"%self.__class__.__name__

    def getParameterDescriptions(self):
        return ['Automatically estimated Area Under the Curve from 0 to t',
                'Automatically estimated Area Under the Curve from 0 to infinity',
                'Automatically estimated Area Under the 1st Moment Curve from 0 to t',
                'Automatically estimated Area Under the 1st Moment Curve from 0 to infinity',
                'Automatically estimated Mean Residence Time',
                'Automatically estimated Apparent Volume of Distribution using measurements from 0 to t',
                'Automatically estimated Apparent Volume of Distribution using measurements from 0 to infinity',
                'Automatically estimated Clearance using measurements from 0 to t',
                'Automatically estimated Clearance using measurements from 0 to infinity',
                'Automatically estimated Half time']

    def getParameterNames(self):
        return ['AUC_0t','AUC_0inf','AUMC_0t','AUMC_0inf','MRT','Vd_0t','Vd_0inf','CL_0t','CL_0inf','thalf']

    def calculateParameterUnits(self, sample):
        tunits = self.experiment.getVarUnits(self.xName)
        Cunits = self.experiment.getVarUnits(self.yName)
        Dunits = sample.getDoseUnits()
        AUCunits = multiplyUnits(tunits,Cunits)
        AUMCunits = multiplyUnits(tunits,AUCunits)
        Vunits = divideUnits(Dunits,Cunits)
        CLunits = divideUnits(Vunits,tunits)
        self.parameterUnits = [AUCunits, AUCunits, AUMCunits, AUMCunits, tunits, Vunits, Vunits, CLunits, CLunits, tunits]

    def calculateParameters(self, show=True):
        t = np.concatenate([[0],self.x[0]]) # From [array(...)] to array(...)
        C = np.concatenate([[0],self.y[0]])

        # AUC0t, AUMC0t
        AUC0t = 0
        AUMC0t = 0
        if self.areaCalc == "Trapezoidal":
            for i in range(len(C)-1):
                dt = (t[i+1]-t[i])
                AUC0t  += dt*(C[i]+C[i+1])
                AUMC0t += dt*(C[i]*t[i]+C[i+1]*t[i+1])
            AUC0t*=0.5
            AUMC0t*=0.5
        elif self.areaCalc == "Mixed":
            for i in range(len(C)-1):
                dt = (t[i+1]-t[i])
                if dt==0:
                    continue
                if C[i+1]>=C[i]: # Trapezoidal in the raise
                    AUC0t  += 0.5*dt*(C[i]+C[i+1])
                    AUMC0t += 0.5*dt*(C[i]*t[i]+C[i+1]*t[i+1])
                else: # Log-trapezoidal in the decay
                    decrement = C[i]/C[i+1]
                    K = math.log(decrement)
                    B = K/dt
                    AUC0t  += dt*(C[i]-C[i+1])/K
                    AUMC0t += (C[i]*t[i]-C[i+1]*t[i+1])/B-(C[i+1]-C[i])/(B*B)
                       # Eq. 2.315 Gabrielsson and Weiner. Pharmacokinetic and Pharmacodynamic data analysis

        # AUC0inf, AUMC0inf
        AUC0inf = AUC0t+C[-1]/self.Ke
        AUMC0inf = AUMC0t+C[-1]*(t[-1]+1/self.Ke)/self.Ke

        # MRT
        MRT = AUMC0inf/AUC0inf

        # Volumes
        Vd0t = self.F*self.D/(AUC0t*self.Ke)
        Vd0inf = self.F*self.D/(AUC0inf*self.Ke)

        # Clearances
        CL0t = self.F*self.D/AUC0t
        CL0inf = self.F*self.D/AUC0inf

        # Thalf
        thalf = math.log(2.0)/self.Ke

        # Finish
        self.parameters = []
        self.parameters.append(AUC0t)
        self.parameters.append(AUC0inf)
        self.parameters.append(AUMC0t)
        self.parameters.append(AUMC0inf)
        self.parameters.append(MRT)
        self.parameters.append(Vd0t)
        self.parameters.append(Vd0inf)
        self.parameters.append(CL0t)
        self.parameters.append(CL0inf)
        self.parameters.append(thalf)
