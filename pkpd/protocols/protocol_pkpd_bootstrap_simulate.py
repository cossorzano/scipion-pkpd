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

import numpy as np
import math
import random

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDDose, PKPDSample, PKPDVariable
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from .protocol_pkpd_ode_base import ProtPKPDODEBase
from pkpd.pkpd_units import createUnit, multiplyUnits, strUnit
from pkpd.utils import find_nearest

# Tested in test_worokflow_deconvolution.py

class ProtPKPDODESimulate(ProtPKPDODEBase):
    """ Simulate a population of ODE parameters.
    These parameters can be specifically given or from a bootstrap population"""

    _label = 'PK simulate'

    PRM_POPULATION = 0
    PRM_USER_DEFINED = 1

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputODE', params.PointerParam, label="Input ODE model",
                      pointerClass='ProtPKPDMonoCompartment, ProtPKPDTwoCompartments, ProtPKPDMonoCompartmentPD, ProtPKPDTwoCompartmentsBothPD, '\
                                   'ProtPKPDODERefine', help='Select a run of an ODE model')
        form.addParam('paramsSource', params.EnumParam, label="Source of parameters", choices=['ODE Bootstrap','User defined'], default=0,
                      help="Choose a population of parameters or your own")
        form.addParam('inputPopulation', params.PointerParam, label="Input population", condition="paramsSource==0",
                      pointerClass='PKPDFitting', pointerCondition="isPopulation", help='It must be a fitting coming from a bootstrap sample')
        form.addParam('prmUser', params.TextParam, label="Simulation parameters", height=8, default="", condition="paramsSource==1",
                      help='Specify the parameters for the simulation. The parameters must be written in the same order as they are written by the protocol '
                           'that generated the ODE model. Example: \n'
                           'prm1, prm2, prm3, prm4\n'
                           'prmA, prmB, prmC, prmD')
        form.addParam('doses', params.TextParam, label="Doses", height=5, width=50,
                      default="RepeatedBolus ; via=Oral; repeated_bolus; t=0:24:120 h; d=60 mg",
                      help="Structure: [Dose Name] ; [Via] ; [Dose type] ; [time] ; [dose] \n"\
                           "The dose name should have no space or special character\n"\
                           "The via should be one present in the input experiment to the ODE model.\n"\
                           "Valid units are: h, mg, ug, ...\n"\
                           "The description is either a bolus or an infusion as shown in the examples\n"\
                           "\nIt is important that there are two semicolons.\n"\
                           "Examples:\n"\
                           "Infusion0; via=Intravenous; infusion; t=0.500000...0.750000 h; d=60 mg\n"\
                           "Bolus0; via=Oral; bolus; t=0.000000 min; d=60 mg\n"\
                           "RepeatedBolus; via=Oral; repeated_bolus; t=0:24:120 h; d=60 mg")
        form.addParam('t0', params.FloatParam, label="Initial time (h)", default=0)
        form.addParam('tF', params.FloatParam, label="Final time (h)", default=24*7)
        form.addParam('Nsimulations', params.IntParam, label="Simulation samples", default=200, condition="paramsSource==0", expertLevel=LEVEL_ADVANCED,
                      help='Number of simulations')
        form.addParam('addStats', params.BooleanParam, label="Add simulation statistics", default=True, condition="paramsSource==0", expertLevel=LEVEL_ADVANCED,
                      help="Mean, lower and upper confidence levels are added to the output")
        form.addParam('confidenceLevel', params.FloatParam, label="Confidence interval", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters', condition="addStats and paramsSource==0")
        form.addParam('addIndividuals', params.BooleanParam, label="Add individual simulations", default=False, condition="paramsSource==0", expertLevel=LEVEL_ADVANCED,
                      help="Individual simulations are added to the output")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulate',self.inputODE.get().getObjId(), self.Nsimulations.get(),
                                 self.confidenceLevel.get(), self.doses.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addSample(self, sampleName, doseName, simulationsX, y):
        newSample = PKPDSample()
        newSample.sampleName = sampleName
        newSample.variableDictPtr = self.outputExperiment.variables
        newSample.doseDictPtr = self.outputExperiment.doses
        newSample.descriptors = {}
        newSample.doseList = [doseName]
        if type(self.varNameY)!=list:
            newSample.addMeasurementPattern([self.varNameY])
            newSample.addMeasurementColumn("t", simulationsX)
            newSample.addMeasurementColumn(self.varNameY,y)
        else:
            for j in range(len(self.varNameY)):
                newSample.addMeasurementPattern([self.varNameY[j]])
            newSample.addMeasurementColumn("t", simulationsX)
            for j in range(len(self.varNameY)):
                newSample.addMeasurementColumn(self.varNameY[j], y[j])
        newSample.descriptors["AUC0t"] = self.AUC0t
        newSample.descriptors["AUMC0t"] = self.AUMC0t
        newSample.descriptors["MRT"] = self.MRT
        newSample.descriptors["Cmax"] = self.Cmax
        newSample.descriptors["Cmin"] = self.Cmin
        newSample.descriptors["Cavg"] = self.Cavg
        newSample.descriptors["Tmax"] = self.Tmax
        newSample.descriptors["Tmin"] = self.Tmin
        self.outputExperiment.samples[sampleName] = newSample

    def NCA(self, t, C):
        AUClist = []
        AUMClist = []
        Cminlist = []
        Cavglist = []
        Cmaxlist = []
        Tmaxlist = []
        Tminlist = []
        Ndoses=len(self.drugSource.parsedDoseList)
        for ndose in range(0,max(Ndoses-1,1)):
            tperiod0 = self.drugSource.parsedDoseList[ndose].t0
            if ndose+1<Ndoses:
                tperiodF = self.drugSource.parsedDoseList[ndose+1].t0-self.model.deltaT
            else:
                tperiodF =  np.max(t)-1
            idx0 = find_nearest(t,tperiod0)
            idxF = find_nearest(t,tperiodF)

            AUC0t = 0
            AUMC0t = 0
            t0 = t[idx0+1]
            for idx in range(idx0,idxF+1):
                dt = (t[idx+1]-t[idx])
                if C[idx+1]>=C[idx]: # Trapezoidal in the raise
                    AUC0t  += 0.5*dt*(C[idx]+C[idx+1])
                    AUMC0t += 0.5*dt*(C[idx]*t[idx]+C[idx+1]*t[idx+1])
                else: # Log-trapezoidal in the decay
                    decrement = C[idx]/C[idx+1]
                    K = math.log(decrement)
                    B = K/dt
                    AUC0t  += dt*(C[idx]-C[idx+1])/K
                    AUMC0t += (C[idx]*(t[idx]-tperiod0)-C[idx+1]*(t[idx+1]-tperiod0))/B-(C[idx+1]-C[idx])/(B*B)

                if idx==idx0:
                    Cmax=C[idx]
                    Tmax=t[idx]-t0
                    Cmin=C[idx]
                    Tmin=t[idx]-t0
                else:
                    if C[idx]<Cmin:
                        Cmin=C[idx]
                        Tmin=t[idx]-t0
                    elif C[idx]>Cmax:
                        Cmax=C[idx]
                        Tmax=t[idx]-t0
                        if ndose==0:
                            Cmin=C[idx]
                            Tmin=t[idx]-t0
            AUClist.append(AUC0t)
            AUMClist.append(AUMC0t)
            Cminlist.append(Cmin)
            Cmaxlist.append(Cmax)
            Tmaxlist.append(Tmax)
            Tminlist.append(Tmin)
            Cavglist.append(AUC0t/(t[idxF]-t[idx0]))

        print("Fluctuation = Cmax/Cmin")
        print("Accumulation(1) = Cavg(n)/Cavg(1) %")
        print("Accumulation(n) = Cavg(n)/Cavg(n-1) %")
        print("Steady state fraction(n) = Cavg(n)/Cavg(last) %")
        for ndose in range(0,len(AUClist)):
            fluctuation = Cmaxlist[ndose]/Cminlist[ndose]
            if ndose>0:
                accumn = Cavglist[ndose]/Cavglist[ndose-1]
            else:
                accumn = 0
            print("Dose #%d: Cavg= %f [%s] Cmin= %f [%s] Tmin= %d [min] Cmax= %f [%s] Tmax= %d [min] Fluct= %f %% Accum(1)= %f %% Accum(n)= %f %% SSFrac(n)= %f %% AUC= %f [%s] AUMC= %f [%s]"%\
                  (ndose+1,Cavglist[ndose], strUnit(self.Cunits.unit), Cminlist[ndose],strUnit(self.Cunits.unit), int(Tminlist[ndose]), Cmaxlist[ndose], strUnit(self.Cunits.unit),
                   int(Tmaxlist[ndose]), fluctuation*100, Cavglist[ndose]/Cavglist[0]*100, accumn*100, Cavglist[ndose]/Cavglist[-1]*100, AUClist[ndose],strUnit(self.AUCunits),
                   AUMClist[ndose],strUnit(self.AUMCunits)))

        self.AUC0t = float(AUClist[-1])
        self.AUMC0t = float(AUMClist[-1])
        self.MRT = self.AUMC0t/self.AUC0t
        self.Cmin = float(Cminlist[-1])
        self.Cmax = float(Cmaxlist[-1])
        self.Tmin = float(Tminlist[-1])
        self.Tmax = float(Tmaxlist[-1])
        self.Cavg = float(Cavglist[-1])
        self.fluctuation = self.Cmax/self.Cmin if self.Cmin>0 else np.nan
        self.percentageAccumulation = Cavglist[-1]/Cavglist[0] if Cavglist[0]>0 else np.nan

        print("   AUC0t=%f [%s]"%(self.AUC0t,strUnit(self.AUCunits)))
        print("   AUMC0t=%f [%s]"%(self.AUMC0t,strUnit(self.AUMCunits)))
        print("   MRT=%f [min]"%self.MRT)
        print("   Cmax=%f [%s]"%(self.Cmax,strUnit(self.Cunits.unit)))
        print("   Cmin=%f [%s]"%(self.Cmin,strUnit(self.Cunits.unit)))
        print("   Cavg=%f [%s]"%(self.Cavg,strUnit(self.Cunits.unit)))
        print("   Tmax=%f [min]"%self.Tmax)
        print("   Tmin=%f [min]"%self.Tmin)

    def runSimulate(self, objId, Nsimulations, confidenceInterval, doses):
        self.protODE = self.inputODE.get()
        if hasattr(self.protODE, "outputExperiment"):
            self.experiment = self.readExperiment(self.protODE.outputExperiment.fnPKPD)
        elif hasattr(self.protODE, "outputExperiment1"):
            self.experiment = self.readExperiment(self.protODE.outputExperiment1.fnPKPD)
        else:
            raise Exception("Cannot find an outputExperiment in the input ODE")
        if self.paramsSource==ProtPKPDODESimulate.PRM_POPULATION:
            self.fitting = self.readFitting(self.inputPopulation.get().fnFitting, cls="PKPDSampleFitBootstrap")
        else:
            if hasattr(self.protODE, "outputFitting"):
                self.fitting = self.readFitting(self.protODE.outputFitting.fnFitting)
            elif hasattr(self.protODE, "outputFitting1"):
                self.fitting = self.readFitting(self.protODE.outputFitting1.fnFitting)
        self.varNameX = self.fitting.predictor.varName
        if type(self.fitting.predicted)!=list:
            self.varNameY = self.fitting.predicted.varName
        else:
            self.varNameY = [var.varName for var in self.fitting.predicted]

        # Create drug source
        self.clearGroupParameters()
        self.createDrugSource()

        # Create output object
        self.outputExperiment = PKPDExperiment()
        tvar = PKPDVariable()
        tvar.varName = "t"
        tvar.varType = PKPDVariable.TYPE_NUMERIC
        tvar.role = PKPDVariable.ROLE_TIME
        tvar.units = createUnit("min")

        self.outputExperiment.variables[self.varNameX] = tvar
        if type(self.fitting.predicted)!=list:
            self.outputExperiment.variables[self.varNameY] = self.experiment.variables[self.varNameY]
        else:
            for varName in self.varNameY:
                self.outputExperiment.variables[varName] = self.experiment.variables[varName]
        self.outputExperiment.general["title"]="Simulated ODE response"
        self.outputExperiment.general["comment"]="Simulated ODE response"
        self.outputExperiment.vias = self.experiment.vias

        # Setup model
        self.model = self.protODE.createModel()
        self.model.setExperiment(self.outputExperiment)
        if hasattr(self.protODE,"deltaT"):
            self.model.deltaT = self.protODE.deltaT.get()
        self.model.setXVar(self.varNameX)
        self.model.setYVar(self.varNameY)
        Nsamples = int(60*math.ceil((self.tF.get()-self.t0.get())/self.model.deltaT))+1
        self.model.x = [self.t0.get()+i*self.model.deltaT for i in range(0,Nsamples)]
        self.modelList.append(self.model)

        # Read the doses
        for line in self.doses.get().replace('\n',';;').split(';;'):
            tokens = line.split(';')
            if len(tokens)<5:
                print("Skipping dose: ",line)
                continue
            dosename = tokens[0].strip()
            self.outputExperiment.doses[dosename] = PKPDDose()
            self.outputExperiment.doses[dosename].parseTokens(tokens,self.experiment.vias)
        auxSample = PKPDSample()
        auxSample.descriptors = {}
        auxSample.doseDictPtr = self.outputExperiment.doses
        auxSample.variableDictPtr = self.outputExperiment.variables
        auxSample.doseList = [dosename]
        auxSample.interpretDose()
        self.drugSource.setDoses(auxSample.parsedDoseList, self.t0.get()-10, self.tF.get()+10)
        self.model.drugSource = self.drugSource

        # Check units
        # Dunits = self.outputExperiment.doses[dosename].dunits
        # Cunits = self.experiment.variables[self.varNameY].units

        # Process user parameters
        if self.paramsSource==ProtPKPDODESimulate.PRM_POPULATION:
            Nsimulations = self.Nsimulations.get()
        else:
            lines = self.prmUser.get().strip().replace('\n',';;').split(';;')
            Nsimulations = len(lines)
            prmUser = []
            for line in lines:
                tokens = line.strip().split(',')
                prmUser.append([float(token) for token in tokens])

        # Simulate the different responses
        simulationsX = self.model.x
        simulationsY = np.zeros((Nsimulations,len(simulationsX),self.getResponseDimension()))
        AUCarray = np.zeros(Nsimulations)
        AUMCarray = np.zeros(Nsimulations)
        MRTarray = np.zeros(Nsimulations)
        CminArray = np.zeros(Nsimulations)
        CmaxArray = np.zeros(Nsimulations)
        CavgArray = np.zeros(Nsimulations)
        TminArray = np.zeros(Nsimulations)
        TmaxArray = np.zeros(Nsimulations)
        fluctuationArray = np.zeros(Nsimulations)
        percentageAccumulationArray = np.zeros(Nsimulations)
        for i in range(0,Nsimulations):
            self.setTimeRange(None)

            if self.paramsSource==ProtPKPDODESimulate.PRM_POPULATION:
                # Take parameters randomly from the population
                nfit = int(random.uniform(0,len(self.fitting.sampleFits)))
                sampleFit = self.fitting.sampleFits[nfit]
                nprm = int(random.uniform(0,sampleFit.parameters.shape[0]))
                parameters = sampleFit.parameters[nprm,:]
            else:
                parameters = np.asarray(prmUser[i],np.double)
            print("Simulated sample %d: %s"%(i,str(parameters)))

            # Prepare source and this object
            self.drugSource.setDoses(auxSample.parsedDoseList, self.model.t0, self.model.tF)
            self.protODE.configureSource(self.drugSource)
            self.model.drugSource = self.drugSource
            parameterNames = self.getParameterNames() # Necessary to count the number of source and PK parameters

            # Prepare the model
            self.setParameters(parameters)
            y = self.forwardModel(parameters, [simulationsX]*self.getResponseDimension())

            # Create AUC, AUMC, MRT variables and units
            if i == 0:
                if type(self.varNameY)!=list:
                    self.Cunits = self.experiment.variables[self.varNameY].units
                else:
                    self.Cunits = self.experiment.variables[self.varNameY[0]].units
                self.AUCunits = multiplyUnits(tvar.units.unit, self.Cunits.unit)
                self.AUMCunits = multiplyUnits(tvar.units.unit, self.AUCunits)

                if self.addStats or self.addIndividuals:
                    AUCvar = PKPDVariable()
                    AUCvar.varName = "AUC0t"
                    AUCvar.varType = PKPDVariable.TYPE_NUMERIC
                    AUCvar.role = PKPDVariable.ROLE_LABEL
                    AUCvar.units = createUnit(strUnit(self.AUCunits))

                    AUMCvar = PKPDVariable()
                    AUMCvar.varName = "AUMC0t"
                    AUMCvar.varType = PKPDVariable.TYPE_NUMERIC
                    AUMCvar.role = PKPDVariable.ROLE_LABEL
                    AUMCvar.units = createUnit(strUnit(self.AUMCunits))

                    MRTvar = PKPDVariable()
                    MRTvar.varName = "MRT"
                    MRTvar.varType = PKPDVariable.TYPE_NUMERIC
                    MRTvar.role = PKPDVariable.ROLE_LABEL
                    MRTvar.units = createUnit("min")

                    Cmaxvar = PKPDVariable()
                    Cmaxvar.varName = "Cmax"
                    Cmaxvar.varType = PKPDVariable.TYPE_NUMERIC
                    Cmaxvar.role = PKPDVariable.ROLE_LABEL
                    Cmaxvar.units = createUnit(strUnit(self.Cunits.unit))

                    Tmaxvar = PKPDVariable()
                    Tmaxvar.varName = "Tmax"
                    Tmaxvar.varType = PKPDVariable.TYPE_NUMERIC
                    Tmaxvar.role = PKPDVariable.ROLE_LABEL
                    Tmaxvar.units = createUnit("min")

                    Cminvar = PKPDVariable()
                    Cminvar.varName = "Cmin"
                    Cminvar.varType = PKPDVariable.TYPE_NUMERIC
                    Cminvar.role = PKPDVariable.ROLE_LABEL
                    Cminvar.units = createUnit(strUnit(self.Cunits.unit))

                    Tminvar = PKPDVariable()
                    Tminvar.varName = "Tmin"
                    Tminvar.varType = PKPDVariable.TYPE_NUMERIC
                    Tminvar.role = PKPDVariable.ROLE_LABEL
                    Tminvar.units = createUnit("min")

                    Cavgvar = PKPDVariable()
                    Cavgvar.varName = "Cavg"
                    Cavgvar.varType = PKPDVariable.TYPE_NUMERIC
                    Cavgvar.role = PKPDVariable.ROLE_LABEL
                    Cavgvar.units = createUnit(strUnit(self.Cunits.unit))

                    self.outputExperiment.variables["AUC0t"] = AUCvar
                    self.outputExperiment.variables["AUMC0t"] = AUMCvar
                    self.outputExperiment.variables["MRT"] = MRTvar
                    self.outputExperiment.variables["Cmax"] = Cmaxvar
                    self.outputExperiment.variables["Tmax"] = Tmaxvar
                    self.outputExperiment.variables["Cmin"] = Cminvar
                    self.outputExperiment.variables["Tmin"] = Tminvar
                    self.outputExperiment.variables["Cavg"] = Cavgvar

            # Evaluate AUC, AUMC and MRT in the last full period
            self.NCA(self.model.x,y[0])

            # Keep results
            for j in range(self.getResponseDimension()):
                simulationsY[i,:,j] = y[j]
            AUCarray[i] = self.AUC0t
            AUMCarray[i] = self.AUMC0t
            MRTarray[i] = self.MRT
            CminArray[i] = self.Cmin
            CmaxArray[i] = self.Cmax
            CavgArray[i] = self.Cavg
            TminArray[i] = self.Tmin
            TmaxArray[i] = self.Tmax
            fluctuationArray[i] = self.fluctuation
            percentageAccumulationArray[i] = self.percentageAccumulation
            if self.addIndividuals or self.paramsSource==ProtPKPDODESimulate.PRM_USER_DEFINED:
                self.addSample("Simulation_%d"%i, dosename, simulationsX, y[0])

        # Report NCA statistics
        fhSummary = open(self._getPath("summary.txt"),"w")
        alpha_2 = (100-self.confidenceLevel.get())/2
        limits = np.percentile(AUCarray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"AUC %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.AUCunits),np.mean(AUCarray)))
        limits = np.percentile(AUMCarray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"AUMC %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.AUMCunits),np.mean(AUMCarray)))
        limits = np.percentile(MRTarray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"MRT %f%% confidence interval=[%f,%f] [min] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],np.mean(MRTarray)))
        limits = np.percentile(CminArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Cmin %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.Cunits.unit),np.mean(CminArray)))
        limits = np.percentile(CmaxArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Cmax %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.Cunits.unit),np.mean(CmaxArray)))
        limits = np.percentile(CavgArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Cavg %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.Cunits.unit),np.mean(CavgArray)))
        limits = np.percentile(TminArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Tmin %f%% confidence interval=[%f,%f] [min] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],np.mean(TminArray)))
        limits = np.percentile(TmaxArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Tmax %f%% confidence interval=[%f,%f] [min] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],np.mean(TmaxArray)))
        aux = fluctuationArray[~np.isnan(fluctuationArray)]
        limits = np.percentile(aux,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Fluctuation %f%% confidence interval=[%f,%f] [%%] mean=%f"%(self.confidenceLevel.get(),limits[0]*100,limits[1]*100,np.mean(aux)*100))
        aux = percentageAccumulationArray[~np.isnan(percentageAccumulationArray)]
        limits = np.percentile(aux,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Accum(1) %f%% confidence interval=[%f,%f] [%%] mean=%f"%(self.confidenceLevel.get(),limits[0]*100,limits[1]*100,np.mean(aux)*100))
        fhSummary.close()

        # Calculate statistics
        if self.addStats:
            if self.paramsSource!=ProtPKPDODESimulate.PRM_USER_DEFINED:
                limits = np.percentile(simulationsY,[alpha_2,100-alpha_2],axis=0)

                print("Lower limit NCA")
                self.NCA(simulationsX,limits[0])
                self.addSample("LowerLimit", dosename, simulationsX, limits[0])

            print("Mean profile NCA")
            if self.getResponseDimension()==1:
                mu = np.mean(simulationsY,axis=0)
                self.NCA(simulationsX, mu)
            else:
                mu = []
                for j in range(self.getResponseDimension()):
                    mu.append(np.mean(simulationsY[:,:,j],axis=0))
                self.NCA(simulationsX,mu[0])
            self.addSample("Mean", dosename, simulationsX, mu)

            if self.paramsSource!=ProtPKPDODESimulate.PRM_USER_DEFINED:
                print("Upper limit NCA")
                self.NCA(simulationsX,limits[1])
                self.addSample("UpperLimit", dosename, simulationsX, limits[1])

        self.outputExperiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        self._defineSourceRelation(self.inputODE.get(), self.outputExperiment)
        if self.paramsSource==ProtPKPDODESimulate.PRM_POPULATION:
            self._defineSourceRelation(self.inputPopulation.get(), self.outputExperiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        msg.append("Dose: %s"%self.doses.get())
        if self.paramsSource ==  ProtPKPDODESimulate.PRM_POPULATION:
            msg.append("Number of simulations: %d"%self.Nsimulations.get())
        else:
            msg.append("Parameters:\n"+self.prmUser.get())
        msg.append(" ")
        self.addFileContentToMessage(msg,self._getPath("summary.txt"))
        return msg

    def _validate(self):
        msg=[]
        if self.paramsSource ==  ProtPKPDODESimulate.PRM_POPULATION and \
            not self.inputPopulation.get().fnFitting.get().endswith("bootstrapPopulation.pkpd"):
            msg.append("Population must be a bootstrap sample")
        return msg
