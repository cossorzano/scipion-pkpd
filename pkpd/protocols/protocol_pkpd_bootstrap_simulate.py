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

import math
import numpy as np
import openpyxl
import random
from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDDose, PKPDSample, PKPDVariable
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from .protocol_pkpd_ode_base import ProtPKPDODEBase
from pkpd.pkpd_units import createUnit, multiplyUnits, divideUnits, strUnit, PKPDUnit, unitFromString
from pkpd.utils import find_nearest, excelWriteRow
from pkpd.biopharmaceutics import PKPDVia
from pkpd.models.pk_models import PK_Monocompartment, PK_Twocompartments, PK_TwocompartmentsClintCl

# Tested in test_workflow_deconvolution.py
# Tested in test_workflow_ivivc.py
# Tested in test_workflow_ivivc2.py

class ProtPKPDODESimulate(ProtPKPDODEBase):
    """ Simulate a population of ODE parameters.
    These parameters can be specifically given, from a bootstrap population, a previous fitting, or an experiment.

    AUC0t and AUMC0t are referred to each dose (that is, t=0 is the beginning of the dose). Ctau is the concentration
    at the end of the dose.

    Tmin, Tmax, and Ttau are referred to the beginning of the dose.

    This protocol writes an Excel file (nca.xlsx) in which all the simulations and doses are written."""

    _label = 'PK simulate'

    PRM_POPULATION = 0
    PRM_USER_DEFINED = 1
    PRM_FITTING = 2
    PRM_EXPERIMENT = 3

    SRC_ODE = 0
    SRC_LIST = 1

    VIATYPE_IV = 0
    VIATYPE_EV0 = 1
    VIATYPE_EV1 = 2

    PKTYPE_COMP1 = 0
    PKTYPE_COMP2 = 1
    PKTYPE_COMP2CLCLINT = 2

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('odeSource', params.EnumParam, label='Source of ODE type',
                      choices=['A previous PK fit', 'List'], default=self.SRC_ODE,
                      help='The input model is used to define the absorbtion via and the type of distribution and elimination. '
                           'Its parameters are not important, only its types')
        form.addParam('inputODE', params.PointerParam, label="Input ODE model", condition='odeSource==0',
                      pointerClass='ProtPKPDMonoCompartment, ProtPKPDTwoCompartments, ProtPKPDMonoCompartmentPD, ProtPKPDTwoCompartmentsBothPD, '\
                                   'ProtPKPDODERefine, ProtPKPDTwoCompartmentsClint, ProtPKPDTwoCompartmentsClintCl',
                      help='Select a run of an ODE model. This input is used to learn the kind of model that is needed. The parameters are '
                           'taken from the sources below')
        form.addParam('viaType',params.EnumParam, label='Input via', condition='odeSource==1', default=0,
                      choices=['Intravenous (iv)','Extra vascular, 0th order absorption (ev0)',
                               'Extra vascular, 1st order absorption (ev1)'],
                      help='Parameters:\n'
                           'iv (vianame=Intravenous): no parameters\n'
                           'ev0 (vianame=Oral): Rin, tlag, F (bioavailability)\n'
                           'ev1 (vianame=Oral): Ka, tlag, F (bioavailability)\n')
        form.addParam('viaPrm', params.TextParam, label="Via parameters", height=8, default="",
                      condition="odeSource==1",
                      help='Specify the parameters for the via separated by commas. Example: \n'
                           'prm1, prm2, prm3\n'
                           '\n'
                           'Parameters:\n'
                           'iv: no parameters\n'
                           'ev0: Rin, tlag, F (bioavailability)\n'
                           'ev1: Ka, tlag, F (bioavailability)\n')
        form.addParam('pkType',params.EnumParam, label='PK model', condition='odeSource==1', default=0,
                      choices=['1 compartment', '2 compartments', '2 compartments Cl+Clint'],
                      help='Parameters:\n'
                           '1 compartment: Cl, V\n'
                           '2 compartments: Cl, V, Clp, Vp\n'
                           '2 compartments Cl+Clint: Vmax, Km, Cl, V, Clp, Vp\n')
        form.addParam('timeUnits',params.StringParam,label='Time units', default='min', condition="odeSource==1",
                      help='min or h')
        form.addParam('volumeUnits',params.StringParam,label='Volume units', default='L', condition="odeSource==1",
                      help='mL or L')
        form.addParam('paramsSource', params.EnumParam, label="Source of parameters", condition="odeSource==0",
                      choices=['ODE Bootstrap','User defined','ODE Fitting','Experiment'], default=0,
                      help="Choose a population of parameters, your own parameters or a previously fitted set of measurements")
        form.addParam('inputPopulation', params.PointerParam, label="Input population",
                      condition="paramsSource==0 and odeSource==0",
                      pointerClass='PKPDFitting', pointerCondition="isPopulation",
                      help='It must be a fitting coming from a bootstrap sample')
        form.addParam('prmUser', params.TextParam, label="Simulation parameters", height=8, default="",
                      condition="paramsSource==1 or odeSource==1",
                      help='Specify the parameters for the simulation separated by commas. '
                           'The parameters must be written in the same order as they are written by the protocol '
                           'that generated the ODE model. Example: \n'
                           'prm1, prm2, prm3, prm4\n'
                           'prmA, prmB, prmC, prmD\n'
                           '\n'
                           'If the parameters are not taken from a previous ODE:\n'
                           '1 compartment: Cl, V\n'
                           '2 compartments: Cl, V, Clp, Vp\n'
                           '2 compartments Cl+Clint: Vmax, Km, Cl, V, Clp, Vp\n')
        form.addParam('inputFitting', params.PointerParam, label="Input ODE fitting",
                      condition="paramsSource==2 and odeSource==0",
                      pointerClass='PKPDFitting', help='It must be a fitting coming from a compartmental PK fitting')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      condition="paramsSource==3 and odeSource==0",
                      pointerClass='PKPDExperiment',
                      help='It must contain the parameters required for the simulation (Cl, V, Clp, Vp, ...)')
        form.addParam('doses', params.TextParam, label="Doses", height=5, width=50,
                      default="RepeatedBolus ; via=Oral; repeated_bolus; t=0:24:120 h; d=60 mg",
                      help="Structure: [Dose Name] ; [Via] ; [Dose type] ; [time] ; [dose] \n"\
                           "The dose name should have no space or special character\n"\
                           "The via should be one present in the input experiment to the ODE model.\n"\
                           "Valid units are: h, mg, ug, ...\n"\
                           "The description is either a bolus or an infusion as shown in the examples\n"\
                           "\nIt is important that there are two semicolons.\n"\
                           "Examples:\n"\
                           "Infusion0; via=Intravenous; infusion; t=0.500000:0.750000 h; d=60 mg/h\n"\
                           "Bolus0; via=Oral; bolus; t=0.000000 min; d=60 mg\n"\
                           "RepeatedBolus; via=Oral; repeated_bolus; t=0:24:120 h; d=60 mg")
        form.addParam('t0', params.FloatParam, label="Initial time (see help)", default=0,
                      help="Same units as input experiment")
        form.addParam('tF', params.FloatParam, label="Final time (see help)", default=24*7,
                      help="Same units as input experiment")
        form.addParam('sampling', params.FloatParam, label="Sampling time (see help)", default=1,
                      help="Same units as input experiment")
        form.addParam('Nsimulations', params.IntParam, label="Simulation samples", default=200,
                      condition="paramsSource==0 and odeSource==0",
                      expertLevel=LEVEL_ADVANCED, help='Number of simulations')
        form.addParam('addStats', params.BooleanParam, label="Add simulation statistics", default=True,
                      condition="paramsSource==0 and odeSource==0",
                      expertLevel=LEVEL_ADVANCED,
                      help="Mean, lower and upper confidence levels are added to the output")
        form.addParam('confidenceLevel', params.FloatParam, label="Confidence interval", default=95,
                      expertLevel=LEVEL_ADVANCED, help='Confidence interval for the fitted parameters',
                      condition="addStats and paramsSource==0 and odeSource==0")
        form.addParam('addIndividuals', params.BooleanParam, label="Add individual simulations", default=False,
                      expertLevel=LEVEL_ADVANCED, help="Individual simulations are added to the output")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulate', self.Nsimulations.get(), self.confidenceLevel.get(), self.doses.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def addSample(self, sampleName, doseName, simulationsX, y):
        newSample = PKPDSample()
        newSample.sampleName = sampleName
        newSample.variableDictPtr = self.outputExperiment.variables
        newSample.doseDictPtr = self.outputExperiment.doses
        newSample.descriptors = {}
        newSample.doseList = [doseName]
        tsample = np.arange(0.0, np.max(simulationsX), self.sampling.get())
        if type(self.varNameY)!=list:
            newSample.addMeasurementPattern([self.varNameY])
            B=InterpolatedUnivariateSpline(simulationsX, y, k=1)
            newSample.addMeasurementColumn("t", tsample)
            newSample.addMeasurementColumn(self.varNameY,B(tsample))
        else:
            for j in range(len(self.varNameY)):
                newSample.addMeasurementPattern([self.varNameY[j]])
            newSample.addMeasurementColumn("t", tsample)
            for j in range(len(self.varNameY)):
                B = InterpolatedUnivariateSpline(simulationsX, y[j], k=1)
                newSample.addMeasurementColumn(self.varNameY[j], B(tsample))
        newSample.descriptors["FromSample"] = self.fromSample
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
        Ctaulist = []
        Tmaxlist = []
        Tminlist = []
        Ttaulist = []
        Ndoses=len(self.drugSource.parsedDoseList)
        for ndose in range(0,max(Ndoses,1)):
            tperiod0 = self.drugSource.parsedDoseList[ndose].t0
            if ndose+1<Ndoses:
                tperiodF = self.drugSource.parsedDoseList[ndose+1].t0-self.model.deltaT
            else:
                tperiodF =  np.max(t)-1
            idx0 = find_nearest(t,tperiod0)
            idxF = find_nearest(t,tperiodF)
            if idxF>=len(t)-1:
                idxF=len(t)-2

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
                        # if ndose==0:
                        #     Cmin=C[idx]
                        #     Tmin=t[idx]-t0
            AUClist.append(AUC0t)
            AUMClist.append(AUMC0t)
            Cminlist.append(Cmin)
            Cmaxlist.append(Cmax)
            Ctaulist.append(C[idxF])
            Ttaulist.append(t[idxF]-t0)
            Tmaxlist.append(Tmax)
            Tminlist.append(Tmin)
            Cavglist.append(AUC0t/(t[idxF]-t[idx0]))

        print("Fluctuation = Cmax/Cmin")
        print("Accumulation(1) = Cavg(n)/Cavg(1) %")
        print("Accumulation(n) = Cavg(n)/Cavg(n-1) %")
        print("Steady state fraction(n) = Cavg(n)/Cavg(last) %")
        for ndose in range(0,len(AUClist)):
            if Cminlist[ndose]!=0:
                fluctuation = Cmaxlist[ndose]/Cminlist[ndose]
            else:
                fluctuation = np.nan
            if ndose>0:
                accumn = Cavglist[ndose]/Cavglist[ndose-1]
            else:
                accumn = 0
            print("Dose #%d t=[%f,%f]: Cavg= %f [%s] Cmin= %f [%s] Tmin= %d [%s] Cmax= %f [%s] Tmax= %d [%s] Ctau= %f [%s] Ttau= %d [%s] Fluct= %f %% Accum(1)= %f %% Accum(n)= %f %% SSFrac(n)= %f %% AUC= %f [%s] AUMC= %f [%s]"%\
                  (ndose+1,t[idx0],t[idxF],Cavglist[ndose], strUnit(self.Cunits.unit), Cminlist[ndose],
                   strUnit(self.Cunits.unit), int(Tminlist[ndose]), strUnit(self.outputExperiment.getTimeUnits().unit),
                   Cmaxlist[ndose], strUnit(self.Cunits.unit),
                   int(Tmaxlist[ndose]), strUnit(self.outputExperiment.getTimeUnits().unit),
                   Ctaulist[ndose], strUnit(self.Cunits.unit), int(Ttaulist[ndose]),
                   strUnit(self.outputExperiment.getTimeUnits().unit),
                   fluctuation*100, Cavglist[ndose]/Cavglist[0]*100, accumn*100, Cavglist[ndose]/Cavglist[-1]*100,
                   AUClist[ndose],strUnit(self.AUCunits),
                   AUMClist[ndose],strUnit(self.AUMCunits)))

        self.AUC0t = float(AUClist[-1])
        self.AUMC0t = float(AUMClist[-1])
        self.MRT = self.AUMC0t/self.AUC0t
        self.Cmin = float(Cminlist[-1])
        self.Cmax = float(Cmaxlist[-1])
        self.Tmin = float(Tminlist[-1])
        self.Tmax = float(Tmaxlist[-1])
        self.Cavg = float(Cavglist[-1])
        self.Ctau = float(Ctaulist[-1])
        self.Ttau = float(Ttaulist[-1])
        self.fluctuation = self.Cmax/self.Cmin if self.Cmin>0 else np.nan
        self.percentageAccumulation = Cavglist[-1]/Cavglist[0] if Cavglist[0]>0 else np.nan

        print("   AUC0t=%f [%s]"%(self.AUC0t,strUnit(self.AUCunits)))
        print("   AUMC0t=%f [%s]"%(self.AUMC0t,strUnit(self.AUMCunits)))
        print("   MRT=%f [%s]"%(self.MRT,strUnit(self.outputExperiment.getTimeUnits().unit)))
        print("   Cmax=%f [%s]"%(self.Cmax,strUnit(self.Cunits.unit)))
        print("   Cmin=%f [%s]"%(self.Cmin,strUnit(self.Cunits.unit)))
        print("   Cavg=%f [%s]"%(self.Cavg,strUnit(self.Cunits.unit)))
        print("   Ctau=%f [%s]"%(self.Ctau,strUnit(self.Cunits.unit)))
        print("   Tmax=%f [%s]"%(self.Tmax,strUnit(self.outputExperiment.getTimeUnits().unit)))
        print("   Tmin=%f [%s]"%(self.Tmin,strUnit(self.outputExperiment.getTimeUnits().unit)))
        print("   Ttau=%f [%s]"%(self.Ttau,strUnit(self.outputExperiment.getTimeUnits().unit)))
        return AUClist, AUMClist, Cminlist, Cavglist, Cmaxlist, Ctaulist, Tmaxlist, Tminlist, Ttaulist

    def runSimulate(self, Nsimulations, confidenceInterval, doses):
        if self.odeSource.get()==self.SRC_ODE:
            self.protODE = self.inputODE.get()
            if hasattr(self.protODE, "outputExperiment"):
                self.experiment = self.readExperiment(self.protODE.outputExperiment.fnPKPD)
            elif hasattr(self.protODE, "outputExperiment1"):
                self.experiment = self.readExperiment(self.protODE.outputExperiment1.fnPKPD)
            else:
                raise Exception("Cannot find an outputExperiment in the input ODE")
            if self.paramsSource.get()==ProtPKPDODESimulate.PRM_POPULATION:
                self.fitting = self.readFitting(self.inputPopulation.get().fnFitting, cls="PKPDSampleFitBootstrap")
            elif self.paramsSource.get()==ProtPKPDODESimulate.PRM_FITTING:
                self.fitting = self.readFitting(self.inputFitting.get().fnFitting)
            else:
                # User defined or experiment
                if hasattr(self.protODE, "outputFitting"):
                    self.fitting = self.readFitting(self.protODE.outputFitting.fnFitting)
                elif hasattr(self.protODE, "outputFitting1"):
                    self.fitting = self.readFitting(self.protODE.outputFitting1.fnFitting)
            self.varNameX = self.fitting.predictor.varName
            if type(self.fitting.predicted)!=list:
                self.varNameY = self.fitting.predicted.varName
            else:
                self.varNameY = [var.varName for var in self.fitting.predicted]
            tunits = self.experiment.getTimeUnits().unit
        else:
            self.varNameX = 't'
            self.varNameY = 'C'
            self.fitting = None
            self.experiment = None
            self.protODE = None
            tunits = unitFromString(self.timeUnits.get())

        # Create drug source
        self.clearGroupParameters()
        self.createDrugSource()

        # Create output object
        self.outputExperiment = PKPDExperiment()
        self.outputExperiment.general["title"]="Simulated ODE response"
        self.outputExperiment.general["comment"]="Simulated ODE response"

        # Create the predictor variable
        tvar = PKPDVariable()
        tvar.varName = "t"
        tvar.varType = PKPDVariable.TYPE_NUMERIC
        tvar.role = PKPDVariable.ROLE_TIME
        tvar.units = createUnit(tunits)
        self.outputExperiment.variables[self.varNameX] = tvar

        # Vias
        if self.odeSource.get() == self.SRC_ODE:
            self.outputExperiment.vias = self.experiment.vias
        else:
            # "[ViaName]; [ViaType]; [tlag]; [bioavailability]"
            viaPrmList = [token for token in self.viaPrm.get().strip().split(',')]
            if self.viaType.get()==self.VIATYPE_IV:
                tokens = ["Intravenous", "iv", "tlag=0 min", "bioavailability=1"]
            elif self.viaType.get()==self.VIATYPE_EV0:
                tokens = ["Oral", "ev0"]+["tlag="+viaPrmList[-2].strip()+" "+self.timeUnits.get()]+["bioavailability="+viaPrmList[-1].strip()]
            elif self.viaType.get()==self.VIATYPE_EV1:
                tokens = ["Oral", "ev1"]+["tlag="+viaPrmList[-2].strip()+" "+self.timeUnits.get()]+["bioavailability="+viaPrmList[-1].strip()]

            vianame = tokens[0]
            self.outputExperiment.vias[vianame] = PKPDVia(ptrExperiment=self.outputExperiment)
            self.outputExperiment.vias[vianame].parseTokens(tokens)

        # Read the doses
        dunits = PKPDUnit.UNIT_NONE
        for line in self.doses.get().replace('\n',';;').split(';;'):
            tokens = line.split(';')
            if len(tokens)<5:
                print("Skipping dose: ",line)
                continue
            dosename = tokens[0].strip()
            self.outputExperiment.doses[dosename] = PKPDDose()
            self.outputExperiment.doses[dosename].parseTokens(tokens,self.outputExperiment.vias)
            dunits = self.outputExperiment.doses[dosename].getDoseUnits()

        # Create predicted variables
        if self.odeSource.get()==self.SRC_ODE:
            if type(self.fitting.predicted)!=list:
                self.outputExperiment.variables[self.varNameY] = self.experiment.variables[self.varNameY]
            else:
                for varName in self.varNameY:
                    self.outputExperiment.variables[varName] = self.experiment.variables[varName]
        else:
            Cvar = PKPDVariable()
            Cvar.varName = "C"
            Cvar.varType = PKPDVariable.TYPE_NUMERIC
            Cvar.role = PKPDVariable.ROLE_MEASUREMENT
            Cvar.units = createUnit(divideUnits(dunits.unit,unitFromString(self.volumeUnits.get())))
            self.outputExperiment.variables[self.varNameY] = Cvar

        # Setup model
        if self.odeSource.get() == self.SRC_ODE:
            self.model = self.protODE.createModel()
            if hasattr(self.protODE, "deltaT"):
                self.model.deltaT = self.protODE.deltaT.get()
        else:
            if self.pkType.get() == self.PKTYPE_COMP1:
                self.model = PK_Monocompartment()
            elif self.pkType.get() == self.PKTYPE_COMP2:
                self.model = PK_Twocompartment()
            elif self.pkType.get() == self.PKTYPE_COMP2CLCLINT:
                self.model = PK_TwocompartmentsClintCl()

        self.model.setExperiment(self.outputExperiment)
        self.model.setXVar(self.varNameX)
        self.model.setYVar(self.varNameY)
        Nsamples = int(math.ceil((self.tF.get()-self.t0.get())/self.model.deltaT))+1
        self.model.x = [self.t0.get()+i*self.model.deltaT for i in range(0,Nsamples)]
        self.modelList.append(self.model)

        auxSample = PKPDSample()
        auxSample.descriptors = {}
        auxSample.doseDictPtr = self.outputExperiment.doses
        auxSample.variableDictPtr = self.outputExperiment.variables
        auxSample.doseList = self.outputExperiment.doses.keys()
        auxSample.interpretDose()
        self.drugSource.setDoses(auxSample.parsedDoseList, self.t0.get()-10, self.tF.get()+10)
        self.model.drugSource = self.drugSource

        # Check units
        # Dunits = self.outputExperiment.doses[dosename].dunits
        # Cunits = self.experiment.variables[self.varNameY].units

        # Process user parameters
        if self.odeSource.get()==self.SRC_ODE:
            if self.paramsSource==ProtPKPDODESimulate.PRM_POPULATION:
                Nsimulations = self.Nsimulations.get()
            elif self.paramsSource == ProtPKPDODESimulate.PRM_USER_DEFINED:
                lines = self.prmUser.get().strip().replace('\n',';;').split(';;')
                Nsimulations = len(lines)
                prmUser = []
                for line in lines:
                    tokens = line.strip().split(',')
                    prmUser.append([float(token) for token in tokens])
            elif self.paramsSource == ProtPKPDODESimulate.PRM_FITTING:
                Nsimulations = len(self.fitting.sampleFits)
            else:
                self.inputExperiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
                Nsimulations = len(self.inputExperiment.samples)
                inputSampleNames = [x for x in self.inputExperiment.samples.keys()]
        else:
            lines = self.prmUser.get().strip().replace('\n', ';;').split(';;')
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
        CtauArray = np.zeros(Nsimulations)
        TminArray = np.zeros(Nsimulations)
        TmaxArray = np.zeros(Nsimulations)
        TtauArray = np.zeros(Nsimulations)
        fluctuationArray = np.zeros(Nsimulations)
        percentageAccumulationArray = np.zeros(Nsimulations)

        wb = openpyxl.Workbook()
        wb.active.title = "Simulations"
        for i in range(0,Nsimulations):
            self.setTimeRange(None)

            if self.odeSource.get() == self.SRC_ODE:
                if self.paramsSource==ProtPKPDODESimulate.PRM_POPULATION:
                    # Take parameters randomly from the population
                    nfit = int(random.uniform(0,len(self.fitting.sampleFits)))
                    sampleFit = self.fitting.sampleFits[nfit]
                    nprm = int(random.uniform(0,sampleFit.parameters.shape[0]))
                    parameters = sampleFit.parameters[nprm,:]
                    self.fromSample="Population %d"%nprm
                elif self.paramsSource==ProtPKPDODESimulate.PRM_USER_DEFINED:
                    parameters = np.asarray(prmUser[i],np.double)
                    self.fromSample="User defined"
                elif self.paramsSource == ProtPKPDODESimulate.PRM_FITTING:
                    parameters = self.fitting.sampleFits[i].parameters
                    self.fromSample = self.fitting.sampleFits[i].sampleName
                else:
                    parameters = []
                    self.fromSample=inputSampleNames[i]
                    sample = self.inputExperiment.samples[self.fromSample]
                    for prmName in self.fitting.modelParameters:
                        parameters.append(float(sample.getDescriptorValue(prmName)))
            else:
                parameters = np.asarray(viaPrmList[:-2] + prmUser[i], np.double)
                self.fromSample = "User defined"

            print("From sample name: %s"%self.fromSample)
            print("Simulated sample %d: %s"%(i,str(parameters)))

            # Prepare source and this object
            self.drugSource.setDoses(auxSample.parsedDoseList, self.model.t0, self.model.tF)
            if self.protODE is not None:
                self.protODE.configureSource(self.drugSource)
            self.model.drugSource = self.drugSource
            parameterNames = self.getParameterNames() # Necessary to count the number of source and PK parameters

            # Prepare the model
            self.setParameters(parameters)
            y = self.forwardModel(parameters, [simulationsX]*self.getResponseDimension())

            # Create AUC, AUMC, MRT variables and units
            if i == 0:
                if type(self.varNameY)!=list:
                    self.Cunits = self.outputExperiment.variables[self.varNameY].units
                else:
                    self.Cunits = self.outputExperiment.variables[self.varNameY[0]].units
                self.AUCunits = multiplyUnits(tvar.units.unit, self.Cunits.unit)
                self.AUMCunits = multiplyUnits(tvar.units.unit, self.AUCunits)

                if self.addStats or self.addIndividuals:
                    fromvar = PKPDVariable()
                    fromvar.varName = "FromSample"
                    fromvar.varType = PKPDVariable.TYPE_TEXT
                    fromvar.role = PKPDVariable.ROLE_LABEL
                    fromvar.units = createUnit("none")

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
                    MRTvar.units = createUnit(self.outputExperiment.getTimeUnits().unit)

                    Cmaxvar = PKPDVariable()
                    Cmaxvar.varName = "Cmax"
                    Cmaxvar.varType = PKPDVariable.TYPE_NUMERIC
                    Cmaxvar.role = PKPDVariable.ROLE_LABEL
                    Cmaxvar.units = createUnit(strUnit(self.Cunits.unit))

                    Tmaxvar = PKPDVariable()
                    Tmaxvar.varName = "Tmax"
                    Tmaxvar.varType = PKPDVariable.TYPE_NUMERIC
                    Tmaxvar.role = PKPDVariable.ROLE_LABEL
                    Tmaxvar.units = createUnit(self.outputExperiment.getTimeUnits().unit)

                    Cminvar = PKPDVariable()
                    Cminvar.varName = "Cmin"
                    Cminvar.varType = PKPDVariable.TYPE_NUMERIC
                    Cminvar.role = PKPDVariable.ROLE_LABEL
                    Cminvar.units = createUnit(strUnit(self.Cunits.unit))

                    Tminvar = PKPDVariable()
                    Tminvar.varName = "Tmin"
                    Tminvar.varType = PKPDVariable.TYPE_NUMERIC
                    Tminvar.role = PKPDVariable.ROLE_LABEL
                    Tminvar.units = createUnit(self.outputExperiment.getTimeUnits().unit)

                    Cavgvar = PKPDVariable()
                    Cavgvar.varName = "Cavg"
                    Cavgvar.varType = PKPDVariable.TYPE_NUMERIC
                    Cavgvar.role = PKPDVariable.ROLE_LABEL
                    Cavgvar.units = createUnit(strUnit(self.Cunits.unit))

                    self.outputExperiment.variables["FromSample"] = fromvar
                    self.outputExperiment.variables["AUC0t"] = AUCvar
                    self.outputExperiment.variables["AUMC0t"] = AUMCvar
                    self.outputExperiment.variables["MRT"] = MRTvar
                    self.outputExperiment.variables["Cmax"] = Cmaxvar
                    self.outputExperiment.variables["Tmax"] = Tmaxvar
                    self.outputExperiment.variables["Cmin"] = Cminvar
                    self.outputExperiment.variables["Tmin"] = Tminvar
                    self.outputExperiment.variables["Cavg"] = Cavgvar

                excelWriteRow(["simulationName", "fromSample", "doseNumber",
                               "AUC [%s]" % strUnit(self.AUCunits),
                               "AUMC [%s]" % strUnit(self.AUMCunits),
                               "Cmin [%s]" % strUnit(self.Cunits.unit),
                               "Cavg [%s]" % strUnit(self.Cunits.unit),
                               "Cmax [%s]" % strUnit(self.Cunits.unit),
                               "Ctau [%s]" % strUnit(self.Cunits.unit),
                               "Tmin [%s]" % strUnit(self.outputExperiment.getTimeUnits().unit),
                               "Tmax [%s]" % strUnit(self.outputExperiment.getTimeUnits().unit),
                               "Ttau [%s]" % strUnit(self.outputExperiment.getTimeUnits().unit)],
                              wb, 1, bold=True)
                wbRow = 2

            # Evaluate AUC, AUMC and MRT in the last full period
            AUClist, AUMClist, Cminlist, Cavglist, Cmaxlist, Ctaulist, Tmaxlist, Tminlist, Ttaulist = self.NCA(self.model.x,y[0])
            for doseNo in range(0,len(AUClist)):
                excelWriteRow(["Simulation_%d"%i, self.fromSample, doseNo, AUClist[doseNo], AUMClist[doseNo],
                               Cminlist[doseNo], Cavglist[doseNo], Cmaxlist[doseNo], Ctaulist[doseNo],
                               Tminlist[doseNo], Tmaxlist[doseNo], Ttaulist[doseNo]], wb, wbRow)
                wbRow+=1

            # Keep results
            for j in range(self.getResponseDimension()):
                simulationsY[i,:,j] = y[j]
            AUCarray[i] = self.AUC0t
            AUMCarray[i] = self.AUMC0t
            MRTarray[i] = self.MRT
            CminArray[i] = self.Cmin
            CmaxArray[i] = self.Cmax
            CavgArray[i] = self.Cavg
            CtauArray[i] = self.Ctau
            TminArray[i] = self.Tmin
            TmaxArray[i] = self.Tmax
            TtauArray[i] = self.Ttau
            fluctuationArray[i] = self.fluctuation
            percentageAccumulationArray[i] = self.percentageAccumulation
            if self.addIndividuals or self.paramsSource!=ProtPKPDODESimulate.PRM_POPULATION:
                self.addSample("Simulation_%d"%i, dosename, simulationsX, y[0])

        # Report NCA statistics
        fhSummary = open(self._getPath("summary.txt"),"w")
        alpha_2 = (100-self.confidenceLevel.get())/2
        limits = np.percentile(AUCarray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"AUC %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.AUCunits),np.mean(AUCarray)))
        limits = np.percentile(AUMCarray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"AUMC %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.AUMCunits),np.mean(AUMCarray)))
        limits = np.percentile(MRTarray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"MRT %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.outputExperiment.getTimeUnits().unit),np.mean(MRTarray)))
        limits = np.percentile(CminArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Cmin %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.Cunits.unit),np.mean(CminArray)))
        limits = np.percentile(CmaxArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Cmax %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.Cunits.unit),np.mean(CmaxArray)))
        limits = np.percentile(CavgArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Cavg %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.Cunits.unit),np.mean(CavgArray)))
        limits = np.percentile(CtauArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Ctau %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.Cunits.unit),np.mean(CtauArray)))
        limits = np.percentile(TminArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Tmin %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.outputExperiment.getTimeUnits().unit),np.mean(TminArray)))
        limits = np.percentile(TmaxArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Tmax %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.outputExperiment.getTimeUnits().unit),np.mean(TmaxArray)))
        limits = np.percentile(TtauArray,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Ttau %f%% confidence interval=[%f,%f] [%s] mean=%f"%(self.confidenceLevel.get(),limits[0],limits[1],strUnit(self.outputExperiment.getTimeUnits().unit),np.mean(TtauArray)))
        aux = fluctuationArray[~np.isnan(fluctuationArray)]
        if len(aux)>0:
            limits = np.percentile(aux,[alpha_2,100-alpha_2])
            self.doublePrint(fhSummary,"Fluctuation %f%% confidence interval=[%f,%f] [%%] mean=%f"%(self.confidenceLevel.get(),limits[0]*100,limits[1]*100,np.mean(aux)*100))
        aux = percentageAccumulationArray[~np.isnan(percentageAccumulationArray)]
        limits = np.percentile(aux,[alpha_2,100-alpha_2])
        self.doublePrint(fhSummary,"Accum(1) %f%% confidence interval=[%f,%f] [%%] mean=%f"%(self.confidenceLevel.get(),limits[0]*100,limits[1]*100,np.mean(aux)*100))
        fhSummary.close()

        # Calculate statistics
        if self.addStats and self.odeSource==0 and self.paramsSource==0:
            if self.paramsSource!=ProtPKPDODESimulate.PRM_USER_DEFINED:
                limits = np.percentile(simulationsY,[alpha_2,100-alpha_2],axis=0)

                print("Lower limit NCA")
                self.NCA(simulationsX,limits[0])
                self.fromSample="LowerLimit"
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
            self.fromSample = "Mean"
            self.addSample("Mean", dosename, simulationsX, mu)

            if self.paramsSource!=ProtPKPDODESimulate.PRM_USER_DEFINED:
                print("Upper limit NCA")
                self.NCA(simulationsX,limits[1])
                self.fromSample="UpperLimit"
                self.addSample("UpperLimit", dosename, simulationsX, limits[1])

        self.outputExperiment.write(self._getPath("experiment.pkpd"))
        wb.save(self._getPath("nca.xlsx"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.outputExperiment)
        if self.odeSource.get()==self.SRC_ODE:
            self._defineSourceRelation(self.inputODE.get(), self.outputExperiment)
            if self.paramsSource==ProtPKPDODESimulate.PRM_POPULATION:
                self._defineSourceRelation(self.inputPopulation.get(), self.outputExperiment)
            elif self.paramsSource==ProtPKPDODESimulate.PRM_POPULATION:
                self._defineSourceRelation(self.inputFitting.get(), self.outputExperiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        msg.append("Dose: %s"%self.doses.get())
        if self.odeSource== ProtPKPDODESimulate.SRC_ODE and self.paramsSource ==  ProtPKPDODESimulate.PRM_POPULATION:
            msg.append("Number of simulations: %d"%self.Nsimulations.get())
        elif  self.odeSource== ProtPKPDODESimulate.SRC_LIST or self.paramsSource ==  ProtPKPDODESimulate.PRM_USER_DEFINED:
            msg.append("Parameters:\n"+self.prmUser.get())
        else:
            msg.append("Parameters from previous fitting")
        msg.append(" ")
        self.addFileContentToMessage(msg,self._getPath("summary.txt"))
        return msg

    def _validate(self):
        msg=[]
        if self.odeSource == self.SRC_ODE and \
            self.paramsSource == ProtPKPDODESimulate.PRM_POPULATION and \
            not self.inputPopulation.get().fnFitting.get().endswith("bootstrapPopulation.pkpd"):
            msg.append("Population must be a bootstrap sample")
        return msg
