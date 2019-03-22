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

from itertools import izip
import numpy as np
import os

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDODEModel
from pkpd.biopharmaceutics import DrugSource, createDeltaDose, createVia


class PKPDLiver(PKPDODEModel):
    def F(self, t, y):
        Vsys=self.parameters[0]
        ClNH=self.parameters[1]
        fb=self.parameters[2]
        Kp=self.parameters[3]
        Vinlet=self.parameters[4]
        Vliver=self.parameters[5]
        Qh=self.parameters[6]
        Clint=self.parameters[7]
        Iliver=y[0]
        Iinlet=y[1]
        Isys=y[2]

        retval = np.array([1/Vliver*(Qh*Iinlet-(Qh+fb*Clint)*Iliver/Kp),
                          Qh/Vinlet*(Isys-Iinlet),
                          1/Vsys*(Qh*(Iliver/Kp-Isys)-ClNH*Isys)],np.double)
        return retval

    def G(self, t, dD):
        Vinlet=self.parameters[4]
        return np.array([0.0,dD/Vinlet,0.0],np.double)

    def imposeConstraints(self, yt):
        if yt[0]<0:
            yt[0]=0
        if yt[1]<0:
            yt[1]=0
        if yt[2]<0:
            yt[2]=0

    def getResponseDimension(self):
        return 3

    def getStateDimension(self):
        return 3


class PKPDLiverEV1():
    def __init__(self):
        self.drugSource = DrugSource()
        self.model = PKPDLiver()
        self.model.drugSource = self.drugSource
        self.via = createVia("Oral; ev1")

    def setTimeRange(self,tF):
        self.model.t0 = 0
        self.model.tF = tF*60

    def setDose(self,doseAmount):
        self.drugSource.setDoses([createDeltaDose(doseAmount,self.via,0,"mg")],self.model.t0,self.model.tF)
        self.NparametersSource = len(self.drugSource.getParameterNames())

    def simulate(self,params,t):
        self.drugSource.setParameters(params[0:self.NparametersSource])
        return self.model.forwardModel(params[self.NparametersSource:],t)


class ProtPKPDSimulateLiverFlow(ProtPKPD):
    """ Simulate the concentration of a compound (typically an enzyme inhibitor) at liver.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'simulate liver flow'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form, fullForm=True):
        form.addSection('Input')
        form.addParam('tF', params.FloatParam, default=8, label='Max. simulation time [h]')

        group = form.addGroup("Absorption")
        group.addParam("weight", params.StringParam, default=70, label="Weight [kg]")
        group.addParam("dose", params.StringParam, default=1, label="Dose [mg/kg]")
        group.addParam("Fa", params.StringParam, default=1, label="Fraction absorbed (0-1)")
        group.addParam("ka", params.StringParam, default=0.05, label="1st order absorption rate [1/min]", help="Typically between 0.0003 and 0.1")

        group = form.addGroup("Central compartment")
        group.addParam("Vsys", params.StringParam, default=200, label="Distribution volume [mL/kg]")
        group.addParam("ClNH", params.StringParam, default=1, label="Non hepatic clearance [mL/min/kg]", help="Typically between 0.6 and 600")
        group.addParam("fb", params.StringParam, default=1, label="Unbound fraction in blood [0-1]")

        group = form.addGroup("Liver")
        group.addParam("Kp", params.StringParam, default=1, label="Liver-to-blood concentration ratio")
        group.addParam("Vinlet", params.StringParam, default=1, label="Liver inlet volume [mL/kg]")
        group.addParam("Vliver", params.StringParam, default=40, label="Liver volume [mL/kg]")
        group.addParam("Qh", params.StringParam, default=30, label="Liver blood flow [mL/min/kg]")
        group.addParam("Clint", params.StringParam, default=10, label="Intrinsic metabolic clearance [mL/min/kg]", help="Typically between 3 and 300")


    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulate')

    #--------------------------- STEPS functions --------------------------------------------
    def parseList(self, strList):
        return [float(v) for v in strList.split(' ')]

    def runSimulate(self):
        model = PKPDLiverEV1()
        model.setTimeRange(self.tF.get())
        t = np.arange(0.0, self.tF.get()*60, 1)

        I = []
        legends = []
        weightList = self.parseList(self.weight.get())
        doseList = self.parseList(self.dose.get())
        FaList = self.parseList(self.Fa.get())
        kaList = self.parseList(self.ka.get())
        VsysList = self.parseList(self.Vsys.get())
        ClNHList = self.parseList(self.ClNH.get())
        fbList = self.parseList(self.fb.get())
        KpList = self.parseList(self.Kp.get())
        VinletList = self.parseList(self.Vinlet.get())
        VliverList = self.parseList(self.Vliver.get())
        QhList = self.parseList(self.Qh.get())
        ClintList = self.parseList(self.Clint.get())
        for weight in weightList:
            for dose in doseList:
                for Fa in FaList:
                    for ka in kaList:
                        for Vsys in VsysList:
                            for ClNH in ClNHList:
                                for fb in fbList:
                                    for Kp in KpList:
                                        for Vinlet in VinletList:
                                            for Vliver in VliverList:
                                                for Qh in QhList:
                                                    for Clint in ClintList:
                                                        legend="Weight=%f Dose=%f Fa=%f ka=%f Vsys=%f ClNH=%f fb=%f Kp=%f Vinlet=%f Vliver=%f Qh=%f Clint=%f"%\
                                                               (weight, dose, Fa, ka, Vsys, ClNH, fb, Kp, Vinlet, Vliver, Qh, Clint)
                                                        print("Simulating %s"%legend)
                                                        legends.append(legend)
                                                        model.setDose(Fa*dose*weight)
                                                        Ii = model.simulate([ka,Vsys*weight,ClNH*weight,fb,Kp,
                                                                             Vinlet*weight,Vliver*weight,Qh*weight*model.model.deltaT,
                                                                             Clint*weight],t)
                                                        I.append((t,Ii))

        if len(I)>0:
            fh=open(self._getPath("profiles.txt"),'w')
            fhSummary=open(self._getPath("summary.txt"),"w")
            for toPlot, legend in izip(I,legends):
                fh.write("SimulateLiver::"+legend+"\n")
                fhSummary.write("Simulated %s\n"%(legend))
                t, I3 = toPlot
                for n in range(t.size):
                    fh.write("%f %f %f %f\n"%(t[n],I3[0][n],I3[1][n],I3[2][n]))
                fh.write("\n")
            fh.close()
            fhSummary.close()

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        if os.path.exists(self._getPath("summary.txt")):
            fh=open(self._getPath("summary.txt"))
            for line in fh:
                msg.append(line.strip())
            fh.close()
        return msg

    def _citations(self):
        retval = ['Kanamitsu2000']
        return retval
