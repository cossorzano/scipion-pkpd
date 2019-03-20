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
from pyworkflow.protocol.constants import LEVEL_ADVANCED

class ProtPKPDSimulateDrugInteractions(ProtPKPD):
    """ Simulate drug interactions as recommended in EMA CHMP/EWP/560/95 \n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'simulate drug interactions'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form, fullForm=True):
        form.addSection('Basic models Liver')
        fromToLiver = form.addLine('Inhibitor range [I] at liver')
        fromToLiver.addParam('I0Liver', params.FloatParam, default=0, label='Min (ng/mL)')
        fromToLiver.addParam('IFLiver', params.FloatParam, default=10, label='Max (ng/mL)')
        form.addParam("MWLiver", params.FloatParam, default=1, label="Molecular weight (g/mol)")

        form.addParam("doReversibleLiver", params.BooleanParam, default=False, label="Reversible inhibition",
                      help="Investigational drug likely to be a reversible inhibitor if R1=1+[I]/Ki>1.02")
        form.addParam('KiReversibleLiver', params.StringParam, default="5", label='Inhibition constant (Ki [uM])', condition="doReversibleLiver",
                      help="Ki is the in vitro unbound reversible inhibition constant. Several constants can be given separated by space, e.g., 5 10")

        form.addParam("doTimeDependentLiver", params.BooleanParam, default=False, label="Time dependent inhibition",
                      help="Investigational drug likely to be a time dependent inhibitor if R2=1+kinact/kdeg*[I]/([I]+Ki)>1.25")
        form.addParam("kinactLiver",params.StringParam, default="0.01", label='Max. inactivation rate (kinact [min^-1])', condition="doTimeDependentLiver",
                      help="Maximal inactivation rate. Several constants can be given separated by space, e.g., 0.01 0.005")
        form.addParam("kdegLiver",params.StringParam, default="0.008", label='Apparent first order degradation rate (kdeg [min^-1])', condition="doTimeDependentLiver",
                      help="kdeg is the apparent first order degradation rate constant of the affected enzyme. Several constants can be given separated by space, e.g., 0.01 0.005")
        form.addParam('KiTimeLiver', params.StringParam, default="5", label='Inhibition constant (Ki [uM])', condition="doTimeDependentLiver",
                      help="KI is the inhibitor concentration which yields 50% of the maximum inactivation rate. Several constants can be given separated by space, e.g., 5 10")

        form.addParam("doInductionLiver", params.BooleanParam, default=False, label="Induction",
                      help="Investigational drug likely to be a time dependent inhibitor if R2=1+kinact/kdeg*[I]/([I]+Ki)>1.25")
        form.addParam("dLiver",params.StringParam, default="1", label='Scaling factor', condition="doInductionLiver",
                      help="Several constants can be given separated by space, e.g., 1 1.05")
        form.addParam("EmaxLiver",params.StringParam, default="1", label='Max. Induction effect (Emax)', condition="doInductionLiver",
                      help="Several constants can be given separated by space, e.g., 1 2")
        form.addParam('EC50Liver', params.StringParam, default="5", label='Half Max. Effect Conc (EC50 [uM])', condition="doInductionLiver",
                      help="EC50 is the concentration causing half maximal effect. Several constants can be given separated by space, e.g., 5 6")

        form.addSection('Basic models Gut')
        fromToGut1 = form.addLine('Basic model Gut Oral dose')
        fromToGut1.addParam('D0Gut', params.FloatParam, default=0, label='Min (mg)')
        fromToGut1.addParam('DFGut', params.FloatParam, default=11, label='Max (mg)')
        form.addParam("MWGut", params.FloatParam, default=1, label="Molecular weight (g/mol)")

        form.addParam("doReversibleGut", params.BooleanParam, default=False, label="Reversible inhibition",
                      help="Investigational drug likely to be a reversible inhibitor if R1=1+[I]/Ki>1.02")
        form.addParam('KiReversibleGut', params.StringParam, default="5", label='Inhibition constant (Ki [uM])', condition="doReversibleGut",
                      help="Ki is the in vitro unbound reversible inhibition constant. Several constants can be given separated by space, e.g., 5 10")

        form.addParam("doTimeDependentGut", params.BooleanParam, default=False, label="Time dependent inhibition",
                      help="Investigational drug likely to be a time dependent inhibitor if R2=1+kinact/kdeg*[I]/([I]+Ki)>1.25")
        form.addParam("kinactGut",params.StringParam, default="0.01", label='Max. inactivation rate (kinact [min^-1])', condition="doTimeDependentGut",
                      help="Maximal inactivation rate. Several constants can be given separated by space, e.g., 0.01 0.005")
        form.addParam("kdegGut",params.StringParam, default="0.008", label='Apparent first order degradation rate (kdeg [min^-1])', condition="doTimeDependentGut",
                      help="kdeg is the apparent first order degradation rate constant of the affected enzyme. Several constants can be given separated by space, e.g., 0.01 0.005")
        form.addParam('KiTimeGut', params.StringParam, default="5", label='Inhibition constant (Ki [uM])', condition="doTimeDependentGut",
                      help="KI is the inhibitor concentration which yields 50% of the maximum inactivation rate. Several constants can be given separated by space, e.g., 5 10")

        form.addParam("doInductionGut", params.BooleanParam, default=False, label="Induction",
                      help="Investigational drug likely to be a time dependent inhibitor if R2=1+kinact/kdeg*[I]/([I]+Ki)>1.25")
        form.addParam("dGut",params.StringParam, default="1", label='Scaling factor', condition="doInductionGut",
                      help="Several constants can be given separated by space, e.g., 1 1.05")
        form.addParam("EmaxGut",params.StringParam, default="1", label='Max. Induction effect (Emax)', condition="doInductionGut",
                      help="Several constants can be given separated by space, e.g., 1 2")
        form.addParam('EC50Gut', params.StringParam, default="5", label='Half Max. Effect Conc (EC50 [uM])', condition="doInductionGut",
                      help="EC50 is the concentration causing half maximal effect. Several constants can be given separated by space, e.g., 5 6")

        form.addSection('Static model')
        form.addParam("doStatic", params.BooleanParam, default=False, label="Static",
                      help="Mechanistic static model, Fahmi2009")
        form.addParam('KiStatic', params.StringParam, default="5", label='Inhibition constant (Ki [uM])', condition="doStatic",
                      help="Ki is the in vitro unbound reversible inhibition constant. Several constants can be given separated by space, e.g., 5 10")
        form.addParam("EmaxStatic",params.StringParam, default="1", label='Max. Induction effect (Emax)', condition="doStatic",
                      help="Several constants can be given separated by space, e.g., 1 2")
        form.addParam('EC50Static', params.StringParam, default="5", label='Half Max. Effect Conc (EC50 [uM])', condition="doStatic",
                      help="EC50 is the concentration causing half maximal effect. Several constants can be given separated by space, e.g., 5 6")
        form.addParam("kinactStatic",params.StringParam, default="0.01", label='Max. inactivation rate (kinact [min^-1])', condition="doStatic",
                      help="Maximal inactivation rate. Several constants can be given separated by space, e.g., 0.01 0.005")
        form.addParam("dStatic",params.StringParam, default="1", label='Scaling factor', condition="doStatic",
                      help="Scaling factor determined with linear regression of the control data set")
        form.addParam("MWStatic", params.FloatParam, default=1, label="Molecular weight (g/mol)", condition="doStatic")

        group1 = form.addGroup("Physiological", condition="doStatic")
        group1.addParam("doPhysiological", params.BooleanParam, default=False, label="Physiological model",
                      help="The physiological model estimates the concentration at liver and gut from "
                      "the input dose, the fraction available, the absorption rate, and the liver and gut blood flow."
                      "The alternative to the physiological model needs a direct estimation of the gut and liver concentration of the inhibitor")
        fromToGut2 = group1.addLine('Static model Oral dose', condition="doPhysiological and doStatic")
        fromToGut2.addParam('D0Phys', params.FloatParam, default=0, label='Min (mg)', condition="doPhysiological and doStatic")
        fromToGut2.addParam('DFPhys', params.FloatParam, default=10, label='Max (mg)', condition="doPhysiological and doStatic")
        group1.addParam("Fa", params.FloatParam, default=1, label="Fraction absorbed", condition="doPhysiological and doStatic",
                       help="Fa is the fraction absorbed after oral administration (a value of 1 can be used if data are not available)")
        group1.addParam("ka", params.FloatParam, default=0.1, label="1st order absorption rate", condition="doPhysiological and doStatic",
                       help="ka is the first order absorption rate constant in vivo and a value of 0.1 min^-1 can be used if data are not available.")
        group1.addParam("Qen", params.FloatParam, default=18, label="Blood flow at enterocytes [L/h/70kg", condition="doPhysiological and doStatic",
                       help="Qen is blood flow through the enterocytes (e.g., 18L/hr/70kg taken from Yang2007a)")
        group1.addParam("fub", params.FloatParam, default=1, label="Fraction unbound in plasma", condition="doPhysiological and doStatic")
        group1.addParam("Imaxb", params.FloatParam, default=1, label="Maximal total [I] conc.", condition="doPhysiological and doStatic",
                       help="[I]max,b is the maximal total (free and bound) inhibitor concentration in the blood at steady state")
        group1.addParam("Qh", params.FloatParam, default=97, label="Blood flow at hepatocytes [L/h/70kg", condition="doPhysiological and doStatic",
                       help="Qh is hepatic blood flow (e.g., 97L/hr/70kg taken from Yang et al., 2007b)")

        group2 = form.addGroup("Liver", condition="doStatic")
        group2.addParam("doStaticLiver", params.BooleanParam, default=False, label="Liver")
        group2.addParam("fm",params.FloatParam, default=0.1, label='fm', condition="doStaticLiver",
                      help="fm is the fraction of systemic clearance of the substrate mediated by the enzyme that is subject to inhibition/induction")
        fromToH = group2.addLine('Inhibitor range Liver [Ih]', condition="doStaticLiver and not doPhysiological", help='[Ih] is [I]liver = Molar Dose/250mL.')
        fromToH.addParam('Ih0', params.FloatParam, default=0, condition="doStaticLiver and not doPhysiological", label='Min (uM)')
        fromToH.addParam('IhF', params.FloatParam, default=10, condition="doStaticLiver and not doPhysiological", label='Max (uM)')
        group2.addParam("kdeghStatic",params.StringParam, default="0.008", label='Apparent first order degradation rate Liver (kdeg [min^-1])', condition="doStaticLiver",
                      help="kdeg,h is the apparent first order degradation rate constant of the affected enzyme at the liver. Several constants can be given separated by space, e.g., 0.01 0.005")

        group3 = form.addGroup("Gut", condition="doStatic")
        group3.addParam("doStaticGut", params.BooleanParam, default=False, label="Gut")
        group3.addParam("fg",params.FloatParam, default=0.1, label='fg', condition="doStaticGut",
                      help="fg is the fraction available after intestinal metabolism")
        fromToG = group3.addLine('Inhibitor range Gut [Ig]', condition="doStaticGut and not doPhysiological", help='[Ig] is [I]gut = Molar Dose/250mL.')
        fromToG.addParam('Ig0', params.FloatParam, default=0, condition="doStaticGut and not doPhysiological", label='Min (uM)')
        fromToG.addParam('IgF', params.FloatParam, default=10, condition="doStaticGut and not doPhysiological", label='Max (uM)')
        group3.addParam("kdeggStatic",params.StringParam, default="0.008", label='Apparent first order degradation rate Gut (kdeg [min^-1])', condition="doStaticGut",
                      help="kdeg,g is the apparent first order degradation rate constant of the affected enzyme at the gut. Several constants can be given separated by space, e.g., 0.01 0.005")

        form.addSection('Transporters')
        form.addParam("doTransporterGut", params.BooleanParam, default=False, label="Gut transporter")
        fromToGut3 = form.addLine('Transporters Oral dose', condition="doTransporterGut")
        fromToGut3.addParam('D0TransporterGut', params.FloatParam, default=0, label='Min (mg)')
        fromToGut3.addParam('DFTransporterGut', params.FloatParam, default=10, label='Max (mg)')
        form.addParam("MWTransporterGut", params.FloatParam, default=1, label="Molecular weight (g/mol)", condition="doTransporterGut")
        form.addParam('KiTransporterGut', params.StringParam, default="5", label='Inhibition constant (Ki [uM])', condition="doTransporterGut",
                      help="Ki is the in vitro unbound reversible inhibition constant. Several constants can be given separated by space, e.g., 5 10")

        form.addParam("doTransporterLiver", params.BooleanParam, default=False, label="Liver transporter")
        fromToGut4 = form.addLine('Unbound hepatic inlet [I]', condition="doTransporterLiver")
        fromToGut4.addParam('I0TransporterLiver', params.FloatParam, default=0, label='Min (uM)')
        fromToGut4.addParam('IFTransporterLiver', params.FloatParam, default=10, label='Max (uM)')
        form.addParam('KiTransporterLiver', params.StringParam, default="5", label='Inhibition constant (Ki [uM])', condition="doTransporterLiver",
                      help="Ki is the in vitro unbound reversible inhibition constant. Several constants can be given separated by space, e.g., 5 10")

        form.addParam("doTransporterRenal", params.BooleanParam, default=False, label="Renal transporter")
        fromToGut5 = form.addLine('Unbound [Imax]', condition="doTransporterRenal")
        fromToGut5.addParam('I0TransporterRenal', params.FloatParam, default=0, label='Min (uM)')
        fromToGut5.addParam('IFTransporterRenal', params.FloatParam, default=10, label='Max (uM)')
        form.addParam('KiTransporterRenal', params.StringParam, default="5", label='Inhibition constant (Ki [uM])', condition="doTransporterRenal",
                      help="Ki is the in vitro unbound reversible inhibition constant. Several constants can be given separated by space, e.g., 5 10")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulate')

    #--------------------------- STEPS functions --------------------------------------------
    def parseList(self, strList):
        return [float(v) for v in strList.split(' ')]

    def runSimulate(self):
        R = []
        Rlegends = []
        Type = []
        if self.doReversibleLiver:
            I = np.arange(self.I0Liver.get(), self.IFLiver.get(), (self.IFLiver.get()-self.I0Liver.get())/100) # I [ng/mL]
            I /= self.MWLiver.get() # I [uM]
            KiList = self.parseList(self.KiReversibleLiver.get())
            for Ki in KiList:
                legend="Liver Rev. Inh. Ki=%f [uM]"%Ki
                print("Simulating %s"%legend)
                R.append((I,1+I/Ki))
                Rlegends.append(legend)
                Type.append('ReversibleLiver')

        if self.doReversibleGut:
            D = np.arange(self.D0Gut.get(), self.DFGut.get(), (self.DFGut.get()-self.D0Gut.get())/100)
            I = D/(250*self.MWGut.get())*1e6 # I [uM]
            KiList = self.parseList(self.KiReversibleGut.get())
            for Ki in KiList:
                legend="Gut Rev. Inh. Ki=%f [uM]"%Ki
                print("Simulating %s"%legend)
                R.append((I,1+I/Ki))
                Rlegends.append(legend)
                Type.append('ReversibleGut')

        if self.doTimeDependentLiver:
            I = np.arange(self.I0Liver.get(), self.IFLiver.get(), (self.IFLiver.get()-self.I0Liver.get())/100) # I [ng/mL]
            I /= self.MWLiver.get() # I [uM]
            KiList = self.parseList(self.KiReversibleLiver.get())
            kdegList = self.parseList(self.kdegLiver.get())
            kinactList = self.parseList(self.kinactLiver.get())
            for Ki in KiList:
                for kdeg in kdegList:
                    for kinact in kinactList:
                        legend="Liver Time Dep. Inh. Ki=%f [uM], kdeg=%f [min^-1], kinact=%f [min^-1]"%(Ki,kdeg,kinact)
                        print("Simulating %s"%legend)
                        R.append((I,1+kinact/kdeg*I/(Ki+I)))
                        Rlegends.append(legend)
                        Type.append('TimeDependentLiver')

        if self.doTimeDependentGut:
            D = np.arange(self.D0Gut.get(), self.DFGut.get(), (self.DFGut.get()-self.D0Gut.get())/100)
            I = D/(250*self.MWGut.get())*1e6 # I [uM]
            KiList = self.parseList(self.KiReversibleGut.get())
            kdegList = self.parseList(self.kdegGut.get())
            kinactList = self.parseList(self.kinactGut.get())
            for Ki in KiList:
                for kdeg in kdegList:
                    for kinact in kinactList:
                        legend="Gut Time Dep. Inh. Ki=%f [uM], kdeg=%f [min^-1], kinact=%f [min^-1]"%(Ki,kdeg,kinact)
                        print("Simulating %s"%legend)
                        R.append((I,1+kinact/kdeg*I/(Ki+I)))
                        Rlegends.append(legend)
                        Type.append('TimeDependentGut')

        if self.doInductionLiver:
            I = np.arange(self.I0Liver.get(), self.IFLiver.get(), (self.IFLiver.get()-self.I0Liver.get())/100) # I [ng/mL]
            I /= self.MWLiver.get() # I [uM]
            EC50List = self.parseList(self.EC50Liver.get())
            EmaxList = self.parseList(self.EmaxLiver.get())
            dList = self.parseList(self.dLiver.get())
            for EC50 in EC50List:
                for Emax in EmaxList:
                    for d in dList:
                        legend="Liver Induction EC50=%f [uM], Emax=%f, d=%f"%(EC50,Emax,d)
                        print("Simulating %s"%legend)
                        R.append((I,1/(1+d*Emax*I/(EC50+I))))
                        Rlegends.append(legend)
                        Type.append('InductionLiver')

        if self.doInductionGut:
            D = np.arange(self.D0Gut.get(), self.DFGut.get(), (self.DFGut.get()-self.D0Gut.get())/100)
            I = D/(250*self.MWGut.get())*1e6 # I [uM]
            EC50List = self.parseList(self.EC50Gut.get())
            EmaxList = self.parseList(self.EmaxGut.get())
            dList = self.parseList(self.dGut.get())
            for EC50 in EC50List:
                for Emax in EmaxList:
                    for d in dList:
                        legend="Induction EC50=%f [uM], Emax=%f, d=%f"%(EC50,Emax,d)
                        print("Simulating %s"%legend)
                        R.append((I,1/(1+d*Emax*I/(EC50+I))))
                        Rlegends.append(legend)
                        Type.append('InductionGut')

        if self.doStatic:
            KiList = self.parseList(self.KiStatic.get())
            EmaxList = self.parseList(self.EmaxStatic.get())
            EC50List = self.parseList(self.EC50Static.get())
            kinactList = self.parseList(self.kinactStatic.get())
            dList = self.parseList(self.dStatic.get())
            if self.doStaticLiver:
                kdeghList = self.parseList(self.kdeghStatic.get())
                fm=self.fm.get()
                if self.doPhysiological:
                    D = np.arange(self.D0Phys.get(), self.DFPhys.get(), (self.DFPhys.get()-self.D0Phys.get())/100)
                    Ih = self.fub.get()*(self.Imaxb.get()+self.Fa.get()*self.ka.get()*D/self.Qh.get())
                else:
                    Ih = np.arange(self.Ih0.get(), self.IhF.get(), (self.IhF.get()-self.Ih0.get())/100)
                Ih/= self.MWStatic.get()
                for Ki in KiList:
                    for Emax in EmaxList:
                        for EC50 in EC50List:
                            for kinact in kinactList:
                                for d in dList:
                                    for kdegh in kdeghList:
                                        Ah = kdegh/(kdegh+Ih*kinact/(Ih+Ki))
                                        Bh = 1+d*Emax*Ih/(Ih+EC50)
                                        Ch = 1/(1+Ih/Ki)

                                        legend="Static Liver Ki=%f [uM], EC50=%f [uM], Emax=%f, kinact=%f [min^-1], d=%f, kdegh=%f [min^-1]"%\
                                               (Ki,EC50,Emax,kinact,d,kdegh)
                                        print("Simulating %s"%legend)
                                        R.append((Ih,1/(Ah*Bh*Ch*fm+(1-fm))))
                                        Rlegends.append(legend)
                                        Type.append('StaticLiver')

            if self.doStaticGut:
                kdeggList = self.parseList(self.kdeggStatic.get())
                if self.doPhysiological:
                    D = np.arange(self.D0Phys.get(), self.DFPhys.get(), (self.DFPhys.get()-self.D0Phys.get())/100)
                    Ig = self.Fa.get()*self.ka.get()*D/self.Qen.get()
                else:
                    Ig = np.arange(self.Ig0.get(), self.IgF.get(), (self.IgF.get()-self.Ig0.get())/100)
                Ig/= self.MWStatic.get()
                fg=self.fg.get()
                for Ki in KiList:
                    for Emax in EmaxList:
                        for EC50 in EC50List:
                            for kinact in kinactList:
                                for d in dList:
                                    for kdegg in kdeggList:
                                        Ag = kdegg/(kdegg+Ig*kinact/(Ig+Ki))
                                        Bg = 1+d*Emax*Ig/(Ig+EC50)
                                        Cg = 1/(1+Ig/Ki)

                                        legend="Static Gut Ki=%f [uM], EC50=%f [uM], Emax=%f, kinact=%f [min^-1], d=%f, kdegg=%f [min^-1]"%\
                                               (Ki,EC50,Emax,kinact,d,kdegg)
                                        print("Simulating %s"%legend)
                                        R.append((Ig,1/(Ag*Bg*Cg*fg+(1-fg))))
                                        Rlegends.append(legend)
                                        Type.append('StaticGut')

        if self.doTransporterGut:
            D = np.arange(self.D0TransporterGut.get(), self.DFTransporterGut.get(), (self.DFTransporterGut.get()-self.D0TransporterGut.get())/100)
            I = D/(250*self.MWTransporterGut.get())*1e6 # I [uM]
            KiList = self.parseList(self.KiTransporterGut.get())
            for Ki in KiList:
                legend="Gut Transporter Ki=%f [uM]"%Ki
                print("Simulating %s"%legend)
                R.append((I,1+I/Ki))
                Rlegends.append(legend)
                Type.append('TransporterGut')

        if self.doTransporterLiver:
            I = np.arange(self.I0TransporterLiver.get(), self.IFTransporterLiver.get(), (self.IFTransporterLiver.get()-self.I0TransporterLiver.get())/100)
            KiList = self.parseList(self.KiTransporterLiver.get())
            for Ki in KiList:
                legend="Liver Transporter Ki=%f [uM]"%Ki
                print("Simulating %s"%legend)
                R.append((I,1+I/Ki))
                Rlegends.append(legend)
                Type.append('TransporterLiver')

        if self.doTransporterRenal:
            I = np.arange(self.I0TransporterRenal.get(), self.IFTransporterRenal.get(), (self.IFTransporterRenal.get()-self.I0TransporterRenal.get())/100)
            KiList = self.parseList(self.KiTransporterRenal.get())
            for Ki in KiList:
                legend="Renal Transporter Ki=%f [uM]"%Ki
                print("Simulating %s"%legend)
                R.append((I,1+I/Ki))
                Rlegends.append(legend)
                Type.append('TransporterRenal')

        if len(R)>0:
            fh=open(self._getPath("profiles.txt"),'w')
            fhSummary=open(self._getPath("summary.txt"),"w")
            for toPlot, legend, plotType in izip(R,Rlegends,Type):
                fh.write(plotType+"::"+legend+"\n")
                fhSummary.write("Simulated %s:: %s\n"%(plotType,legend))
                if len(toPlot)==2:
                    I, Ri = toPlot
                    for n in range(I.size):
                        fh.write("%f %f\n"%(I[n],Ri[n]))
                else:
                    Ig, Ih, Ri = toPlot
                    for n in range(Ig.size):
                        fh.write("%f %f %f\n"%(Ig[n],Ih[n],Ri[n]))
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
        retval = ['CHMPEWP56095','Fahmi2009']
        # if self.doPhysiological:
        #     retval+=['Yang2007a','Yang2007b','Rostami2004']
        return retval

    def _methods(self):
        retval = []
        if self.doTransporterGut:
            retval.append("We studied the inhibition of an intestinal transporter with a dose between %f and %f [mg] "
                          "of a compound whose molecular weight was %f [g/mol]. The inhibition constant of the transporter "
                          "was assumed to be %s [uM]"%(self.D0TransporterGut,self.DFTransporterGut,self.MWTransporterGut,self.KiTransporterGut))
        if self.doTransporterLiver:
            retval.append("We studied the inhibition of a hepatic transporter with a unbound hepatic concentration of the inhibitor at the liver inlet between %f and %f [uM]. "
                          "The inhibition constant of the transporter "
                          "was assumed to be %s [uM]"%(self.I0TransporterLiver,self.IFTransporterLiver,self.KiTransporterLiver))
        if self.doTransporterRenal:
            retval.append("We studied the inhibition of a renal transporter with a unbound hepatic concentration of the inhibitor at the liver inlet between %f and %f [uM]. "
                          "The inhibition constant of the transporter "
                          "was assumed to be %s [uM]"%(self.I0TransporterRenal,self.IFTransporterRenal,self.KiTransporterRenal))
        return retval