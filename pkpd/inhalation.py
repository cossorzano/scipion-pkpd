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
Classes for inhalation PK simulations.

See Hartung and Borghardt. A mechanistic framework for a priori pharmacokinetic
predictions of orally inhaled drugs. PLOS Computational Biology, 16: e1008466 (2020)
"""
import math
import numpy as np
from pwem.objects import EMObject
from pyworkflow.object import String

class PKPhysiologyLungParameters(EMObject):
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.fnPhys = String()

    def write(self, fnOut):
        fh=open(fnOut,"w")
        fh.write("%f # cardiac output [mL/min]\n"%self.Qco)
        fh.write("%f # lung tissue weight, without blood [g]\n"%self.OWlun)
        fh.write("%f # alveolar fraction of cardiac output\n"%self.falvCO)
        fh.write("%f # alveolar ELF volume [mL]\n"%self.Velf_alv)
        fh.write("%f # alveolar surface area [cm2]\n"%self.Surf_alv)
        fh.write("%f # bronchial fraction of cardiac output\n"%self.fbrCO)
        fh.write("%f # bronchial fraction of lung tissue volume\n"%self.fbrVlun)
        fh.write("%f # bronchial ELF heights for interpolation, trachea [cm]\n" % self.helf_trach)
        fh.write("%f # bronchial ELF heights for interpolation, terminal bronchioles [cm]\n" % self.helf_termbr)
        fh.write("%s # trachea length [cm]\n" % self.tracheaLength)
        fh.write("%s # trachea diameter [cm]\n" % self.tracheaDiameter)
        fh.write("%s # main bronchi length [cm]\n" % self.bronchi1Length)
        fh.write("%s # main bronchi diameter [cm]\n" % self.bronchi1Diameter)
        fh.write("%s # bronchi length [cm]\n" % self.bronchi2Length)
        fh.write("%s # bronchi diameter [cm]\n" % self.bronchi2Diameter)
        fh.write("%s # bronchiole length [cm]\n" % self.bronchi3Length)
        fh.write("%s # bronchiole diameter [cm]\n" % self.bronchi3Diameter)
        fh.write("%s # terminal bronchiole length [cm]\n" % self.bronchi4Length)
        fh.write("%s # terminal bronchiole diameter [cm]\n" % self.bronchi4Diameter)

        fh.close()
        self.fnPhys.set(fnOut)

    def read(self, fnIn):
        fh=open(fnIn)
        self.Qco=float(fh.readline().split()[0])
        self.OWlun=float(fh.readline().split()[0])
        self.falvCO=float(fh.readline().split()[0])
        self.Velf_alv=float(fh.readline().split()[0])
        self.Surf_alv=float(fh.readline().split()[0])
        self.fbrCO=float(fh.readline().split()[0])
        self.fbrVlun=float(fh.readline().split()[0])
        self.helf_trach=float(fh.readline().split()[0])
        self.helf_termbr=float(fh.readline().split()[0])
        self.tracheaLength=fh.readline().split()[0]
        self.tracheaDiameter=fh.readline().split()[0]
        self.bronchi1Length=fh.readline().split()[0]
        self.bronchi1Diameter=fh.readline().split()[0]
        self.bronchi2Length=fh.readline().split()[0]
        self.bronchi2Diameter=fh.readline().split()[0]
        self.bronchi3Length=fh.readline().split()[0]
        self.bronchi3Diameter=fh.readline().split()[0]
        self.bronchi4Length=fh.readline().split()[0]
        self.bronchi4Diameter=fh.readline().split()[0]
        fh.close()
        self.fnPhys.set(fnIn)

    def getSystemic(self):
        # Hartung2020_MATLAB/physiol_subst/physiology_systemic.m
        data = {}
        data['Qco']=self.Qco
        data['OWlung']=self.OWlun
        return data

    def getBronchial(self):
        # Hartung2020_MATLAB/physiol_subst/physiology_bronchial.m
        def listSplit(valuesStr):
            return [float(x) for x in valuesStr.split(',')]

        tracheaLength    = float(self.tracheaLength)
        tracheaDiameter  = float(self.tracheaDiameter)
        bronchi1Length   = float(self.bronchi1Length)
        bronchi1Diameter = float(self.bronchi1Diameter)
        bronchi2Length   = listSplit(self.bronchi2Length)
        bronchi2Diameter = listSplit(self.bronchi2Diameter)
        bronchi3Length   = listSplit(self.bronchi3Length)
        bronchi3Diameter = listSplit(self.bronchi3Diameter)
        bronchi4Length   = float(self.bronchi4Length)
        bronchi4Diameter = float(self.bronchi4Diameter)

        segmentLengths   = [tracheaLength,   bronchi1Length]  +bronchi2Length  +bronchi3Length  +[bronchi4Length]
        segmentDiameters = [tracheaDiameter, bronchi1Diameter]+bronchi2Diameter+bronchi3Diameter+[bronchi4Diameter]
        segmentType      = ['T','B1'] + ['B2']*len(bronchi2Length) + ['B3']*len(bronchi3Length) + ['B4']

        length_cm = np.asarray(segmentLengths)
        diam_cm = np.asarray(segmentDiameters)
        generation = np.arange(0,len(segmentLengths))+1
        number = np.power(2.0,generation-1)

        xArea_cm2 = math.pi*np.multiply(np.power(diam_cm/2,2.0),number) # cm2
        vol_cm3 = np.multiply(xArea_cm2,length_cm) # cm3
        start_cm = np.insert(np.cumsum(segmentLengths[:-1]),0,0,axis=0)
        end_cm =  np.cumsum(segmentLengths)
        pos =  0.5*(start_cm+end_cm)
        surf_cm2 = np.multiply(np.multiply(length_cm,(math.pi * diam_cm)),number)
        h_elf_cm = np.interp(pos, [0, np.sum(length_cm)], [self.helf_trach, self.helf_termbr])
        elf_cm3 = math.pi/4 * (np.power(diam_cm,2)-np.power((diam_cm-2*h_elf_cm),2))
        elf_cm3 = np.multiply(np.multiply(elf_cm3,length_cm),number)

        Vtis_ctr = self.OWlun*self.fbrVlun; # central lung tissue volume
        voltis_cm3 = elf_cm3 * Vtis_ctr / sum(elf_cm3);

        data = {}
        data['generation']=generation
        data['type']=segmentType
        data['number']=number
        data['length_cm']=length_cm
        data['diam_cm']=diam_cm
        data['xArea_cm2']=xArea_cm2
        data['vol_cm3']=vol_cm3
        data['start_cm']=start_cm
        data['end_cm']=end_cm
        data['pos']=pos
        data['surf_cm2']=surf_cm2

        data['h_elf_trach'] = self.helf_trach
        data['h_elf_termbr'] = self.helf_termbr
        data['h_elf_cm'] = h_elf_cm
        data['elf_cm3'] = elf_cm3

        data['fVol'] = self.fbrVlun
        data['fQco'] = self.fbrCO

        return data

    def getAlveolar(self):
        # Hartung2020_MATLAB/physiol_subst/physiology_alveolar.m
        data = {}
        data['fQco']=self.falvCO
        data['fVol']=self.fbrVlun
        data['ELF_cm3'] = self.Velf_alv
        data['Vol_cm3'] = self.OWlun * data['fVol']  # density = 1 [g/cm^3]
        data['Surf_cm2'] = self.Surf_alv;
        return data

class PKSubstanceLungParameters(EMObject):
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.fnSubst = String()

    def write(self, fnOut):
        fh=open(fnOut,"w")
        fh.write("%s # name\n"%self.name)
        fh.write("%f # maximum dissolution rate in alveolar space in units [nmol/(cm*min)]\n"%self.kdiss_alv)
        fh.write("%f # maximum dissolution rate in conducting airways in units [nmol/(cm*min)]\n"%self.kdiss_br)
        fh.write("%f # steady-state permeability in alveolar space in [cm/min]\n"%self.kp_alv)
        fh.write("%f # steady-state permeability in conducting airways in [cm/min]\n"%self.kp_br)
        fh.write("%f # solubility in alveolar space in [nmol/cm3]=[uM]\n" % self.Cs_alv)
        fh.write("%f # solubility in conducting airways in [nmol/cm3]=[uM]\n" % self.Cs_br)
        fh.write("%f # density in [nmol/cm3] = [uM]\n" % self.rho)
        fh.write("%f # molecular weight [g/mol]\n"%self.MW)
        fh.write("%f # plasma to lung partition coefficient in alveolar space [unitless]\n" % self.Kpl_alv)
        fh.write("%f # plasma to lung partition coefficient in conducting airways [unitless]\n" % self.Kpl_br)
        fh.write("%f # fraction unbound in plasma [unitless]\n" % self.fu)
        fh.write("%f # blood to plasma ratio [unitless]\n" % self.R)

        fh.close()
        self.fnSubst.set(fnOut)

    def read(self, fnIn):
        fh=open(fnIn)
        self.name=fh.readline().split()[0]
        self.kdiss_alv=float(fh.readline().split()[0])
        self.kdiss_br=float(fh.readline().split()[0])
        self.kp_alv=float(fh.readline().split()[0])
        self.kp_br=float(fh.readline().split()[0])
        self.Cs_alv=float(fh.readline().split()[0])
        self.Cs_br=float(fh.readline().split()[0])
        self.rho=float(fh.readline().split()[0])
        self.MW=float(fh.readline().split()[0])
        self.Kpl_alv=float(fh.readline().split()[0])
        self.Kpl_br=float(fh.readline().split()[0])
        self.fu=float(fh.readline().split()[0])
        self.R=float(fh.readline().split()[0])
        fh.close()
        self.fnSubst.set(fnIn)

    def getData(self):
        data = {}
        data['name'] = self.name
        data['kdiss_alv'] = self.kdiss_alv
        data['kdiss_br'] = self.kdiss_br
        data['kp_alv'] = self.kp_alv
        data['kp_br'] = self.kp_br
        data['Cs_alv'] = self.Cs_alv
        data['Cs_br'] = self.Cs_br
        data['rho'] = self.rho
        data['MW'] = self.MW
        data['Kpl_alv'] = self.Kpl_alv
        data['Kpl_br'] = self.Kpl_br
        data['fu'] = self.fu
        data['R'] = self.R
        return data

class PKDepositionParameters(EMObject):
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.fnSubstance = String()
        self.fnLung = String()
        self.fnDeposition = String()

    def setFiles(self, fnSubstance, fnLung, fnDeposition):
        self.fnSubstance.set(fnSubstance)
        self.fnLung.set(fnLung)
        self.fnDeposition.set(fnDeposition)

    def readDepositionFile(self, alvlim):
        fh=open(self.fnDeposition.get())
        state = 0

        for line in fh.readlines():
            if line.strip()=="":
                continue

            if state == 0: # Header
                tokens = line.split(':')
                if tokens[0]=='Total dose [ug]':
                    self.dose = float(tokens[1].strip())
                elif tokens[0]=='Diameter [um]':
                    self.diameterMode = tokens[1].strip()
                elif 'FractionDeposited' in tokens[0]:
                    diameters = []
                    oldGeneration = 0
                    fractionMatrix = []
                    fractionRow = []
                    state = 1
            elif state == 1: # table of numbers
                tokens = [float(x.strip()) for x in line.split()]
                diam = tokens[0]
                generation = tokens[1]
                fractionDeposited = tokens[2]
                if not diam in diameters:
                    diameters.append(diam)
                if oldGeneration!=generation:
                    if oldGeneration!=0:
                        fractionMatrix.append(fractionRow)
                    fractionRow=[]
                    oldGeneration=generation
                fractionRow.append(fractionDeposited)
        fractionMatrix.append(fractionRow)
        fh.close()

        diameters = np.asarray(diameters) # [um] for D
        if self.diameterMode=="aerodynamic":
            # Hartung2020_MATLAB/functions/aero2geom.m
            lambdaVar = 1; # spherical shape
            rho_water = 1; # [g / mL]
            rho_subst = self.substance.rho * self.substance.MW * 1e-9; # [nmol / mL] *[g / mol] ->[g / mL]
            diameters = diameters * np.sqrt( lambdaVar * rho_water / rho_subst);
            print("diameters",diameters)

        print("Bronchi",np.asarray(fractionMatrix[0:(alvlim-1)]))
        print("Alveolar",np.asarray(fractionMatrix[alvlim-1:]))
        self.dose_nmol = self.dose * 1e3 / self.substance.MW # % [ug] *1e3 / ([g/mol]) = [nmol]
        self.bronchiDose_nmol = np.asarray(fractionMatrix[0:(alvlim-1)])*self.dose_nmol
           # This is a matrix with as many rows as generations and as many columns as diameters
           # The content is what is the dose deposited at that generation in nmol
        self.alveolarDose_nmol = np.sum(np.asarray(fractionMatrix[alvlim-1:])*self.dose_nmol)
           # This is a matrix with as many rows as generations in the alveoli and as many columns as diameters
           # The content is what is the dose deposited at that generation in nmol

        self.throatDose = self.dose_nmol - np.sum(self.bronchiDose_nmol) - self.alveolarDose_nmol

        # Hartung2020_MATLAB/functions/diam2vol.m
        self.particleSize = math.pi/6 * np.power(1e-4 * diameters,3) # [cm^3]

    def read(self):
        self.substance = PKSubstanceLungParameters()
        self.substance.read(self.fnSubstance.get())

        self.lung = PKPhysiologyLungParameters()
        self.lung.read(self.fnLung.get())

        alveolarGeneration = len(self.lung.getBronchial()['type'])+1
        self.readDepositionFile(alveolarGeneration)

    def getData(self):
        data = {}
        data['bronchial'] = self.bronchiDose_nmol
        data['alveolar'] = self.alveolarDose_nmol
        data['throat'] = self.throatDose
        data['size'] = self.particleSize
        data['dose_nmol'] = self.dose_nmol
        return data
