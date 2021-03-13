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
from pwem.objects import *

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

    def load(self, fnIn):
        fh.open(fnIn)
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
        self.fnPhys = String()

    def write(self, fnOut):
        fh=open(fnOut,"w")
        fh.write("%s # name\n"%self.name)
        fh.write("%f # maximum dissolution rate in alveolar space in units [nmol/(cm*min)]\n"%self.kdiss_alv)
        fh.write("%f # maximum dissolution rate in conducting airways in units [nmol/(cm*min)]\n"%self.kdiss_br)
        fh.write("%f # steady-state permeability in alveolar space in [cm/min]\n"%self.kp_alv)
        fh.write("%f # steady-state permeability in conducting airways in [cm/min]\n"%self.kp_br)
        fh.write("%f # solubility in alveolar space in [nmol/cm3]=[uM]\n" % self.Cs_alv)
        fh.write("%f # solubility in conducting airways in [nmol/cm3]=[uM]\n" % self.Cs_br)
        fh.write("%s # density in [nmol/cm3] = [uM]\n" % self.rho)
        fh.write("%f # molecular weight [g/mol]\n"%self.MW)
        fh.write("%s # plasma to lung partition coefficient in alveolar space [unitless]\n" % self.Kpl_alv)
        fh.write("%s # plasma to lung partition coefficient in conducting airways [unitless]\n" % self.Kpl_br)
        fh.write("%s # fraction unbound in plasma [unitless]\n" % self.fu)
        fh.write("%s # blood to plasma ratio [unitless]\n" % self.R)

        fh.close()
        self.fnPhys.set(fnOut)

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