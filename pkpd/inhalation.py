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
from scipy.interpolate import interp1d, interp2d
from pwem.objects import EMObject
from pyworkflow.object import String, Integer

def diam2vol(diameters):
    # Hartung2020_MATLAB/functions/diam2vol.m
    # Convert diameter to volume, assuming spherical shape
    # Units: [um] for diameters, [cm^3] for V
    return math.pi / 6 * np.power(1e-4 * diameters, 3)

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
            # ELF = Epithelial Lining Fluid
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
        data['voltis_cm3'] = voltis_cm3

        data['fVol'] = self.fbrVlun
        data['fQco'] = self.fbrCO

        return data

    def getAlveolar(self):
        # Hartung2020_MATLAB/physiol_subst/physiology_alveolar.m
        data = {}
        data['fQco']=self.falvCO
        data['fVol']=1-self.fbrVlun
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
        self.particleSize = diam2vol(diameters) # cm3

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

class PKCiliarySpeed(EMObject):
    # Hartung2020_MATLAB/functions/get_cilspeed.m
    # Only 'interp' model is implemented here

    EXPON = 0
    FIT = 1
    INTERP = 2

    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.type = Integer()
        self.type.set(self.INTERP)

        self.lungParams = None

    def prepare(self, lungParams):
        self.lungParams = lungParams
        lungData = lungParams.getBronchial()

        self.Lx = np.sum(lungData['length_cm'])
        self.pos = lungData['pos']
        self.diam_cm = lungData['diam_cm']

        print("Ciliary speed ================")
        print("Total length of the lung [cm]",self.Lx)
        print("Location of branches [cm]",self.pos)
        print("Diameter of the branches [cm]",self.diam_cm)

        if self.type.get()==self.INTERP:
            # Hofmann - Sturm model for mucociliary transport velocity
            # Reference: Hofmann / Sturm (2004) J Aerosol Med, Vol 17(1)
            self.v = 0.1 * 1.2553 * np.power(lungData['diam_cm'],2.808) # in [cm / min]
            self.cilspeed = interp1d(self.pos, self.v, bounds_error=False, fill_value=(self.v[0], self.v[-1]))

class PKInhalationDissolution(EMObject):
    # Hartung2020_MATLAB/functions/get_disol.m
    # All implemented dissolution models are based on an adapted version of
    # the Noyes - Whitney equation, but differ in whether they
    # 1) describe saturable or unsaturable dissolution
    # 2) truncate or not the dissolution speed for small particles (for  numeric stability)
    # 3) truncate or not the dissolution speed for large particles
    #    (based on the assumption that only a part of the particle surface
    #    is in contact with the dissolution medium)

    UNSAT       = 0
    TRUNC_UNSAT = 1
    SAT         = 2 # Only this one is implemented
    TRUNC_SAT   = 3
    TRUNC2_SAT  = 4
    CAP         = 5
    XCAP        = 6
    UNSOL       = 7

    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.type = Integer()
        self.type.set(self.SAT)
        self.substanceParams = None

    def prepare(self, substanceParams, part):
        self.substanceParams = substanceParams
        substanceData = substanceParams.getData()

        if part=="bronchi":
            self.kdiss = substanceData['kdiss_br']
            self.Cs = substanceData['Cs_br']
        else:
            self.kdiss = substanceData['kdiss_alv']
            self.Cs = substanceData['Cs_alv']
        self.rho = substanceData['rho']

        print("Substance in %s ==================="%part)
        print("Maximum dissolution rate [nmol/(cm*min)]",self.kdiss)
        print("Solubility in [nmol/cm3]=[uM]",self.Cs)
        print("Density in [nmol/cm3]",self.rho)

        D = self.kdiss / self.Cs
        K = 4 * math.pi * D / (self.rho * np.power(4.0/3.0 * math.pi, 1.0/3.0))

        if self.type.get()==self.SAT:
            self.dissol = lambda s, Cf: np.where(s>0, K*(self.Cs-Cf) * np.power(s,1.0/3.0), 0.0) # [cm^3/min]

class PKLung(EMObject):
    # Hartung2020_MATLAB/functions/get_bronchial_kinetics.m
    # Hartung2020_MATLAB/functions/get_alveolar_kinetics.m
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.substanceParams = None
        self.lungParams = None
        self.ciliarySpeed = None
        self.inhalationDissolutionBronchi = None
        self.inhalationDissolutionAlveoli = None

    def prepare(self, substanceParams, lungParams):
        self.substanceParams = substanceParams
        self.lungParams = lungParams

        self.ciliarySpeed = PKCiliarySpeed()
        self.ciliarySpeed.prepare(lungParams)

        self.inhalationDissolutionBronchi = PKInhalationDissolution()
        self.inhalationDissolutionBronchi.prepare(substanceParams,"bronchi")

        # Hartung2020_MATLAB/functions/get_bronchial_kinetics.m
        bronchialData = lungParams.getBronchial()
        pos = bronchialData['pos']
        diam_cm = bronchialData['diam_cm']

        self.substanceData = substanceParams.getData()
        kp = self.substanceData['kp_br']
        ka = math.pi * kp * diam_cm
        print("Absorption rate per bronchial segment [1/min]",ka)
        self.bronchialAbsorb = interp1d(pos, ka, bounds_error=False, fill_value=(ka[0], ka[-1]))

        # Hartung2020_MATLAB/functions/get_alveolar_kinetics.m
        self.inhalationDissolutionAlveoli = PKInhalationDissolution()
        self.inhalationDissolutionAlveoli.prepare(substanceParams,"alveoli")

    def cilspeed(self, x):
        return self.ciliarySpeed.cilspeed(x)

def conserving_projection(X1, V1, X2):
    # Hartung2020_MATLAB/functions/conserving_projection.m
    #   V2 = CONSERVING_PROJECTION(X1,V1,X2) projects quantity V1 from location
    #   grid X1 to grid X2, ensuring int_x^y(v2(z)dz) = int_x^y(v1(z)dz) holds
    #   for any x,y, where vI(z) := vI(j) for Xi(j) < z < Xi(j+1) (i=1,2)
    #   (constant between grid cells)
    cumV1 = np.concatenate(([0],np.cumsum(V1)));
    interpolator = interp1d(X1, cumV1, bounds_error=False, fill_value=(cumV1[0], cumV1[-1]))
    return np.diff(interpolator(X2))

def P_aELF(lungData, Xbnd):
    # Hartung2020_MATLAB/functions/P_aELF.m
    #   Aim: locally and globally conserve ELF volumes
    #   Input:  boundary gridpoints of computational location grid
    #   Output: cross-sectional areas at (ctr) grid points (a_Xctr), satisfying
    #           int_dx(Xbnd, a_Xctr) == sum(aw.elf_cm3) (and local volume conservation)
    V_Xctr = conserving_projection(np.concatenate(([0],lungData['end_cm'])), lungData['elf_cm3'], Xbnd)
    return np.divide(V_Xctr, np.diff(Xbnd))

def P_aTis(lungData, Xbnd):
    # Hartung2020_MATLAB/functions/P_aTis.m
    #   Aim: locally and globally conserve tissue volumes
    #   Input:  boundary gridpoints of computational location grid
    #   Output: cross-sectional areas at (ctr) grid points (a_Xctr), satisfying
    #           int_dx(Xbnd, a_Xctr) == sum(aw.voltis_cm3) (and local volume conservation)
    V_Xctr = conserving_projection(np.concatenate(([0],lungData['end_cm'])), lungData['voltis_cm3'], Xbnd)
    return np.divide(V_Xctr, np.diff(Xbnd))

def P_hELF(lungData, Xctr):
    # Hartung2020_MATLAB/functions/P_hELF.m
    #   Input:  center gridpoints of computational location grid
    #   Output: ELF heights at (ctr) grid points (h_Xctr), in cm.
    interpolator = interp1d(lungData['pos'], lungData['h_elf_cm'], fill_value='extrapolate')
    return interpolator(Xctr)

def P_Qbr(lungData, systemicData, Xbnd):
    # Hartung2020_MATLAB/functions/P_Qbr.m
    #   Project bronchial (central lung) blood flow onto computational grid
    #   Total blood flow is conserved
    Qcl = lungData['fQco'] * systemicData['Qco']

    # Modelling assumption: blood flow proportional to tissue volume
    Qgen = Qcl * lungData['voltis_cm3'] / np.sum(lungData['voltis_cm3'])
    Q_Xctr = conserving_projection(np.concatenate(([0],lungData['end_cm'])), Qgen, Xbnd)
    return np.divide(Q_Xctr,np.diff(Xbnd));

def project_deposition_2D(depositionData,X,S,lungData):
    #   Project deposition data on 2D computational grid
    #   Strategy for projection on 2D grid:
    #     - Uniform distribution of dose in data location-size gridcells
    #       (per generation and per size category) --> f(x,s)
    #     - rho0 at solver location-size gridcell C = average of f(x,s) in C
    xbnddat = np.concatenate(([0],lungData['end_cm']))

    smax = diam2vol(50)
    sdil_default = diam2vol(0.1)

    # Projection onto bronchi
    # step 1: distribute delta peaks in size over short size intervals
    sdattmp = depositionData['size']
    sdil = np.min([sdil_default, 0.5*np.min(sdattmp)])
    if sdattmp.size>1:
        sdil=np.min([sdil,np.min(np.diff(sdil))])
    sbnddat = np.concatenate(([0],np.kron(sdattmp,[1,1]) + sdil*np.kron(np.ones(sdattmp.shape),[-1,1]), [smax]))

    amtxs_br = depositionData['bronchial']
    amtxs_br_ext = np.zeros((amtxs_br.shape[0],2*amtxs_br.shape[1]+1))
    amtxs_br_ext[:,1::2] = amtxs_br

    # step 2: compute cumulative amount matrix per data location-size grid
    camtdat_br_x = np.cumsum(np.pad(amtxs_br_ext,((1,0),(0,0)),'constant',constant_values=0),axis=0)
    camtdat_br_xs = np.cumsum(np.pad(camtdat_br_x,((0,0),(1,0)),'constant',constant_values=0),axis=1)

    # step 3: linear interpolation projects onto solver location-size grid
    interpolator = interp2d(sbnddat, xbnddat, camtdat_br_xs, bounds_error=False, fill_value=0)
    camtbnd_br_xs = interpolator(S, X)

    amtgrd_br_x = np.diff(camtbnd_br_xs, axis=0)
    amtgrd_br_xs = np.diff(amtgrd_br_x, axis=1)

    amtgrd_br_xs = np.where(amtgrd_br_xs<0,0,amtgrd_br_xs) # Make sure that everything is positive

    # step 4: for each grid cell: amount -> location-resolved amount per size
    dx = np.diff(X)
    ds = np.diff(S)
    s = 0.5*(S[0:-1] + S[1:])
    rho0br = np.divide(amtgrd_br_xs,np.reshape(dx,(dx.size,1))*np.reshape(np.multiply(s,ds),(1,s.size)))

def saturable_2D_upwind_IE(lungParams, pkLung, depositionParams, tt, Sbnd):
    # Hartung2020_MATLAB/models/saturable_2D_upwind_IE.m
    #   Algorithm features:
    #   - Conducting airways, peripheral airways and systemic circulation fully coupled
    #   - Upwind discretisation of PDE and compatible integration rules for mass conservation
    #   - implicit Euler discretisation of linear processes
    #   - saturable dissolution model
    #   - initial particle size treated as a delta distribution; particle size is treat followed over time.
    #   - along the location axis, an arbitrary grid can be used
    #
    #   Since dissolution is saturable, particle sizes at different airway
    #   positions may differ from each other for t > 0.
    #
    #   The idea of this approach is that the size resolution in the data may
    #   be much coarser than that the location grid
    #   (extreme case: monodisperse particle distribution)
    dt = np.diff(tt)

    lungData = lungParams.getBronchial()
    alveolarData = lungParams.getAlveolar()
    systemicData = lungParams.getSystemic()
    depositionData = depositionParams.getData()

    Xbnd = np.sort([0] + lungData['end_cm'].tolist() + lungData['pos'].tolist())
    Nx = Xbnd.size - 1
    Ns = Sbnd.size - 1

    # midpoints in discretisation cell (where rho, Cflu, Ctis are modelled)
    Xctr = Xbnd[:-1] + np.diff(Xbnd)/2
    Sctr = Sbnd[:-1] + np.diff(Sbnd)/2

    # transport velocity
    lambdaX = pkLung.cilspeed(Xbnd)

    # location/size discretisation steps
    dx = np.diff(Xbnd)
    ds = np.diff(Sbnd)

    # absorption into tissue
    kaX = pkLung.bronchialAbsorb(Xctr)

    # cross-sectional areas
    aFluX = P_aELF(lungData, Xbnd)
    aTisX = P_aTis(lungData, Xbnd)

    # ELF heights in bronchi / alveolar space
    hFluX = P_hELF(lungData, Xctr)
    hFlualv = alveolarData['ELF_cm3']/alveolarData['Surf_cm2']

    # blood flow
    Qalv = alveolarData['fQco'] * systemicData['Qco'] # Alveolar
    qX = P_Qbr(lungData, systemicData, Xbnd)
    QbrX = np.multiply(qX, dx) # bronchial (location-resolved)

    # plasma to lung partition coefficients
    Kpl_br = pkLung.substanceData['Kpl_br']
    Kpl_alv = pkLung.substanceData['Kpl_alv']
    Kpl_u_br = Kpl_br / pkLung.substanceData['fu']
    Kpl_u_alv = Kpl_alv / pkLung.substanceData['fu']

    # permeability-surface area products [cm3/min]
    PS_alv = pkLung.substanceData['kp_alv'] * alveolarData['Surf_cm2'];  # alveolar
    PS_br = np.multiply(kaX , dx) # bronchial

    # assign volumes to variables
    alvELF = alveolarData['ELF_cm3']  # alveolar fluid
    alvTis = alveolarData['Vol_cm3']  # alveolar tissue
    brELF  = np.multiply(aFluX,dx)    # bronchial fluid
    brTis  = np.multiply(aTisX,dx)    # bronchial tissue

    # Allocate lung amounts: A_flu(alv/br), A_tis(alv/br)
    Nt=tt.size-1
    Aalvflu = np.zeros((Nt+1,1));
    Aalvtis = np.zeros((Nt+1,1));

    Abrflu  = np.zeros((Nt+1,Nx));
    Abrtis  = np.zeros((Nt+1,Nx));

    # Allocate systemic amounts: Asysgut, Asysctr, Asysper
    Asysgut   = np.zeros((Nt+1,1));
    Asysctr   = np.zeros((Nt+1,1));
    Asysper   = np.zeros((Nt+1,1));

    # Allocate amount cleared:
    Aclear = np.zeros((Nt+1,1));     # for mass balance
    Amcc   = np.zeros((Nt+1,1));     # to track cumulative amount cleared by MCC (not in mass balance)

    # Initialize PSPM densities (br/alv) and gut compartment with dosing
    rhobr  = np.zeros((Nt+1,Nx,Ns));
    rhoalv = np.zeros((Nt+1, 1,Ns));

    project_deposition_2D(depositionData, Xbnd, Sbnd, lungData)