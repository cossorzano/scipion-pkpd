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
from pwem.objects import *

class PKPhysiologyLungParameters(EMObject):
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.fnPhys = String()

    def write(self, fnOut):
        fh=open(fnOut,"w")
        fh.write("Qco=%f   [mL/min]   # Cardiac output\n"%self.Qco)
        fh.write("OWlun=%f [g]        # lung tissue weight, without blood\n"%self.OWlun)
        fh.write("falvCO=%f           # alveolar fraction of cardiac output\n"%self.falvCO)
        fh.write("Velf_alv=%f [mL]    # alveolar ELF volume\n"%self.Velf_alv)
        fh.write("Surf_alv=%f [cm2]   # alveolar surface area\n"%self.Surf_alv)
        fh.write("fbrCO=%f            # bronchial fraction of cardiac output\n"%self.fbrCO)
        fh.write("fbrVlun=%f          # bronchial fraction of lung tissue volume\n"%self.fbrVlun)
        fh.write("helf_trach=%f [cm]  # bronchial ELF heights for interpolation, trachea\n" % self.helf_trach)
        fh.write("helf_termbr=%f [cm] # bronchial ELF heights for interpolation, terminal bronchioles\n" % self.helf_termbr)
        fh.close()
        self.fnPhys.set(fnOut)

    def load(self, fnIn):
        def readParam(fh):
            return float(fh.readline().split()[0].split('=')[1])
        fh.open(fnIn)
        self.Qco=readParam(fh)
        self.OWlun=readParam(fh)
        self.falvCO=readParam(fh)
        self.Velf_alv=readParam(fh)
        self.Surf_alv=readParam(fh)
        self.fbrCO=readParam(fh)
        self.fbrVlun=readParam(fh)
        self.helf_trach=readParam(fh)
        self.helf_termbr=readParam(fh)
        fh.close()
        self.fnPhys.set(fnIn)
