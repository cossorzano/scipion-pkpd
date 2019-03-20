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
# *  All comments concerning this program package may be sent to thes
# *  e-mail address 'info@kinestat.com'
# *
# **************************************************************************

import pyworkflow.protocol.params as params
from .protocol_pkpd_ode_base import ProtPKPDODEBase
from pkpd.models.pk_models import PK_TwocompartmentsClintMetabolite

# TESTED in test_workflow_gabrielsson_pk19.py

class ProtPKPDTwoCompartmentsClintMetabolite(ProtPKPDODEBase):
    """ Fit a two-compartmentx model to a set of measurements (any arbitrary dosing regimen is allowed)\n
        The central compartment is referred to as C, while the peripheral compartment as Cp.
        The differential equation is V dC/dt = -(Clint+Clp) * C + Clp * Cp + dD/dt, Vp dCp/dt = Clp * C - Clp * Cp, dCm/dt=Clint*C/Vm-Clm*Cm/Vm\n
        and Clint=Vmax/(Km+C)\n
        where C is the concentration of the central compartment, V, Vp, Vm the distribution volume of the central, peripheral and metabolite compartment,
        Clp is the distribution rate between the central and the peripheral compartments, \n
        Clm is the clearance of metabolite at the metabolite compartments, \n
        Vmax is the maximum processing capability of the metabolic pathway degrading the drug into the metabolite, \n
        Km is the Michaelis-Menten constant,and D the input dosing regime. \n
        Note that the intrinsic clearance occurs at the central volume.\n
        Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
        are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'pk two-compartments intrinsic, metabolite'

    def __init__(self,**kwargs):
        ProtPKPDODEBase.__init__(self,**kwargs)

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineParams1(form)
        form.addParam('predictor', params.StringParam, label="Time variable", default="t")
        form.addParam('predicted', params.StringParam, label="Plasma concentration", default="Cp")
        form.addParam('metabolite', params.StringParam, label="Metabolite concentration", default="Cm")
        form.addParam('bounds', params.StringParam, label="Parameter bounds ([tlag], sourceParameters, Vmax, Km, V, Clp, Vp, Clm, Vm)", default="",
                      help="Bounds for time delay, maximum processivity, Michaelis constant, volume and peripheral clearance and volume, clearance and volume of the metabolite. "\
                      'Make sure that the bounds are expressed in the expected units (estimated from the sample itself).'\
                      'Be careful that Cl bounds must be given here. If you have an estimate of the elimination rate, this is Ke=Cl/V. Consequently, Cl=Ke*V ')

    def getXYvars(self):
        self.varNameX=self.predictor.get()
        self.varNameY=[self.predicted.get(),self.metabolite.get()]

    def createModel(self):
        return PK_TwocompartmentsClintMetabolite()
