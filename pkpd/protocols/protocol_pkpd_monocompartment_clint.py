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
from pkpd.models.pk_models import PK_MonocompartmentClint

# TESTED in test_workflow_gabrielsson_pk17.py
# TESTED in test_workflow_gabrielsson_pk20.py

class ProtPKPDMonoCompartmentClint(ProtPKPDODEBase):
    """ Fit a monocompartmental model to a set of measurements obtained by oral doses (any arbitrary dosing regimen is allowed)\n
        The differential equation is dC/dt = -Clint * C/V + 1/V * dD/dt and Clint=Vmax/(Km+C)\n
        where C is the concentration, Clint the intrinsic clearance, V the distribution volume, Vmax is the maximum\n
        processing capability of the metabolic pathway degrading the drug, Km is the Michaelis-Menten constant, \n
        and D the input dosing regime. Confidence intervals calculated by this fitting may be pessimistic because\n
        it assumes that all model parameters are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'pk monocompartment intrinsic'

    def __init__(self,**kwargs):
        ProtPKPDODEBase.__init__(self,**kwargs)

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineParams1(form, True, "t", "Cp")
        form.addParam('bounds', params.StringParam, label="Parameter bounds ([tlag], sourceParameters, Vmax, Km, V)", default="",
                      help="Bounds for the tlag (if it must be estimated), parameters for the source, maximum processivity, Michaelis constant and volume. Example: (0.01,0.04);(0,10);(0.2,0.4);(10,20). "\
                      'Make sure that the bounds are expressed in the expected units (estimated from the sample itself).')

    def createModel(self):
        return PK_MonocompartmentClint()
