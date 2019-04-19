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
from pkpd.models.pk_models import PK_TwocompartmentsUrine

# TESTED in test_workflow_gabrielsson_pk16.py
# TESTED in test_workflow_gabrielsson_pk25.py


class ProtPKPDTwoCompartmentsUrine(ProtPKPDODEBase):
    """ Fit a two-compartments model to a set of plasma and urine (cumulated) measurements ((any arbitrary dosing regimen is allowed)\n
        The differential equation is dC/dt = -Cl * C/V -Clp *(C-Cp)/V + 1/V * dD/dt, dCp/dt=Cl*C/Vp+Clp*(C-Cp)/Vp and dA/dt = fe * Cl * C\n
        where C is the concentration, Cl the total clearance (metabolic and excretion), V the distribution volume, Clp is the Clearance \n
        to the peripheric compartment, Vp is the volume of the peripheric compartment, D the input dosing regime and fe the fraction excreted.
Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'two-compartments urine'

    def __init__(self,**kwargs):
        ProtPKPDODEBase.__init__(self,**kwargs)

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineParams1(form)
        form.addParam('predictor', params.StringParam, label="Time variable", default="t")
        form.addParam('predicted', params.StringParam, label="Plasma concentration", default="Cp")
        form.addParam('Au', params.StringParam, label="Cumulated amount in urine", default="Au")
        form.addParam('bounds', params.StringParam, label="Parameter bounds ([tlag], Cl, V, Clp, Vp, fe)", default="",
                      help="Bounds for the tlag (if it must be estimated), clearance, volume and fraction excreted."\
                      'Make sure that the bounds are expressed in the expected units (estimated from the sample itself).'\
                      'If tlag must be estimated, its bounds must always be specified')

    def getXYvars(self):
        self.varNameX=self.predictor.get()
        self.varNameY=[self.predicted.get(),self.Au.get()]

    def createModel(self):
        return PK_TwocompartmentsUrine()
