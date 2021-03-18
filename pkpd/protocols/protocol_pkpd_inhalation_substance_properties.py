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

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD

from pkpd.objects import PKSubstanceLungParameters

# Tested in test_workflow_inhalation1

class ProtPKPDInhSubstanceProperties(ProtPKPD):
    """ Produce a description of the properties of a substance related to inhalation\n
        See Hartung2020.
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'substance parameters'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        # Default parameters from
        # Hartung2020_MATLAB/physiol_subst/properties_subst.m

        form.addSection('Input')
        form.addParam('name', params.StringParam, label="Name")

        groupS = form.addGroup('Substance parameters')
        groupS.addParam('rho', params.FloatParam, label="Density [g/cm3]", default=1.3)
        groupS.addParam('MW', params.FloatParam, label="Molecular weight [g/mol]", default=430.53)

        groupA = form.addGroup('Parameters of the substance in the lung')
        groupA.addParam('kdiss_alv', params.FloatParam,
                        label="Maximum dissolution rate in alveolar space [nmol/(cm*min)]", default=3.3e-4)
        groupA.addParam('kdiss_br', params.FloatParam,
                        label="Maximum dissolution rate in alveolar space [nmol/(cm*min)]", default=3.3e-4)
        groupA.addParam('kp_alv', params.FloatParam,
                        label="Steady-state permeability in alveolar space [cm/min]", default=5.33e-6*60)
        groupA.addParam('kp_br', params.FloatParam,
                        label="Steady-state permeability in conducting airways [cm/min]", default=5.33e-6*60)
        groupA.addParam('Cs_alv', params.FloatParam,
                        label="Solubility in alveolar space [nmol/cm3]=[uM]", default=69.797)
        groupA.addParam('Cs_br', params.FloatParam,
                        label="Solubility in conducting airways [nmol/cm3]=[uM]", default=69.797)

        groupB = form.addGroup('Lung-blood partition parameters')
        groupB.addParam('Kpl_alv', params.FloatParam,
                        label="Plasma to lung partition coefficient in alveolar space [unitless]", default=8)
        groupB.addParam('Kpl_br', params.FloatParam,
                        label="Plasma to lung partition coefficient in conducting airways [unitless]", default=8)

        groupC = form.addGroup('Substance parameters in blood')
        groupC.addParam('fu', params.FloatParam,
                        label="Fraction unbound in plasma [unitless]", default=0.161)
        groupC.addParam('R', params.FloatParam,
                        label="Blood to plasma ratio [unitless]", default=0.8)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runParams')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runParams(self):
        self.substanceParams = PKSubstanceLungParameters()
        self.substanceParams.name = self.name.get()
        self.substanceParams.rho=self.rho.get()/self.MW.get()*1e9
        self.substanceParams.MW=self.MW.get()
        self.substanceParams.kdiss_alv=self.kdiss_alv.get()
        self.substanceParams.kdiss_br=self.kdiss_br.get()
        self.substanceParams.kp_alv=self.kp_alv.get()
        self.substanceParams.kp_br=self.kp_br.get()
        self.substanceParams.Cs_alv=self.Cs_alv.get()
        self.substanceParams.Cs_br=self.Cs_br.get()
        self.substanceParams.Kpl_alv=self.Kpl_alv.get()

        self.substanceParams.Kpl_br=self.Kpl_br.get()
        self.substanceParams.fu=self.fu.get()
        self.substanceParams.R=self.R.get()

        self.substanceParams.write(self._getPath("substance_parameters.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputSubstanceParameters=self.substanceParams)

    def _citations(self):
        return ['Hartung2020']