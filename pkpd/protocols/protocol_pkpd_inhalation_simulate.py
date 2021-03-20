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

from pkpd.objects import PKDepositionParameters, PKCiliarySpeed, PKInhalationDissolution

# Tested in test_workflow_inhalation1

class ProtPKPDInhSimulate(ProtPKPD):
    """ Simulate inhalation PK\n
        See Hartung2020.
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'simulate inhalation'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('ptrDeposition', params.PointerParam, pointerClass='PKDepositionParameters', label="Deposition")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulation')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runSimulation(self):
        self.deposition = PKDepositionParameters()
        self.deposition.setFiles(self.ptrDeposition.get().fnSubstance.get(),
                                 self.ptrDeposition.get().fnLung.get(),
                                 self.ptrDeposition.get().fnDeposition.get())

        self.ciliarySpeed = PKCiliarySpeed()
        self.ciliarySpeed.setFiles(self.ptrDeposition.get().fnLung.get())
        self.ciliarySpeed.prepare()

        self.inhalationDissolution = PKInhalationDissolution()
        self.inhalationDissolution.setFiles(self.ptrDeposition.get().fnSubstance.get())
        self.inhalationDissolution.prepare()

    def createOutputStep(self):
        pass

    def _citations(self):
        return ['Hartung2020']