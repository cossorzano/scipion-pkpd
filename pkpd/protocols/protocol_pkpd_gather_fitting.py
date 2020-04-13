# **************************************************************************
# *
# * Authors:  Carlos Oscar Sorzano (info@kinestat.com)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************


from pyworkflow.protocol.params import PointerParam, MultiPointerParam
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDExperiment, PKPDFitting


class ProtPKPDGatherFitting(ProtPKPD):
    """ Gather several fittings coming from a split experiment into a single fitting.\n
        Protocol created by http://www.kinestatpharma.com
    """
    _label = 'gather fitting'

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputFittings', MultiPointerParam, label="Fittings to gather", pointerClass='PKPDFitting')
        form.addParam('inputExperiment', PointerParam, label="Experiment this fitting is linked to", pointerClass='PKPDExperiment')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        newFitting = PKPDFitting()
        newFitting.fnExperiment.set(self.inputExperiment.get().fnPKPD)
        for ptrFitting in self.inputFittings:
            newFitting.gather(self.readFitting(ptrFitting.get().fnFitting))
        self.writeExperiment(newFitting, self._getPath("fitting.pkpd"))
        self._defineOutputs(outputFitting=newFitting)
        for ptrFitting in self.inputFittings:
            self._defineSourceRelation(ptrFitting, newFitting)
