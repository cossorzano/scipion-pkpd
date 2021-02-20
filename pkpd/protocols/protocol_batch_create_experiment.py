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


from pyworkflow.protocol.params import PointerParam, StringParam
from pwem.protocols import BatchProtocol
from .protocol_pkpd import ProtPKPD
import copy


class BatchProtCreateExperiment(BatchProtocol, ProtPKPD):
    """ Create experiment.\n
        Protocol created by http://www.kinestatpharma.com
    """
    _label = 'create experiment'

    def _defineParams(self, form):
        form.addHidden('inputExperiment', PointerParam, pointerClass='EMObject')
        form.addHidden('listOfSamples', StringParam)
        form.addHidden('newTitle', StringParam)
        form.addHidden('newComment', StringParam)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------

    def createOutputStep(self):
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        newExperiment = experiment.subset(self.listOfSamples.get().split(';'))
        newExperiment.general["title"]=self.newTitle.get().strip()
        newExperiment.general["comment"]=self.newComment.get().strip()

        self.writeExperiment(newExperiment,self._getPath("experiment.pkpd"))
        self._defineOutputs(outputExperiment=newExperiment)
        self._defineSourceRelation(self.inputExperiment, newExperiment)
