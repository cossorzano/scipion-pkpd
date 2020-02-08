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


from pyworkflow.protocol.params import PointerParam, MultiPointerParam, IntParam, EnumParam
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDExperiment


class ProtPKPDSplitGather(ProtPKPD):
    """ Split an experiment into small pieces, or gather these small pieces into an experiment.\n
        Protocol created by http://www.kinestatpharma.com
    """
    _label = 'split/gather experiment'

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('splitOrGather', EnumParam, label="Operation", choices=['Split','Gather'], default=0)
        form.addParam('inputExperiment', PointerParam, label="Experiment to split", pointerClass='PKPDExperiment', condition='splitOrGather==0')
        form.addParam('groupSize', IntParam, label="Group size", default=5, condition='splitOrGather==0')
        form.addParam('inputExperiments', MultiPointerParam, label="Experiments to gather", pointerClass='PKPDExperiment', condition='splitOrGather==1')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def createSubGroup(self, g, experiment, listOfSamples):
        newExperiment = experiment.subset(listOfSamples)
        self.writeExperiment(newExperiment, self._getPath("experiment%05d.pkpd"%g))
        self._defineOutputs(**{"outputExperiment%d"%g: newExperiment})
        self._defineSourceRelation(self.inputExperiment, newExperiment)

    def createOutputStep(self):
        if self.splitOrGather.get()==0:
            experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
            listOfSamples=[]
            g=1
            for sampleName in experiment.samples:
                listOfSamples.append(sampleName)
                if len(listOfSamples)==self.groupSize.get():
                    self.createSubGroup(g, experiment, listOfSamples)
                    listOfSamples=[]
                    g+=1
            if len(listOfSamples)!=0:
                self.createSubGroup(g, experiment, listOfSamples)
        else:
            newExperiment = PKPDExperiment()
            for ptrExperiment in self.inputExperiments:
                newExperiment.gather(self.readExperiment(ptrExperiment.get().fnPKPD))
            self.writeExperiment(newExperiment, self._getPath("experiment.pkpd"))
            self._defineOutputs(outputExperiment=newExperiment)
            for ptrExperiment in self.inputExperiments:
                self._defineSourceRelation(ptrExperiment, newExperiment)
