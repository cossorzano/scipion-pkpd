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

from pyworkflow.protocol.params import PointerParam, StringParam, TextParam
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDExperiment

class ProtPKPDCreateExperiment(ProtPKPD):
    """ Create experiment.\n
        Protocol created by http://www.kinestatpharma.com
    """
    _label = 'create experiment'

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('newTitle', StringParam, label="Title", default="New experiment")
        form.addParam('newComment', StringParam, label="Comment", default="")
        form.addParam('newVariables', TextParam, label="Variables", default="",
                      help="varName ; units ; numeric[%f]/text[%s] ; label/time/measurement ; comment")
        form.addParam('newVias', TextParam, label="Vias", default="",
                      help="viaName; iv/ev0/ev1/ev01/ev0tlag1/evFractional/ev1-ev1/doubleWeibull/splineN/splineXYN; tlag[=0.000000] min; bioavailability[=1.000000]\n"
                           "It can be empty")
        form.addParam('newDoses', TextParam, label="Doses", default="",
                      help="doseName; via=viaName; bolus/repeated_bolus/infusion; t=0.000000 min; d=10 mg\n"
                           "It can be empty")
        form.addParam('newGroups', TextParam, label="Groups", default="",
                      help="groupName\n"
                           "It can be empty")
        form.addParam('newSamples', TextParam, label="Samples", default="",
                      help="sampleName; [dose=doseName]; [group=groupName]; [varName=value]")
        form.addParam('newMeasurements', TextParam, label="Measurements", default="",
                      help="sampleName1; t; varMeasurement\n"
                           "t0; v0\n"
                           "t1; v1\n"
                           "\n"
                           "sampleName2; t; varMeasurement\n"
                           "t0; v0\n"
                           "t1; v1\n")

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------

    def createOutputStep(self):
        fnTmp = self._getExtraPath("aux.pkpd")
        fh = open(fnTmp,"w")

        fh.write("[EXPERIMENT] ===========================\n")
        fh.write("title=%s\n"%self.newTitle.get())
        fh.write("comment=%s\n"%self.newComment.get())
        fh.write("\n")

        fh.write("[VARIABLES] ============================\n")
        fh.write("%s\n"%self.newVariables.get())
        fh.write("\n")

        fh.write("[VIAS] ================================\n")
        fh.write("%s\n"%self.newVias.get())
        fh.write("\n")

        fh.write("[DOSES] ================================\n")
        fh.write("%s\n"%self.newDoses.get())
        fh.write("\n")

        fh.write("[GROUPS] ================================\n")
        fh.write("%s\n"%self.newGroups.get())
        fh.write("\n")

        fh.write("[SAMPLES] ================================\n")
        fh.write("%s\n"%self.newSamples.get())
        fh.write("\n")

        fh.write("[MEASUREMENTS] ===========================\n")
        fh.write("%s\n"%self.newMeasurements.get())
        fh.write("\n")
        fh.close()

        experiment = PKPDExperiment()
        experiment.load(fnTmp, verifyIntegrity=False)

        self.writeExperiment(experiment,self._getPath("experiment.pkpd"))
        self._defineOutputs(outputExperiment=experiment)
