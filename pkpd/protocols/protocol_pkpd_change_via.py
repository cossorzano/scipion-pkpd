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
from pkpd.pkpd_units import PKPDUnit

# TESTED in test_workflow_gabrielsson_pk03.py, test_workflow_deconvolution.py

class ProtPKPDChangeVia(ProtPKPD):
    """ Change via of administration\n
        This protocol may also be used to change the bioavailability or the tlag
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'change via'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('viaName', params.StringParam, label="Via name",
                      help='Name of the via you want to change')
        form.addParam('newViaType', params.StringParam, label="New via type",
                      help='New via type, leave it empty to keep current via. Valid vias are iv, ev0, ev1, ev01, evFractional, '
                           'ev0tlag1, spline2, spline3, ..., spline20, splinexy2, splinexy3, ..., splinexy20')
        form.addParam('tlag', params.StringParam, label="New tlag",
                      help='New tlag of the dose, leave it empty to let it free so that it can be optimized by an ODE model')
        form.addParam('bioavailability', params.StringParam, label="New bioavailability",
                      help='New bioavailability of the dose, leave it empty to let it free so that it can be optimized by an ODE model')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runChange',self.inputExperiment.get().getObjId(),self.viaName.get(), self.newViaType.get(),
                                 self.tlag.get(),self.bioavailability.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runChange(self, objId, viaName, newViaType, tlag, bioavailability):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        via = self.experiment.vias[viaName]
        if newViaType!="":
            via.via = newViaType

        if not tlag:
            if not 'tlag' in via.paramsToOptimize:
                via.paramsToOptimize.append("tlag")
                via.paramsUnitsToOptimize.append(PKPDUnit.UNIT_TIME_MIN)
        else:
            if 'tlag' in via.paramsToOptimize:
                idx=via.paramsToOptimize.index('tlag')
                del via.paramsToOptimize[idx]
                del via.paramsUnitsToOptimize[idx]
            via.tlag=float(tlag)

        if not bioavailability:
            if not 'bioavailability' in via.paramsToOptimize:
                via.paramsToOptimize.append("bioavailability")
                via.paramsUnitsToOptimize.append(PKPDUnit.UNIT_NONE)
        else:
            if 'bioavailability' in via.paramsToOptimize:
                idx=via.paramsToOptimize.index('bioavailability')
                del via.paramsToOptimize[idx]
                del via.paramsUnitsToOptimize[idx]
            via.bioavailability=float(bioavailability)
        self.experiment.write(self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        if not self.viaName.get() in experiment.vias:
            errors.append("%s is not a via of the experiment"%self.viaName.get())
        return errors

    def _summary(self):
        msg=[]
        return msg
