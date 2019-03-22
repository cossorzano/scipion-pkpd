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
from pkpd.objects import PKPDVariable
from pyworkflow.protocol.constants import LEVEL_ADVANCED


class ProtPKPDScaleToCommonDose(ProtPKPD):
    """ Scale to common dose\n
        If the system is linear, then we may scale all measurements as if all individuals had been given the same dose.
        In this way we may construct a cleaner version of the response (by averaging) and use this cleaner version
        to find the initial parameters of the rest of samples.\n
        \n
        This protocol calculates the total dose of each individual and constructs new individuals such that
        Cnew = Cold * newDose/oldDose\n
        \n
        The old dose is evaluated in a period (by default, 1 week) that should include all doses
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'scale dose'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('measurementsToChange', params.StringParam, label="Measurements to change", default="Cp",
                      help='Example: Cp, Cu')
        form.addParam('t0', params.FloatParam, label="Initial time (h)", default=0, expertLevel=LEVEL_ADVANCED, help="Period in which the old dose will be evaluated")
        form.addParam('tF', params.FloatParam, label="Final time (h)", default=24*7, expertLevel=LEVEL_ADVANCED, help="Period in which the old dose will be evaluated")
        form.addParam('newDose', params.FloatParam, label="New total dose", default=10,
                      help='It must have the same units as the dose in the input experiment')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runChange',self.inputExperiment.get().getObjId(),self.measurementsToChange.get(),self.newDose.get(),
                                 self.t0.get(), self.tF.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runChange(self, objId, labels, newDose, t0, tF):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        self.printSection("Changing measurements")
        varNames = [token.strip() for token in self.measurementsToChange.get().split(',')]
        for sampleName, sample in self.experiment.samples.iteritems():
            sample.interpretDose()
            oldDose = sample.getCumulatedDose(self.t0.get(),self.tF.get())
            K = self.newDose.get()/oldDose

            for varName in varNames:
                values = sample.getValues(varName)
                for i in range(0,len(values)):
                    try:
                        value = values[i]
                        floatval = float(value)
                        values[i] = str(K*floatval)
                    except:
                        pass
                sample.setValues(varName,values)
        self.experiment.write(self._getPath("experiment.pkpd"))

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = "Measurements to change: %s"%self.measurementsToChange.get()
        msg.append("New dose: %f"%self.newDose.get())
        return msg

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    def _validate(self):
        msg = []
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD,False)
        tokens = self.measurementsToChange.get().split(',')
        for token in tokens:
            if not token.strip() in experiment.variables:
                msg.append("%s is not a variable of the experiment"%token)
            else:
                if experiment.variables[token].role != PKPDVariable.ROLE_MEASUREMENT:
                    msg.append("%s is not a measurement variable"%token)
        return msg

    def filterVarForWizard(self, v):
        """ Define the type of variables required (used in wizard). """
        return v.isMeasurement()