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

import copy
import numpy as np

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDVariable

class ProtPKPDOperateExperiment(ProtPKPD):
    """ Create a new label or measurement to an existing experiment.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'operate experiment'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('newVarLine', params.StringParam, default="", label="New variable",
                      help="Examples:\n" \
                      "Cp2 ; ug/mL ; numeric ; measurement ; plasma concentration divided by 2\n" \
                      "weight2 ; none ; numeric; label ; Square of the weight of the animal\n" \
                      "sex ; none ; text ; label ; sex of the animal\n")

        form.addParam('operation', params.StringParam, default="", label="Operation",
                      help='You must use available labels or measurements and a valid Python expression. \n'
                        'You have all numpy operations available as np. Examples:\n'
                        '"Female"\n'
                        '0.0\n'
                        'np.log($(C))\n'
                        '2*$(weight)')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runOperate',self.inputExperiment.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runOperate(self, objId):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        tokens = self.newVarLine.get().split(';')
        newVarName = tokens[0].strip()
        self.experiment.variables[newVarName] = PKPDVariable()
        self.experiment.variables[newVarName].parseTokens(tokens)

        varList = []
        operation =self.operation.get()
        idx0 = operation.find("$(")
        while idx0>=0:
            idxF = operation.find(")",idx0)
            varName = operation[(idx0+2):idxF]
            if not varName in varList:
                varList.append(varName)
            idx0 = operation.find("$(",idxF+1)
        parsedOperation=copy.copy(operation)
        for varName in varList:
            exec("parsedOperation=parsedOperation.replace('$(%s)','%s')"%(varName,varName))

        self.printSection("Operating")
        print("Operation performed: %s"%parsedOperation)
        for sampleName in self.experiment.samples:
            print("   Sample " + sampleName)
            sample = self.experiment.samples[sampleName]

            for varName in varList:
                variable = self.experiment.variables[varName]
                if variable.isLabel():
                    exec("%s=sample.getDescriptorValue('%s')"%(varName,varName))
                    if variable.isNumeric():
                        exec("%s=float(%s)"%(varName,varName))
                else:
                    # Measurement or time
                    exec("%s=np.asarray(sample.getValues('%s'),dtype=np.float)"%(varName,varName))
                # exec("print('%s=',%s)"%(varName,varName))

            exec("%s=%s"%(newVarName,parsedOperation))
            # exec ("print('%s=',%s)" % (newVarName, newVarName))

            if self.experiment.variables[newVarName].isLabel():
                exec('sample.setDescriptorValue("%s",%s)'%(newVarName,newVarName))
            else:
                # Measurement or time
                exec('sample.addMeasurementColumn("%s",%s)'%(newVarName,newVarName))

        # Print and save
        self.writeExperiment(self.experiment,self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

    def validate(self):
        retval = []
        tokens = self.newVarLine.get().split(';')
        if len(tokens) != 5:
            retval.append("The new variable description is not well formed: %s"%self.newVarLine)
        return retval