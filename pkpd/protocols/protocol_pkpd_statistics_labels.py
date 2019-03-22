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
from pkpd.pkpd_units import strUnit
import os
import numpy as np
import scipy.stats

class ProtPKPDStatisticsLabel(ProtPKPD):
    """ Calculate statistics of the labels\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'statistics labels'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runStatistcs',self.inputExperiment.get().getObjId())

    #--------------------------- STEPS functions --------------------------------------------
    def runStatistcs(self, objId):
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        self.printSection("Calculating statistical descriptors")
        fnStatistics = self._getPath("statistics.txt")
        fhOut = open(fnStatistics,'w')

        for varName, variable in experiment.variables.iteritems():
            if variable.role == PKPDVariable.ROLE_LABEL:
                varValues = experiment.getSubGroupLabels("True",varName)
                if variable.varType == PKPDVariable.TYPE_TEXT:
                    self.doublePrint(fhOut,"%s ------------"%varName)
                    counter = {}
                    for value in varValues:
                        if value in counter:
                            counter[value]+=1
                        else:
                            counter[value]=1
                    for value in counter:
                        self.doublePrint(fhOut,"Value=%s Total count=%d (%f%%)"%(value,counter[value],100*float(counter[value])/len(varValues)))
                elif variable.varType == PKPDVariable.TYPE_NUMERIC:
                    self.doublePrint(fhOut,"%s [%s] ------------"%(varName,strUnit(variable.units.unit)))
                    varValues = np.asarray(varValues,np.double)
                    self.doublePrint(fhOut,"Number of observations= %d"%len(varValues))
                    self.doublePrint(fhOut,"Range=    [%f,%f]"%(np.min(varValues),np.max(varValues)))
                    self.doublePrint(fhOut,"Mean=     %f"%np.mean(varValues))
                    self.doublePrint(fhOut,"StdDev=   %f"%np.std(varValues))
                    self.doublePrint(fhOut,"Skewness= %f"%scipy.stats.skew(varValues))
                    self.doublePrint(fhOut,"Kurtosis= %f"%scipy.stats.kurtosis(varValues))
                    self.doublePrint(fhOut,"Quantile 5%%= %f"%np.percentile(varValues,5))
                    self.doublePrint(fhOut,"Quantile 25%%= %f"%np.percentile(varValues,25))
                    self.doublePrint(fhOut,"Quantile 50%%= %f"%np.percentile(varValues,50))
                    self.doublePrint(fhOut,"Quantile 75%%= %f"%np.percentile(varValues,75))
                    self.doublePrint(fhOut,"Quantile 95%%= %f"%np.percentile(varValues,95))
                self.doublePrint(fhOut," ")

        fhOut.close()

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg = []
        fnStatistics = self._getPath("statistics.txt")
        if os.path.exists(fnStatistics):
            self.addFileContentToMessage(msg,fnStatistics)
        return msg
