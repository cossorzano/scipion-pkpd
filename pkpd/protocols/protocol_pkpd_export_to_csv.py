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
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.data import PKPDExperiment, PKPDVariable
from pyworkflow.protocol.constants import LEVEL_ADVANCED


class ProtPKPDExportToCSV(ProtPKPD):
    """ Export experiment to CSV.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'export to csv'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('suffix', params.StringParam, label="File suffix", default="", expertLevel=LEVEL_ADVANCED,
                      help='The output filename is called experiment[Suffix].csv. Do not use spaces. Examples: _new')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('exportToCSV',self.inputExperiment.get().getObjId())

    #--------------------------- STEPS functions --------------------------------------------
    def getFilenameOut(self):
        preprocessedSuffix = self.suffix.get().replace(' ','_')
        return self._getPath("experiment%s.csv"%preprocessedSuffix)

    def exportToCSV(self, objId):
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)

        self.printSection("Writing "+self.getFilenameOut())
        fhOut = open(self.getFilenameOut(),"w")

        # Prepare header
        header="SampleID; SampleName"
        headerDefaultDict={}
        headerDefaultDict["SampleID"]="NA"
        headerDefaultDict["SampleName"]="NA"
        linePattern="%(SampleID)s; %(SampleName)s"
        listOfVariables = []
        for varRole in [PKPDVariable.ROLE_LABEL, PKPDVariable.ROLE_TIME, PKPDVariable.ROLE_MEASUREMENT]:
            for varName,var in experiment.variables.iteritems():
                if var.role == varRole:
                    header+="; %s"%varName
                    headerDefaultDict[varName]="NA"
                    linePattern+="; %%(%s)s"%varName
                    if var.role == PKPDVariable.ROLE_TIME:
                        listOfVariables.append(varName)
                    elif var.role == PKPDVariable.ROLE_MEASUREMENT:
                        listOfVariables.append(varName)
        print(header)
        fhOut.write(header+"\n")

        # Print all samples
        counter=1
        for sampleName,sample in experiment.samples.iteritems():
            sampleDict=headerDefaultDict.copy()
            sampleDict["SampleID"]=counter
            sampleDict["SampleName"]=sampleName
            for descriptor,value in sample.descriptors.iteritems():
                sampleDict[descriptor]=str(value)
            for i in range(sample.getNumberOfMeasurements()):
                lineDict=sampleDict.copy()
                for varName in listOfVariables:
                    aux = getattr(sample,"measurement_%s"%varName)
                    lineDict[varName]=str(aux[i])
                lineToPrint=linePattern%lineDict
                fhOut.write(lineToPrint+"\n")
                print(lineToPrint)
            counter+=1
        fhOut.close()


    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        return ["Output file: %s"%self.getFilenameOut()]