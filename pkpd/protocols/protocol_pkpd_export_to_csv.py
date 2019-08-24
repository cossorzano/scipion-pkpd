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

import numpy as np
import os
import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDExperiment, PKPDVariable
from pyworkflow.protocol.constants import LEVEL_ADVANCED


class ProtPKPDExportToCSV(ProtPKPD):
    """ Export experiment to CSV.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'export to csv'

    ONEMEASUREMENTPERROW=0
    TABULAR=1

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('suffix', params.StringParam, label="File suffix", default="", expertLevel=LEVEL_ADVANCED,
                      help='The output filename is called experiment[Suffix].csv. Do not use spaces. Examples: _new')
        form.addParam('format', params.EnumParam, label="Output format", default=ProtPKPDExportToCSV.ONEMEASUREMENTPERROW,
                      choices=['One measurement per row','Tabular form'])
        form.addParam("tVar", params.StringParam, label="Time variable", default="t", condition="format==1",
                      help="This will be the first variable in the table, it does not need to be a time variable, but it typically is")
        form.addParam("xVar", params.StringParam, label="Measurement variable", default="C", condition="format==1",
                      help="This is the variable whose content is in the table. It does not need to be a measurement variable, but it typically is")

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

        if self.format.get()==ProtPKPDExportToCSV.ONEMEASUREMENTPERROW:
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
                if sample.descriptors:
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
        else:
            allT=[]
            for sampleName,sample in experiment.samples.iteritems():
                allT+=getattr(sample,"measurement_%s"%self.tVar.get())
            allT=set(allT) # Get unique values
            allT=[float(t) for t in allT]
            allT=sorted(allT)

            table=np.full((len(allT),len(experiment.samples)),np.nan)
            j=0
            sortedNames=sorted(experiment.samples)
            for sampleName in sortedNames:
                sample=experiment.samples[sampleName]
                sampleT = getattr(sample, "measurement_%s" % self.tVar.get())
                sampleX = getattr(sample, "measurement_%s" % self.xVar.get())
                for t,x in zip(sampleT,sampleX):
                    tf=float(t)
                    xf=float(x)
                    i=allT.index(tf)
                    table[i][j]=xf
                j+=1

            header=self.tVar.get()+"; "
            for sampleName in sortedNames:
                header+=sampleName+"; "
            print(header)
            fhOut.write(header+"\n")

            for i in range(len(allT)):
                lineToPrint=str(allT[i])+"; "
                for j in range(len(experiment.samples)):
                    if np.isfinite(table[i][j]):
                        lineToPrint += str(table[i][j])
                    lineToPrint+="; "
                print(lineToPrint)
                fhOut.write(lineToPrint+"\n")

        fhOut.close()


    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        retval=[]
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        if format==1 and not self.tVar.get() in experiment.variables.keys():
            retval.append("Cannot find %s among the experiment variables"%self.tVar.get())
        if format==1 and not self.xVar.get() in experiment.variables.keys():
            retval.append("Cannot find %s among the experiment variables"%self.xVar.get())
        return retval

    def _summary(self):
        return ["Output file: %s"%os.path.abspath(self.getFilenameOut())]