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

import os
import sys

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD, addDoseToForm
from pyworkflow.em.data import PKPDExperiment, PKPDVariable, PKPDDose, PKPDVia, PKPDSample
from pyworkflow.utils import copyFile


class ProtPKPDImportFromText(ProtPKPD):
    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form, type):
        form.addSection('Input')
        inputFileHelp = "Specify a path to desired %s file.\n"%type
        inputFileHelp += "You may specify missing values with NA. You may also use LLOQ and ULOQ (Lower and Upper limit of quantification)\n"
        inputFileHelp += "to specify measurements that are below or above these limits"
        if type=="CSV":
            inputFileHelp += "The field separator must be a semicolon (;), decimal point must be a dot (.).\n"
            inputFileHelp += "The first row must contain the variable names and one of them must be SampleName which will serve as identifier."
        form.addParam('inputFile', params.PathParam,
                      label="File path", allowsNull=False, help=inputFileHelp)
        form.addParam('title', params.StringParam, label="Title", default="My experiment")
        form.addParam('comment', params.StringParam, label="Comment", default="")
        form.addParam('variables', params.TextParam, height=8, width=80, label="Variables", default="",
                      help="Structure: [Variable Name] ; [Units] ; [Type] ; [Role] ; [Comment]\n"\
                           "The variable name should have no space or special character\n"\
                           "Valid units are: h, mg, ug, ug/mL, ...\n"\
                           "Type is either numeric or text\n"\
                           "Role is either time, label or measurement\n"\
                           "The comment may be empty\n"\
                           "\nIt is important that there are three semicolons (none of them may be missing even if the comment is not present).\n"\
                           "Examples:\n"\
                           "t ; h ; numeric ; time ; \n"\
                           "Cp ; ug/mL ; numeric ; measurement ; plasma concentration\n"\
                           "weight ; g ; numeric; label ; weight of the animal\n"\
                           "sex ; none ; text ; label ; sex of the animal\n")
        form.addParam('vias', params.TextParam, height=8, width=80, label="Vias",
                      help="[ViaName]; [ViaType]; [tlag]; [bioavailability]"\
                       "Valid ViaTypes are: iv (intravenous), ev0 (extra-vascular order 0), ev1 (extra-vascular order 1), \n"\
                       "     ev01 (extra-vascular first order 0 and then order 1), evFractional (extra-vascular fractional order)\n"\
                       "Optional parameters are tlag (e.g. tlag=0)\n"\
                       "   and bioavailability (e.g. bioavailability=0.8)\n"\
                       "Examples:\n"\
                       "Intravenous; iv\n"\
                       "Oral; ev1; tlag; bioavailability=1\n")
        addDoseToForm(form)
        form.addParam('dosesToSamples', params.TextParam, height=5, width=70, label="Assign doses to samples", default="",
                      help="Structure: [Sample Name] ; [DoseName1,DoseName2,...] \n"\
                           "The sample name should have no space or special character\n"\
                           "\nIt is important that there is one semicolon.\n"\
                           "Examples:\n"\
                           "FemaleRat1 ; Bolus0,Bolus1,Infusion0\n")


    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep',self.inputFile.get())

    #--------------------------- STEPS functions --------------------------------------------
    def readTextFile(self):
        pass

    def createOutputStep(self, objId):
        fnFile = os.path.basename(self.inputFile.get())
        copyFile(self.inputFile.get(),self._getPath(fnFile))

        self.experiment = PKPDExperiment()
        self.experiment.general["title"]=self.title.get()
        self.experiment.general["comment"]=self.comment.get()

        ok = True

        # Read the variables
        self.listOfVariables = []
        for line in self.variables.get().replace('\n',';;').split(';;'):
            tokens = line.split(';')
            if len(tokens)!=5:
                print("Skipping variable: ",line)
                ok = False
                continue
            varname = tokens[0].strip()
            self.listOfVariables.append(varname)
            self.experiment.variables[varname] = PKPDVariable()
            self.experiment.variables[varname].parseTokens(tokens)

        # Read vias
        for line in self.vias.get().replace('\n',';;').split(';;'):
            if line!="":
                tokens = line.split(';')
                if len(tokens)<2:
                    print("Skipping via: ",line)
                    ok = False
                    continue
                vianame = tokens[0].strip()
                self.experiment.vias[vianame] = PKPDVia()
                self.experiment.vias[vianame].parseTokens(tokens)

        # Read the doses
        for line in self.doses.get().replace('\n',';;').split(';;'):
            if line!="":
                tokens = line.split(';')
                if len(tokens)<5:
                    print("Skipping dose: ",line)
                    ok = False
                    continue
                dosename = tokens[0].strip()
                self.experiment.doses[dosename] = PKPDDose()
                self.experiment.doses[dosename].parseTokens(tokens,self.experiment.vias)

        # Read the sample doses
        for line in self.dosesToSamples.get().replace('\n',';;').split(';;'):
            try:
                tokens = line.split(';')
                samplename = tokens[0].strip()
                if len(tokens)>1:
                    tokens[1]="dose="+tokens[1]
                self.experiment.samples[samplename] = PKPDSample()
                self.experiment.samples[samplename].parseTokens(tokens, self.experiment.variables, self.experiment.doses,
                                                                self.experiment.groups)
            except Exception as e:
                ok = False
                print("Problem with line: ",line,str(e))

        if ok:
            # Read the measurements
            self.readTextFile()
            self.experiment.write(self._getPath("experiment.pkpd"))
            self.experiment._printToStream(sys.stdout)
            self._defineOutputs(outputExperiment=self.experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        return ["Input file: %s"%self.inputFile.get()]

class ProtPKPDImportFromCSV(ProtPKPDImportFromText):
    """ Import experiment from CSV.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'import from csv'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        ProtPKPDImportFromText._defineParams(self,form,"CSV")

    def readTextFile(self):
        fh=open(self.inputFile.get())
        lineNo = 1
        for line in fh.readlines():
            tokens = line.split(';')
            if len(tokens)==0:
                continue
            if lineNo==1:
                listOfVariables=[]
                listOfSkips=[]
                iSampleName=-1
                varNo=0
                for token in tokens:
                    varName=token.strip()
                    if varName=="SampleName":
                        iSampleName=varNo
                    listOfVariables.append(varName)
                    listOfSkips.append(not (varName in self.experiment.variables))
                    varNo+=1
                if iSampleName==-1:
                    raise Exception("Cannot find the SampleName in: %s\n"%line)
            else:
                if len(tokens)!=len(listOfSkips):
                    print("Skipping line: %s"%line)
                    print("   It does not have the same number of values as the header")
                varNo = 0
                sampleName = tokens[iSampleName].strip()
                if not sampleName in self.experiment.samples:
                    print("Skipping sample: The sample %s does not have a dose"%sampleName)
                    continue
                samplePtr=self.experiment.samples[sampleName]
                for skip in listOfSkips:
                    if not skip:
                        varName = listOfVariables[varNo]
                        varRole = self.experiment.variables[listOfVariables[varNo]].role
                        if varRole == PKPDVariable.ROLE_LABEL:
                            samplePtr.descriptors[listOfVariables[varNo]] = tokens[varNo]
                        else:
                            ok=True
                            if tokens[varNo]=="NA":
                                ok = (varRole != PKPDVariable.ROLE_TIME)
                            if ok:
                                if not hasattr(samplePtr,"measurement_%s"%varName):
                                    setattr(samplePtr,"measurement_%s"%varName,[])
                                    samplePtr.measurementPattern.append(varName)
                                exec("samplePtr.measurement_%s.append('%s')"%(varName,tokens[varNo].strip()))
                            else:
                                raise Exception("Time measurements cannot be NA")
                    varNo+=1
            lineNo+=1

def getSampleNamesFromCSVfile(fnCSV):
    sampleNames = []
    fh=open(fnCSV)
    lineNo = 1
    for line in fh.readlines():
        tokens = line.split(';')
        if len(tokens)==0:
            continue
        if lineNo==1:
            iSampleName=-1
            varNo = 0
            for token in tokens:
                varName=token.strip()
                if varName=="SampleName":
                    iSampleName=varNo
                    break
                varNo += 1
            if iSampleName==-1:
                fh.close()
                return
        else:
            if len(tokens)>iSampleName:
                sampleName = tokens[iSampleName]
                if not sampleName in sampleNames:
                    sampleNames.append(sampleName.strip())
        lineNo+=1
    fh.close()
    return sampleNames

def getVarNamesFromCSVfile(fnCSV):
    varNames = []
    fh=open(fnCSV)
    for line in fh.readlines():
        tokens = line.split(';')
        if len(tokens)==0:
            continue
        for token in tokens:
            varNames.append(token.strip())
        break
    fh.close()
    return varNames
