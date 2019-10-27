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
from itertools import izip
import math
import numpy as np
import os
import openpyxl
from .pkpd_units import (PKPDUnit, convertUnits, changeRateToMinutes,
                        changeRateToWeight)
import pyworkflow as pw
import pyworkflow.utils as pwutils
from pyworkflow.em.data import *
from pyworkflow.install.funcs import *
from .utils import writeMD5, verifyMD5, excelWriteRow, excelFillCells, excelAdjustColumnWidths, computeXYmean
from .biopharmaceutics import PKPDDose, PKPDVia, DrugSource, createDeltaDose, createVia

class PKPDVariable:
    TYPE_NUMERIC = 1000
    TYPE_TEXT = 1001
    TYPE_NAMES = {TYPE_NUMERIC: 'numeric',
                  TYPE_TEXT: 'text'}

    ROLE_TIME = 1010
    ROLE_MEASUREMENT = 1011
    ROLE_LABEL = 1012
    ROLE_NAMES = {ROLE_TIME: 'time',
                  ROLE_MEASUREMENT: 'measurement',
                  ROLE_LABEL: 'label'}

    def __init__(self):
        self.varName = None
        self.varType = None
        self.displayString = ""
        self.role = None
        self.comment = ""
        self.units = None

    def parseTokens(self, tokens):
        # t ; h ; numeric[%f] ; time ;

        strippedTokens=[]
        for token in tokens:
            strippedTokens.append(token.replace(';',''))
        tokens=strippedTokens

        # Get name
        self.varName = tokens[0].strip()

        # Get units
        unitString = tokens[1].strip()
        self.units = PKPDUnit(unitString)
        if not self.units.unit:
            raise Exception("Unrecognized unit: %s"%unitString)

        # Get type and display
        typeString = tokens[2].strip()
        self.displayString = ""
        leftBracket=typeString.find('[')
        rightBracket=typeString.find(']')
        if leftBracket!=-1 and rightBracket!=-1:
            self.displayString = typeString[leftBracket+1:rightBracket]
            typeString = typeString[0:leftBracket]
        typeString = typeString.lower()
        if typeString=="numeric":
            self.varType = PKPDVariable.TYPE_NUMERIC
            if self.displayString=="":
                self.displayString="%f"
        elif typeString=="text":
            self.varType = PKPDVariable.TYPE_TEXT
            if self.displayString=="":
                self.displayString="%s"
        else:
            raise Exception("Unrecognized type: %s"%typeString)

        # Get role
        roleString = tokens[3].strip().lower()
        if roleString=="time":
            self.role = PKPDVariable.ROLE_TIME
        elif roleString=="measurement":
            self.role = PKPDVariable.ROLE_MEASUREMENT
        elif roleString=="label":
            self.role = PKPDVariable.ROLE_LABEL
        else:
            raise Exception("Unrecognized role: %s"%roleString)

        # Get comment
        self.comment = tokens[4].strip()

    def _printToStream(self,fh):
        displayString = self.displayString.replace("%%","%%%%")
        if displayString=="":
            if self.varType == PKPDVariable.TYPE_NUMERIC:
                displayString="%f"
            elif self.varType == PKPDVariable.TYPE_TEXT:
                displayString="%s"

        fh.write("%s ; %s ; %s[%s] ; %s ; %s\n" % (self.varName,
                                                   self.getUnitsString(),
                                                   self.getTypeString(),
                                                   displayString,
                                                   self.getRoleString(),
                                                   self.comment))

    def _printToExcel(self,wb,row, col=1):
        displayString = self.displayString.replace("%%","%%%%")
        if displayString=="":
            if self.varType == PKPDVariable.TYPE_NUMERIC:
                displayString="%f"
            elif self.varType == PKPDVariable.TYPE_TEXT:
                displayString="%s"
        excelWriteRow([self.varName, self.getUnitsString(),
                       "%s[%s]"%(self.getTypeString(),displayString),
                       self.getRoleString(), self.comment], wb, row, col)
        return row+1

    def getTypeString(self):
        return self.TYPE_NAMES.get(self.varType, '')

    def getRoleString(self):
        return self.ROLE_NAMES.get(self.role, '')

    def getUnitsString(self):
        return self.units._toString()

    def isNumeric(self):
        return self.varType == self.TYPE_NUMERIC

    def isLabel(self):
        return self.role == self.ROLE_LABEL

    def isMeasurement(self):
        return self.role == self.ROLE_MEASUREMENT

    def isTime(self):
        return self.role == self.ROLE_TIME

    def getLabel(self):
        return "%s [%s]" % (self.varName, self.getUnitsString())


class PKPDSample:
    def __init__(self):
        self.sampleName = ""
        self.variableDictPtr = None
        self.doseDictPtr = None
        self.doseList = []
        self.groupList = []
        self.descriptors = None
        self.measurementPattern = None

    def parseTokens(self,tokens,variableDict,doseDict,groupDict):
        # FemaleRat1; dose=Dose1[,Dose2]; weight=207; [group=Group1,Group2]

        # Keep a pointer to variableDict and doseDict
        self.variableDictPtr = variableDict
        self.doseDictPtr = doseDict

        # Get name
        self.sampleName = tokens[0].strip()

        if len(tokens)>1:
            # Get rest of variables
            self.descriptors = {}
            for n in range(1,len(tokens)):
                if '=' in tokens[n]:
                    varTokens = tokens[n].split('=')
                    varName  = varTokens[0].strip()
                    varValue = varTokens[1].strip()
                    if varName=="dose":
                        for doseName in varValue.split(','):
                            doseName=doseName.strip()
                            if doseName in doseDict:
                                self.doseList.append(doseName)
                            else:
                                raise Exception("Unrecognized dose %s"%doseName)
                    elif varName=="group":
                        for groupName in varValue.split(','):
                            groupName=groupName.strip()
                            if groupName in groupDict.keys():
                                self.groupList.append(groupName)
                            else:
                                raise Exception("Unrecognized group %s"%groupName)
                    else:
                        if varName in variableDict:
                            varPtr = variableDict[varName]
                            if varPtr.role != PKPDVariable.ROLE_LABEL:
                                raise Exception("Samples can only use role variables")
                            self.descriptors[varName] = varValue

        self.measurementPattern = []

    def getSampleName(self):
        return self.sampleName

    def interpretDose(self):
        self.parsedDoseList = []
        firstUnit = None
        for doseName in sorted(self.doseList):
            dose = copy.copy(self.doseDictPtr[doseName])
            dose.doseAmount = self.evaluateExpression(dose.doseAmount)
            dose.changeTimeUnitsToMinutes()
            dose.prepare()

            if dose.doseType == PKPDDose.TYPE_INFUSION:
                dose.doseAmount, dose.dunits.unit = changeRateToMinutes(dose.doseAmount, dose.dunits.unit)

            if firstUnit==None:
                firstUnit = dose.dunits.unit
            else:
                dose.doseAmount = convertUnits(dose.doseAmount, dose.dunits.unit, firstUnit)
            self.parsedDoseList.append(dose)
        # if len(self.parsedDoseList)==0:
        #     raise Exception("Cannot find any useful dose")

    def isDoseABolus(self):
        if len(self.parsedDoseList)!=1:
            return False
        return self.parsedDoseList[0].isDoseABolus()

    def getDoseAt(self,t0,dt=0.5):
        doseAmount = 0.0
        for dose in self.parsedDoseList:
            doseAmount += dose.getDoseAt(t0,dt)
        return doseAmount

    def getCumulatedDose(self,t0,tF):
        return self.getDoseAt(t0,tF-t0)

    def getDoseUnits(self):
        if len(self.parsedDoseList)>0:
            if self.parsedDoseList[0].doseType==PKPDDose.TYPE_INFUSION:
                return changeRateToWeight(self.parsedDoseList[0].dunits.unit)
            else:
                return self.parsedDoseList[0].dunits.unit
        else:
            return PKPDUnit.UNIT_NONE

    def addMeasurementPattern(self,tokens):
        self.measurementPattern = []
        for n in range(1,len(tokens)):
            varName = tokens[n].strip()
            if varName in self.variableDictPtr:
                self.measurementPattern.append(varName)
                setattr(self,"measurement_%s"%varName,[])
            else:
                raise Exception("Unrecognized variable %s"%varName)

    def addMeasurement(self,line):
        line = line.strip()
        if line=="":
            return
        tokens = line.split()
        if len(tokens)<len(self.measurementPattern):
            raise Exception("Not enough values to fill measurement pattern")
        for n in range(0,len(tokens)):
            ok=True
            varName = self.measurementPattern[n]
            if tokens[n]=="NA" or tokens[n]=="ULOQ" or tokens[n]=="LLOQ":
                ok = (self.variableDictPtr[varName].role != PKPDVariable.ROLE_TIME)
            if ok:
                exec("self.measurement_%s.append('%s')"%(varName,tokens[n]))
            else:
                raise Exception("Time measurements cannot be NA")

    def addMeasurementColumn(self,varName,values):
        if self.measurementPattern is None:
            self.measurementPattern = []
        if not varName in self.measurementPattern:
            self.measurementPattern.append(varName)
        setattr(self, "measurement_%s"%varName, [])
        if type(values)==list:
            for value in values:
                exec("self.measurement_%s.append('%s')"%(varName,str(float(value))))
        elif type(values)==np.ndarray:
            for i in range(values.size):
                exec("self.measurement_%s.append('%s')"%(varName,str(float(values[i]))))

    def getNumberOfVariables(self):
        return len(self.measurementPattern)

    def getNumberOfMeasurements(self):
        return len(getattr(self,"measurement_%s"%self.measurementPattern[0]))

    def _printToStream(self,fh):
        fh.write("%s"%self.sampleName)
        if self.doseList:
            fh.write("; dose=%s"%(",".join(self.doseList)))
        if self.groupList:
            fh.write("; group=%s"%(",".join(self.groupList)))
        if self.descriptors:
            descriptorString = ""
            for key in sorted(self.descriptors.keys()):
                descriptorString +="; %s=%s"%(key,self.descriptors[key])
            fh.write(" %s"%descriptorString)
        fh.write("\n")

    def _printToExcel(self,wb,row,col=1):
        toPrint=[self.sampleName]
        if self.doseList:
            toPrint.append("dose=%s"%(",".join(self.doseList)))
        if self.groupList:
            toPrint.append("group=%s"%(",".join(self.groupList)))
        if self.descriptors:
            for key in sorted(self.descriptors.keys()):
                toPrint.append("%s=%s"%(key,self.descriptors[key]))
        excelWriteRow(toPrint,wb,row,col)
        return row+1

    def _printMeasurements(self,fh):
        patternString = ""
        for n in range(0,len(self.measurementPattern)):
            patternString += "; %s"%self.measurementPattern[n]
        fh.write("%s %s\n"%(self.sampleName,patternString))
        if len(self.measurementPattern)>0:
            aux=getattr(self,"measurement_%s"%self.measurementPattern[0])
            N = len(aux)
            for i in range(0,self.getNumberOfMeasurements()):
                lineString = ""
                for n in range(0,len(self.measurementPattern)):
                    aux=getattr(self,"measurement_%s"%self.measurementPattern[n])
                    lineString += aux[i]+" "
                fh.write("%s\n"%lineString)
        fh.write("\n")

    def _printMeasurementsToExcel(self,wb,row):
        toPrint = [self.sampleName]
        for n in range(0,len(self.measurementPattern)):
            toPrint.append(self.measurementPattern[n])
        excelWriteRow(toPrint,wb,row); row+=1
        if len(self.measurementPattern)>0:
            aux=getattr(self,"measurement_%s"%self.measurementPattern[0])
            N = len(aux)
            for i in range(0,self.getNumberOfMeasurements()):
                toPrint = [""]
                for n in range(0,len(self.measurementPattern)):
                    aux=getattr(self,"measurement_%s"%self.measurementPattern[n])
                    try:
                        toPrint.append(float(aux[i]))
                    except:
                        if isinstance(aux[i],basestring):
                            toPrint.append(aux[i])
                        else:
                            toPrint.append("")
                excelWriteRow(toPrint,wb,row); row+=1
        return row+1

    def getRange(self, varName):
        if varName not in self.measurementPattern:
            return [None, None]
        else:
            aux = getattr(self,"measurement_%s"%varName)
            aux = [x for x in aux if x != "NA" and x!="LLOQ" and x!="ULOQ" and x!="None"]
            x = np.asarray(aux, dtype=np.double)
            return [x.min(),x.max()]

    def getValues(self, varName):
        if type(varName)==list:
            retval=[]
            for vName in varName:
                if vName not in self.measurementPattern:
                    retval.append(None)
                else:
                    retval.append(getattr(self,"measurement_%s"%vName))
            return retval
        else:
            if varName not in self.measurementPattern:
                return None
            else:
                return getattr(self,"measurement_%s"%varName)

    def setValues(self, varName, varValues):
        setattr(self,"measurement_%s"%varName,varValues)

    def getXYValues(self,varNameX,varNameY):
        xl = []
        yl = []
        xs = self.getValues(varNameX)
        if type(varNameY)==list:
            ys = self.getValues(varNameY)
            for ysi in ys:
                xPartial =[]
                yPartial = []
                for x, y in izip(xs, ysi):
                    if x != "NA" and x!="LLOQ" and y!="ULOQ" and y != "NA" and y!= "LLOQ" and y!="ULOQ":
                        xPartial.append(float(x))
                        yPartial.append(float(y))
                xl.append(np.array(xPartial))
                yl.append(np.array(yPartial))
        else:
            ys = self.getValues(varNameY)
            xPartial =[]
            yPartial = []
            for x, y in izip(xs, ys):
                if x != "NA" and x!="LLOQ" and y!="ULOQ" and y != "NA" and y!= "LLOQ" and y!="ULOQ" and \
                   x!= "None" and y!="None":
                    xPartial.append(float(x))
                    if "[" in y:
                        y=y.replace("[","").replace("]","")
                    yPartial.append(float(y))
            xl.append(np.array(xPartial))
            yl.append(np.array(yPartial))
        return xl, yl

    def getSampleMeasurements(self):
        return [PKPDSampleMeasurement(self,n) for n in range(0,self.getNumberOfMeasurements())]

    def substituteValuesInExpression(self, expression, prefix=""):
        expressionPython = copy.copy(expression)
        if self.descriptors is not None:
            for key, variable in self.variableDictPtr.iteritems():
                if key in self.descriptors:
                    value = self.descriptors[key]
                    if value=="NA" or value=="LLOQ" or value=="ULOQ":
                        expressionPython="None"
                        break
                    else:
                        if variable.varType == PKPDVariable.TYPE_NUMERIC:
                            expressionPython = expressionPython.replace("$%s(%s)"%(prefix,key),"%f"%float(value))
                        else:
                            expressionPython = expressionPython.replace("$%s(%s)"%(prefix,key),"'%s'"%value)
        return expressionPython

    def evaluateExpression(self, expression, prefix=""):
        expressionPython=self.substituteValuesInExpression(expression,prefix)
        return eval(expressionPython, {"__builtins__" : {"True": True, "False": False, "None": None} }, {})

    def evaluateParsedExpression(self, parsedOperation, varList):
        for varName in varList:
            variable = self.variableDictPtr[varName]
            if variable.isLabel():
                exec ("%s=self.getDescriptorValue('%s')" % (varName, varName))
                if variable.isNumeric():
                    exec ("%s=float(%s)" % (varName, varName))
            else:
                # Measurement or time
                exec ("%s=np.asarray(self.getValues('%s'),dtype=np.float)" % (varName, varName))
        aux = None
        exec ("aux=%s" % parsedOperation)
        return aux

    def getVariableValues(self, varList):
        varDict = {}
        for varName in varList:
            var = self.variableDictPtr[varName]
            if var.isLabel():
                varDict[varName]=self.descriptors[varName]
            elif var.isMeasurement():
                varDict[varName]=getattr(self,"measurement_%s"%varName)
        return varDict

    def getDescriptorValue(self,descriptorName):
        if self.descriptors is None:
            return None
        if descriptorName in self.descriptors.keys():
            return self.descriptors[descriptorName]
        else:
            return None

    def setDescriptorValue(self, descriptorName, descriptorValue):
        self.descriptors[descriptorName] = descriptorValue


class PKPDSampleMeasurement():
    def __init__(self, sample, n):
        self.sample = sample
        self.n = n

    def getValues(self):
        values = []
        for i in range(0,self.sample.getNumberOfVariables()):
            aux=getattr(self.sample,"measurement_%s"%self.sample.measurementPattern[i])
            values.append(aux[self.n])
        return values


class PKPDGroup():
    def __init__(self, groupName):
        self.groupName = groupName
        self.sampleList = []

    def getSamplesString(self):
        return ",".join(self.sampleList)


class PKPDExperiment(EMObject):
    READING_GENERAL = 1
    READING_VARIABLES = 2
    READING_VIAS = 3
    READING_DOSES = 4
    READING_GROUPS = 5
    READING_SAMPLES = 6
    READING_MEASUREMENTS = 7
    READING_A_MEASUREMENT = 8

    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.fnPKPD = String()
        self.infoStr = String()
        self.general = {}
        self.variables = {}
        self.samples = {}
        self.doses = {}
        self.vias = {}
        self.groups = {}

    def __str__(self):
        if not self.infoStr.hasValue():
            self.load(fullRead=False)
            self.infoStr.set("variables: %d, samples: %d"
                             % (len(self.variables), len(self.samples)))
        return self.infoStr.get()

    def load(self, fnExperiment="", verifyIntegrity=True, fullRead=True):
        if fnExperiment!="":
            self.fnPKPD.set(fnExperiment)
        if verifyIntegrity and not verifyMD5(self.fnPKPD.get()):
            raise Exception("The file %s has been modified since its creation"%self.fnPKPD.get())
        if self.fnPKPD.get() is None:
            return
        fh=open(self.fnPKPD.get(),'r')
        if not fh:
            raise Exception("Cannot open the file "+self.fnPKPD)

        state=None
        for line in fh.readlines():
            line=line.strip()
            if line=="":
                if state==PKPDExperiment.READING_A_MEASUREMENT:
                    state=PKPDExperiment.READING_MEASUREMENTS
                continue
            if line[0]=='[':
                section = line.split('=')[0].strip().lower()
                if section=="[experiment]":
                    state=PKPDExperiment.READING_GENERAL
                elif section=="[variables]":
                    state=PKPDExperiment.READING_VARIABLES
                elif section=="[vias]":
                    state=PKPDExperiment.READING_VIAS
                elif section=="[doses]":
                    state=PKPDExperiment.READING_DOSES
                elif section=="[groups]":
                    state=PKPDExperiment.READING_GROUPS
                elif section=="[samples]":
                    state=PKPDExperiment.READING_SAMPLES
                elif section=="[measurements]":
                    if fullRead:
                        state=PKPDExperiment.READING_MEASUREMENTS
                    else:
                        break
                else:
                    print("Skipping: ",line)

            elif state==PKPDExperiment.READING_GENERAL:
                tokens = line.split('=')
                lhs = tokens[0].strip().lower()
                rhs = tokens[1].strip()
                self.general[lhs]=rhs

            elif state==PKPDExperiment.READING_VARIABLES:
                tokens = line.split(';')
                if len(tokens)!=5:
                    print("Skipping variable: ",line)
                    continue
                varname = tokens[0].strip()
                self.variables[varname] = PKPDVariable()
                self.variables[varname].parseTokens(tokens)

            elif state==PKPDExperiment.READING_VIAS:
                if line!="":
                    tokens = line.split(';')
                    if len(tokens)<2:
                        print("Skipping via: ",line)
                        continue
                    vianame = tokens[0].strip()
                    self.vias[vianame] = PKPDVia(ptrExperiment=self)
                    self.vias[vianame].parseTokens(tokens)

            elif state==PKPDExperiment.READING_DOSES:
                if line!="":
                    tokens = line.split(';')
                    if len(tokens)!=5:
                        print("Skipping dose: ",line)
                        continue
                    dosename = tokens[0].strip()
                    self.doses[dosename] = PKPDDose()
                    self.doses[dosename].parseTokens(tokens,self.vias)

            elif state==PKPDExperiment.READING_GROUPS:
                if line!="":
                    groupName = line.strip()
                    self.groups[groupName] = PKPDGroup(groupName)

            elif state==PKPDExperiment.READING_SAMPLES:
                if line!="":
                    tokens = line.split(';')
                    samplename = tokens[0].strip()
                    newSample = PKPDSample()
                    newSample.parseTokens(tokens,self.variables, self.doses, self.groups)
                    if not newSample.groupList: # If there is no group, create one for this sample
                        groupName = "__"+newSample.sampleName
                        self.groups[groupName]=PKPDGroup(groupName)
                        newSample.groupList.append(groupName)
                    for groupName in newSample.groupList: # A sample may belong to several groups
                        self.groups[groupName].sampleList.append(newSample.sampleName)
                    self.samples[samplename] = newSample

            elif state==PKPDExperiment.READING_MEASUREMENTS:
                tokens = line.split(';')
                if len(tokens)<3:
                    print("Skipping measurement: ",line)
                    continue
                samplename = tokens[0].strip()
                if samplename in self.samples:
                    self.samples[samplename].addMeasurementPattern(tokens)
                    state=PKPDExperiment.READING_A_MEASUREMENT
                else:
                    print("Skipping measurement: %s"%line)
            elif state==PKPDExperiment.READING_A_MEASUREMENT:
                self.samples[samplename].addMeasurement(line)

        fh.close()

    def write(self, fnExperiment, writeToExcel=True):
        fh=open(fnExperiment,'w')
        self._printToStream(fh)
        fh.close()
        self.fnPKPD.set(fnExperiment)
        writeMD5(fnExperiment)
        self.infoStr.set("variables: %d, samples: %d" % (len(self.variables), len(self.samples)))
        if writeToExcel:
            self.writeToExcel(os.path.splitext(fnExperiment)[0]+".xlsx")

    def _printToStream(self,fh):
        fh.write("[EXPERIMENT] ===========================\n")
        for key, value in self.general.iteritems():
            fh.write("%s = %s\n"%(key,value))
        fh.write("\n")

        fh.write("[VARIABLES] ============================\n")
        for key in sorted(self.variables.keys()):
            self.variables[key]._printToStream(fh)
        fh.write("\n")

        fh.write("[VIAS] ================================\n")
        for key in sorted(self.vias.keys()):
            self.vias[key]._printToStream(fh)
        fh.write("\n")

        fh.write("[DOSES] ================================\n")
        for key in sorted(self.doses.keys()):
            self.doses[key]._printToStream(fh)
        fh.write("\n")

        fh.write("[GROUPS] ================================\n")
        for groupName in sorted(self.groups.keys()):
            fh.write("%s\n"%groupName)
        fh.write("\n")

        fh.write("[SAMPLES] ================================\n")
        for key in sorted(self.samples.keys()):
            self.samples[key]._printToStream(fh)
        fh.write("\n")

        fh.write("[MEASUREMENTS] ===========================\n")
        for key in sorted(self.samples.keys()):
            self.samples[key]._printMeasurements(fh)
        fh.write("\n")

    def writeToExcel(self, fnXls):
        wb = openpyxl.Workbook()
        wb.active.title = "Experiment"

        currentRow = 1
        excelWriteRow("EXPERIMENT",wb,currentRow,bold=True); excelFillCells(wb,currentRow); currentRow+=1
        for key, value in self.general.iteritems():
            excelWriteRow([key,value],wb,currentRow); currentRow+=1
        currentRow+=1

        excelWriteRow("VARIABLES",wb,currentRow,bold=True); excelFillCells(wb,currentRow); currentRow+=1
        for key in sorted(self.variables.keys()):
            currentRow=self.variables[key]._printToExcel(wb,currentRow)
        currentRow+=1

        excelWriteRow("VIAS",wb,currentRow,bold=True); excelFillCells(wb,currentRow); currentRow+=1
        for key in sorted(self.vias.keys()):
            currentRow=self.vias[key]._printToExcel(wb,currentRow)
        currentRow+=1

        excelWriteRow("DOSES",wb,currentRow,bold=True); excelFillCells(wb,currentRow); currentRow+=1
        for key in sorted(self.doses.keys()):
            currentRow=self.doses[key]._printToExcel(wb,currentRow)
        currentRow+=1

        excelWriteRow("GROUPS",wb,currentRow,bold=True); excelFillCells(wb,currentRow); currentRow+=1
        for groupName in sorted(self.groups.keys()):
            excelWriteRow(groupName,wb, currentRow); currentRow+=1
        currentRow+=1

        excelWriteRow("SAMPLES",wb,currentRow,bold=True); excelFillCells(wb,currentRow); currentRow+=1
        for key in sorted(self.samples.keys()):
            currentRow=self.samples[key]._printToExcel(wb,currentRow)
        currentRow+=1

        excelWriteRow("MEASUREMENTS",wb,currentRow,bold=True); excelFillCells(wb,currentRow); currentRow+=1
        for key in sorted(self.samples.keys()):
            currentRow=self.samples[key]._printMeasurementsToExcel(wb,currentRow)
        currentRow+=1

        excelAdjustColumnWidths(wb)
        wb.save(fnXls)

    def getRange(self,varName):
        vmin = None
        vmax = None
        for key, value in self.samples.iteritems():
            vmini, vmaxi = value.getRange(varName)
            if vmin==None or vmini<vmin:
                vmin = vmini
            if vmax==None or vmaxi>vmax:
                vmax = vmaxi
        return [vmin,vmax]

    def sampleSummary(self):
        summary=[]
        for varName, var in self.variables.iteritems():
            if var.role == PKPDVariable.ROLE_LABEL:
                toAdd = varName+": "
                if var.varType==PKPDVariable.TYPE_NUMERIC:
                    listOfValues=[]
                else:
                    listOfValues={}
                for sampleName, sample in self.samples.iteritems():
                    value = sample.descriptors[varName]
                    if var.varType==PKPDVariable.TYPE_NUMERIC:
                        listOfValues.append(float(value))
                    else:
                        if value in listOfValues:
                            listOfValues[value]+=1
                        else:
                            listOfValues[value]=1
                if var.varType==PKPDVariable.TYPE_NUMERIC:
                    listOfValuesNp = np.array(listOfValues)
                    toAdd += " mean=%f std=%f 5%%=%f 25%%=%f 50%%=%f 75%%=%f 95%%=%f"%\
                             (np.mean(listOfValuesNp),np.std(listOfValuesNp),np.percentile(listOfValuesNp,5),\
                              np.percentile(listOfValuesNp,25),np.percentile(listOfValuesNp,50),\
                              np.percentile(listOfValuesNp,75),np.percentile(listOfValuesNp,95))
                else:
                    for value in listOfValues:
                        toAdd += value + "(" + str(listOfValues[value]) + ") "
                summary.append(toAdd)
        return summary

    def getXYMeanValues(self,varNameX,varNameY):
        XYlist = []
        for sampleName, sample in self.samples.iteritems():
            xValues, yValues = sample.getXYValues(varNameX,varNameY)
            XYlist.append((xValues,yValues))
        return computeXYmean(XYlist)

    def getVarUnits(self,varName):
        if varName in self.variables:
            return self.variables[varName].units.unit
        else:
            return PKPDUnit.UNIT_NONE

    def getDoseUnits(self):
        if len(self.samples)==0:
            return PKPDUnit.UNIT_NONE
        listOfSamples = list(self.samples.values())
        listOfSamples[0].interpretDose()
        return listOfSamples[0].getDoseUnits()

    def addParameterToSample(self, sampleName, varName, varUnits, varDescr, varValue, rewrite=False):
        if not varName in self.variables:
            varX = PKPDVariable()
            varX.varName = varName
            varX.varType = PKPDVariable.TYPE_NUMERIC
            varX.displayString = "%f"
            varX.role = PKPDVariable.ROLE_LABEL
            varX.comment = varDescr
            varX.units = PKPDUnit()
            varX.units.unit = varUnits
            self.variables[varName] = varX
        else:
            varPresent = self.variables[varName]
            if varPresent.role!=PKPDVariable.ROLE_LABEL:
                raise Exception("Only labels can be reused (%s)"%varName)

            if rewrite:
                varPresent.comment = varDescr
                varPresent.units.unit = varUnits
            else:
                if varPresent.comment!=varDescr or varPresent.units.unit!=varUnits:
                    raise Exception("%s is already a variable in the experiment with a different purpose"%varName)

        if sampleName in self.samples:
            sample = self.samples[sampleName]
            if sample.descriptors==None:
                sample.descriptors={}
            sample.descriptors[varName] = varValue

    def addLabelToSample(self, sampleName, varName, varDescr, varValue, rewrite=False):
        if not varName in self.variables:
            varX = PKPDVariable()
            varX.varName = varName
            varX.varType = PKPDVariable.TYPE_TEXT
            varX.displayString = "%s"
            varX.role = PKPDVariable.ROLE_LABEL
            varX.comment = varDescr
            varX.units = PKPDUnit()
            varX.units.unit = PKPDUnit.UNIT_NONE
            self.variables[varName] = varX
        else:
            varPresent = self.variables[varName]
            if varPresent.role!=PKPDVariable.ROLE_LABEL:
                raise Exception("Only labels can be reused (%s)"%varName)

            if rewrite:
                varPresent.comment = varDescr
                varPresent.units.unit = PKPDUnit.UNIT_NONE
            else:
                if varPresent.comment!=varDescr:
                    raise Exception("%s is already a variable in the experiment with a different purpose"%varName)

        if sampleName in self.samples:
            sample = self.samples[sampleName]
            if sample.descriptors==None:
                sample.descriptors={}
            sample.descriptors[varName] = varValue

    def getSubGroup(self,condition):
        if condition=="":
            return self.samples
        samplesSubGroup = {}
        for sampleName, sample in self.samples.iteritems():
            if sample.evaluateExpression(condition):
                samplesSubGroup[sampleName] = sample
        return samplesSubGroup

    def getSubGroupLabels(self,condition,labelName):
        subgroupLabels = []
        for sampleName, sample in self.samples.iteritems():
            if condition!="" and sample.evaluateExpression(condition) or condition=="":
                subgroupLabels.append(sample.descriptors[labelName])
        return subgroupLabels

    def getNonBolusDoses(self):
        nonBolusList = []
        for sampleName, sample in self.samples.iteritems():
            sample.interpretDose()
            if not sample.isDoseABolus():
                nonBolusList.append(sampleName)
        return nonBolusList

    def addSampleToGroup(self,groupName,sample):
        if not groupName in self.groups.keys():
            self.groups[groupName] = PKPDGroup(groupName)
        self.groups[groupName].sampleList.append(sample.sampleName)

    def addSample(self,sample):
        self.samples[sample.sampleName] = sample
        for groupName in sample.groupList:
            self.addSampleToGroup(groupName,sample)

    def getTimeVariable(self):
        for varName in self.variables:
            if self.variables[varName].isTime():
                return varName

    def getMeasurementVariables(self):
        retval=[]
        for varName in self.variables:
            if self.variables[varName].isMeasurement():
                retval.append(varName)
        return retval


class PKPDModelBase(object):
    def __init__(self):
        self.fnExperiment = None
        self.parameters = None
        self.parameterUnits = None
        self.xName = None
        self.yName = None
        self.experiment = None

    def setExperiment(self, experiment):
        self.experiment = experiment
        if experiment!=None:
            self.fnExperiment = experiment.fnPKPD

    def setXVar(self, x):
        if not x in self.experiment.variables:
            raise Exception("Cannot find %s as a variable in the experiment"%x)
        self.xName = x
        self.xRange = self.experiment.getRange(x)

    def setYVar(self, y):
        if type(y)==list:
            self.yName = []
            self.yRange = []
            for yi in y:
                if not yi in self.experiment.variables:
                    raise Exception("Cannot find %s as a variable in the experiment"%yi)
                self.yName.append(yi)
                self.yRange.append(self.experiment.getRange(yi))
        else:
            if not y in self.experiment.variables:
                raise Exception("Cannot find %s as a variable in the experiment"%y)
            self.yName = y
            self.yRange = self.experiment.getRange(y)

    def setXYValues(self, x, y):
        self.x = []
        self.y = []
        self.ylog = []
        for n in range(len(x)):
            idx = np.logical_and(np.isfinite(x[n]), np.isfinite(y[n]))
            xidx=x[n][idx]
            yidx=y[n][idx]
            self.x.append(xidx)
            self.y.append(yidx)
            self.ylog.append(np.array([math.log10(yidxi) if yidxi>0 else float("inf") for yidxi in yidx]))

    def getNumberOfParameters(self):
        return len(self.getParameterNames())

    def getDescription(self):
        pass

    def getParameterNames(self):
        pass

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form %s'%self.getModelEquation()]*self.getNumberOfParameters()

    def calculateParameterUnits(self,sample):
        pass

    def rearrange(self,parameters):
        return parameters

    def setParameters(self, parameters):
        self.parameters = self.rearrange(parameters)


class PKPDModelBase2(PKPDModelBase):
    def __init__(self):
        PKPDModelBase.__init__(self)
        self.bounds = None

    def forwardModel(self, parameters, x=None):
        pass

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Variables: "+str(self.getParameterNames()))
        print("Bounds: "+str(self.getBounds()))

    def setSample(self, sample):
        self.sample = sample
        self.Dunits = sample.getDoseUnits()

    def getEquation(self):
        pass

    def getModelEquation(self):
        pass

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form %s'%self.getModelEquation()]*self.getNumberOfParameters()
        pass

    def areParametersSignificant(self, lowerBound, upperBound):
        """
        :param lowerBound and upperBound: a numpy array of parameters
        :return: a list of string with "True", "False", "NA", "Suspicious"
        """
        pass

    def areParametersValid(self, p):
        pass

    def setBounds(self, boundsString):
        self.bounds = None
        if boundsString!="" and boundsString!=None:
            tokens=boundsString.split(';')
            if len(tokens)!=self.getNumberOfParameters():
                raise Exception("The number of bound intervals (%d) does not match the number of parameters (%d)"%\
                                (len(tokens),self.getNumberOfParameters()))
            self.bounds=[]
            for token in tokens:
                values = token.strip().split(',')
                self.bounds.append((float(values[0][1:]),float(values[1][:-1])))

    def getBounds(self):
        return self.bounds

    def setConfidenceInterval(self,lowerBound,upperBound):
        yPredictedBackup = copy.copy(self.yPredicted)
        self.yPredictedLower=[]
        self.yPredictedUpper=[]
        for j in range(len(self.yPredicted)):
            self.yPredictedLower.append(np.copy(self.yPredicted[j]))
            self.yPredictedUpper.append(np.copy(self.yPredicted[j]))
        for i in range(0,int(math.pow(2,self.getNumberOfParameters()))):
            pattern = ("{0:0%db}"%(self.getNumberOfParameters())).format(i)
            p = np.where(np.array(list(pattern))=="1",upperBound,lowerBound)
            p = p.astype(np.float)
            if not self.areParametersValid(p):
                continue
            y = self.forwardModel(p)
            for j in range(len(y)):
                yj=y[j]
                for n in range(len(yj)):
                    if yj[n]<(self.yPredictedLower[j][n]):
                        if yj[n]<0:
                            self.yPredictedLower[j][n]=0
                        else:
                            (self.yPredictedLower[j][n])=yj[n]
                    if yj[n]>(self.yPredictedUpper[j][n]):
                        (self.yPredictedUpper[j][n])=yj[n]
        self.yPredicted = yPredictedBackup

    def setConfidenceIntervalNA(self):
        self.yPredictedUpper = []
        self.yPredictedLower = []
        for y in self.yPredicted:
            self.yPredictedUpper.append(["NA"]*y.shape[0])
            self.yPredictedLower.append(["NA"]*y.shape[0])


class PKPDModel(PKPDModelBase2):
    def prepare(self):
        pass


class PKPDODEModel(PKPDModelBase2):
    def __init__(self):
        PKPDModelBase2.__init__(self)
        self.t0 = None # (min)
        self.tF = None # (min)
        self.deltaT = 0.25 # (min)
        self.drugSource = None
        self.drugSourceImpulse = None
        self.tFImpulse = None
        self.thImpulse = None
        # self.show = False

    def setXYValues(self, x, y):
        if type(x)!=list or (type(x) and type(x[0])!=np.ndarray):
            x = [np.array(x)]*self.getResponseDimension()

        if type(y)!=list or (type(y)==list and type(y[0])!=np.ndarray):
            y = [np.array(y)]

        self.x=[]
        self.y=[]
        self.ylog=[]
        for n in range(self.getResponseDimension()):
            idx = np.logical_and(np.isfinite(x[n]), np.isfinite(y[n]))
            xidx=x[n][idx]
            yidx=y[n][idx]
            self.x.append(xidx)
            self.y.append(yidx)
            self.ylog.append(np.array([math.log10(yidxi) if yidxi>0 else float("inf") for yidxi in yidx]))

    def F(self, t, y):
        return 0

    def G(self, t, dD):
        return 0

    def imposeConstraints(self, yt):
        if type(yt)==np.ndarray:
            yt[yt<0]=0
        elif type(yt)==np.float64:
            if yt<0:
                yt=0

    def H(self, y):
        pass

    def getResponseDimension(self):
        return None

    def getStateDimension(self):
        return None

    def forwardModel(self, parameters, x=None, drugSource=None):
        self.parameters = parameters
        if drugSource is None:
            drugSource=self.drugSource

        # Simulate the system response
        t = self.t0
        Nsamples = int(math.ceil((self.tF-self.t0)/self.deltaT))+1
        if self.getStateDimension()>1:
            yt = np.zeros(self.getStateDimension(),np.double)
            Yt = np.zeros((Nsamples,self.getStateDimension()),np.double)
        else:
            yt = 0.0
            Yt = np.zeros(Nsamples)
        Xt = np.zeros(Yt.shape[0])
        delta_2 = 0.5*self.deltaT
        K = self.deltaT/3
        for i in range(0,Nsamples):
            t = self.t0 + i*self.deltaT # More accurate than t+= self.deltaT
            Xt[i]=t

            # Internal evolution
            # Runge Kutta's 4th order (http://lpsa.swarthmore.edu/NumInt/NumIntFourth.html)
            k1 = self.F(t,yt)
            dD1 = drugSource.getAmountReleasedAt(t,delta_2)
            dyD1 = self.G(t, dD1)
            y1 = yt+k1*delta_2+dyD1
            # print("t=",t," y0=",yt," k1=",k1," dD1=",dD1," dyD1=",dyD1," y1=",y1)

            t_delta_2=t+delta_2
            k2 = self.F(t_delta_2,y1)
            y2 = yt+k2*delta_2+dyD1
            # print("k2=",k2," y2=",y2)

            dD = drugSource.getAmountReleasedAt(t,self.deltaT)
            dyD = self.G(t, dD)
            k3 = self.F(t_delta_2,y2)
            y3 = yt+k3*self.deltaT+dyD
            # print("k3=",k3," dD=",dD," dyD=",dyD," y3=",y3)

            k4 = self.F(t+self.deltaT,y3)
            # y4 = yt+k4*self.deltaT+dyD
            # print("k4=",k4," y4=",y4)

            # Update state
            yt += (0.5*(k1+k4)+k2+k3)*K+dyD
            # print("yt=",yt)
            # print(" ")

            # Make sure it makes sense
            self.imposeConstraints(yt)

            # Apply measurement transformation
            self.H(yt)

            # if self.show:
            #     print("t=%f dD=%s dyD=%s dy=%s"%(t,str(dD),str(dyD),str((0.5*(k1+k4)+k2+k3)*K)))

            # Keep this result and go to next iteration
            if self.getStateDimension()>1:
                Yt[i,:]=yt
            else:
                Yt[i]=yt

        # Get the values at x
        if x is None:
            x = self.x

        self.yPredicted = []
        for j in range(0,self.getResponseDimension()):
            if self.getStateDimension()==1:
                self.yPredicted.append(np.interp(x[j],Xt,Yt))
            else:
                self.yPredicted.append(np.interp(x[j],Xt,Yt[:,j]))
        return self.yPredicted

    def getImpulseResponse(self, parameters, tImpulse):
        if self.tFImpulse is None:
            self.tFImpulse = self.tF

        # Create unit dose
        if self.drugSourceImpulse is None:
            self.drugSourceImpulse = DrugSource()
            dose = createDeltaDose(1.0, via=createVia("Intravenous; iv",self.experiment),
                                   dunits=self.drugSource.getDoseUnits())
            self.drugSourceImpulse.setDoses([dose], 0.0, self.tFImpulse)
        y=self.forwardModel(parameters, [tImpulse], self.drugSourceImpulse)
        return y[0]

    def forwardModelByConvolution(self, parameters, x=None):
        self.parameters = parameters

        # Simulate the system response
        Nsamples = int(math.ceil(self.tF/self.deltaT))+1
        Xt = np.zeros(Nsamples)

        # Get the drug input
        D = copy.copy(Xt)
        for i in range(0,Nsamples):
            t = i*self.deltaT # More accurate than t+= self.deltaT
            Xt[i]=t
            D[i]=self.drugSource.getAmountReleasedAt(t,self.deltaT)

        # Get the model impulse response
        if self.thImpulse is None:
            Nsamples = int(math.ceil(self.tFImpulse/self.deltaT))+1
            self.thImpulse = np.zeros(Nsamples)
            for i in range(0,Nsamples):
                self.thImpulse[i] = i*self.deltaT # More accurate than t+= self.deltaT
        h = self.getImpulseResponse(parameters, self.thImpulse)
        Yt = np.convolve(D,h,'full')[0:len(D)]

        # Get the values at x
        if x is None:
            x = self.x

        self.yPredicted = []
        for j in range(0,self.getResponseDimension()):
            self.yPredicted.append(np.interp(x[j],Xt,Yt))
        return self.yPredicted

    def printOtherParameterization(self):
        pass


class PKPDOptimizer:
    def __init__(self,model,fitType,goalFunction="RMSE"):
        self.model = model
        self.fitType = fitType
        self.Nevaluations = 0
        self.bestRmse=1e38

        self.yTarget = [np.array(yi, dtype=np.float32) for yi in model.y]
        self.yTargetLogs = [np.log10(yi) for yi in self.yTarget]

        if fitType=="linear":
            self.takeYLogs = False
            self.takeRelative = False
        elif fitType=="log":
            self.yTarget = self.yTargetLogs
            self.takeYLogs = True
            self.takeRelative = False
        elif fitType=="relative":
            self.takeYLogs = False
            self.takeRelative = True

        self.bounds = model.getBounds()

        if goalFunction=="RMSE":
            self.goalFunction = self.goalRMSE
        else:
            raise Exception("Unknown goal function")

        self.verbose = 1

    def inBounds(self,parameters):
        if self.bounds==None or len(self.bounds)!=len(parameters):
            return True
        for n in range(0,len(parameters)):
            if parameters[n]<self.bounds[n][0] or parameters[n]>self.bounds[n][1]:
                return False
        return True

    def hugeError(self):
        allDiffs = None
        for yTarget in self.yTarget:
            diff = 1e38*np.ones(yTarget.shape)
            if allDiffs is None:
                allDiffs = diff
            else:
                allDiffs = np.concatenate([allDiffs, diff])
        return allDiffs

    def getResiduals(self,parameters):
        if not self.inBounds(parameters):
            return self.hugeError()
        yPredicted = self.model.forwardModel(parameters)

        allDiffs = None
        for y, yTarget, yTargetLog in izip(yPredicted,self.yTarget,self.yTargetLogs):
            if self.takeYLogs:
                diff = np.full(yTarget.shape,np.nan)
                idx = np.logical_and(np.isfinite(y),y>=1e-20)
                if np.sum(idx)<0.8*diff.size:
                    return self.hugeError()
                diff[idx] = yTargetLog[idx]-np.log10(y[idx])
            else:
                diff = yTarget - y
            if self.takeRelative:
                diff = diff/yTarget

            if allDiffs is None:
                allDiffs = diff
            else:
                allDiffs = np.concatenate([allDiffs, diff])

        idx = np.logical_not(np.isfinite(allDiffs))
        allDiffs[idx]=np.nan
        e = allDiffs
        if e.size<parameters.size:
            return self.hugeError()

        rmse = math.sqrt(np.nanmean(np.power(e,2)))
        if rmse<self.bestRmse:
            print("   Best rmse so far=%f"%rmse)
            print("      at x=%s"%str(parameters))
            print("      e=%s"%str(e))
            sys.stdout.flush()
            self.bestRmse=rmse
        elif self.Nevaluations%100==0:
            print("   Neval=%d RMSE=%f"%(self.Nevaluations,rmse))
            sys.stdout.flush()
        self.Nevaluations+=1
        return e

    def goalRMSE(self,parameters):
        e = self.getResiduals(parameters)
        rmse = math.sqrt(np.nanmean(np.power(e,2)))
        return rmse

    def _evaluateQuality(self, x, y, yp):
        # Spiess and Neumeyer, BMC Pharmacology 2010, 10:6
        self.e = None
        yToUse = None
        for yi, ypi in izip(y,yp):
            diff = []
            for yii, ypii in izip(yi,ypi):
                if np.isfinite(yii) and np.isfinite(ypii):
                    diff.append(yii-ypii)
            diff = np.asarray(diff)
            if self.e is None:
                self.e = diff
                yToUse = np.asarray(yi)
            else:
                self.e = np.concatenate([self.e, diff])
                yToUse = np.concatenate([yToUse, yi])


        yToUse[np.logical_not(np.isfinite(yToUse))]=np.nan # Remove infinites
        self.R2 = (1-np.nanvar(diff)/np.nanvar(yToUse))
        n=len(diff) # Number of samples
        p=self.model.getNumberOfParameters()
        if n-p>0:
            self.R2adj = 1-self.R2*(n-1)/(n-p)*(1-self.R2)
        else:
            self.R2adj = -1
        d2=np.nansum(np.multiply(diff,diff))
        if d2>0:
            logL = 0.5*(-n*(math.log(2*math.pi)+1-math.log(n)+math.log(d2)))
        else:
            logL=np.nan
        self.AIC = 2*p-2*logL
        if n-p-1>0:
            self.AICc = self.AIC+2*p*(p+1)/(n-p-1)
        else:
            self.AICc = -1
        self.BIC = p*math.log(n)-2*logL

    def _printFitting(self, x, y, yp):
        self._evaluateQuality(x, y, yp)

        for j in range(len(x)):
            if len(x)>1:
                print("Series %d ---------"%j)
            xj=x[j]
            yj=y[j]
            ypj=yp[j]
            for n in range(0,xj.shape[0]):
                print("%f %f %f %f"%(xj[n],yj[n],ypj[n],yj[n]-ypj[n]))
        print("------------------------")
        print("Mean error = %f"%np.mean(self.e))
        print("Std error = %f"%np.std(self.e))
        print("R2 = %f"%self.R2)
        print("R2adj = %f"%self.R2adj)
        print("AIC = %f"%self.AIC)
        print("AICc(Recommended) = %f"%self.AICc)
        print("BIC = %f"%self.BIC)
        print("------------------------")

    def evaluateQuality(self):
        yPredicted=self.model.forwardModel(self.model.parameters)
        x = copy.copy(self.model.x)
        y = copy.copy(self.model.y)
        if self.fitType=="linear" or self.fitType=="relative":
            yp = yPredicted
        elif self.fitType=="log":
            if type(self.model.y[0])!=list and type(self.model.y[0])!=np.ndarray:
                if type(yPredicted[0])==list or type(yPredicted[0])==np.ndarray:
                    yPredicted=yPredicted[0]
                ylog = smartLog(self.model.y)
                yplog = smartLog(yPredicted)
                idx = np.array(np.where(np.logical_and(np.isfinite(ylog),np.isfinite(yplog))))
                if type(idx[0])==np.ndarray:
                    idx=idx[0]
                xToUse = self.model.x
                if type(self.model.x)==list:
                    xToUse=np.asarray(xToUse)
                if type(ylog)==list:
                    ylog=np.asarray(ylog)
                if type(yplog)==list:
                    yplog=np.asarray(yplog)
                x= xToUse[idx].ravel()
                y= ylog[idx].ravel()
                yp=yplog[idx].ravel()
            else:
                y = [smartLog(yi) for yi in self.model.y]
                yp = [smartLog(yi) for yi in yPredicted]
        self._evaluateQuality(x,y,yp)

    def printFitting(self):
        yPredicted = self.model.forwardModel(self.model.parameters)
        print("==========================================")
        print("X     Y    Ypredicted  Error=Y-Ypredicted ")
        print("==========================================")
        self._printFitting(self.model.x, self.model.y, yPredicted)
        print("==================================================================")
        print("X    log10(Y)  log10(Ypredicted)  Error=log10(Y)-log10(Ypredicted)")
        print("==================================================================")
        if type(self.model.y[0])!=list and type(self.model.y[0])!=np.ndarray:
            if type(yPredicted[0])==list or type(yPredicted[0])==np.ndarray:
                yPredicted=yPredicted[0]
            ylog = smartLog(self.model.y)
            yplog = smartLog(yPredicted)
            idx = np.array(np.where(np.logical_and(np.isfinite(ylog),np.isfinite(yplog))))
            if type(idx[0])==np.ndarray:
                idx=idx[0]
            xToUse = self.model.x
            if type(self.model.x)==list:
                xToUse=np.asarray(xToUse)
            if type(ylog)==list:
                ylog=np.asarray(ylog)
            if type(yplog)==list:
                yplog=np.asarray(yplog)
            self._printFitting(xToUse[idx].ravel(), ylog[idx].ravel(), yplog[idx].ravel())
        else:
            logY = [smartLog(y) for y in self.model.y]
            logYp = [smartLog(y) for y in yPredicted]
            self._printFitting(self.model.x, logY, logYp)


class PKPDDEOptimizer(PKPDOptimizer):
    def optimize(self):
        from scipy.optimize import differential_evolution
        if self.verbose>0:
            print("Optimizing with Differential Evolution (DE), a global optimizer")
        self.optimum = differential_evolution(self.goalFunction, self.model.getBounds(), maxiter=30)
        if self.verbose>0:
            print("Best DE function value: "+str(self.optimum.fun))
            print("Best DE parameters: "+str(self.optimum.x))
        self.model.setParameters(self.optimum.x)
        if self.verbose>0:
            print(self.model.getEquation())
            self.printFitting()
            print(" ")
        return self.optimum


class PKPDLSOptimizer(PKPDOptimizer):
    def optimize(self, ftol=1.49012e-8, xtol=1.49012e-8): # Same values as in minpack.py
        from scipy.optimize import leastsq
        if self.verbose>0:
            print("Optimizing with Least Squares (LS), a local optimizer")
            print("Initial parameters: "+str(self.model.parameters))
        self.optimum, J, self.info, mesg, _ = leastsq(self.getResiduals, self.model.parameters, full_output=True,
                                                      ftol=ftol, xtol=xtol)
            # J is the jacobian C=MSE*inv(J'*J)
        if self.verbose>0:
            print("Best LS function value: "+str(self.goalFunction(self.optimum)))
            print("Best LS parameters: "+str(self.optimum))
            print("Covariance matrix:")
            if J is not None:
                e = self.getResiduals(self.optimum)
                JtJ=np.matmul(np.transpose(J),J)
                if np.linalg.det(JtJ)>1e-6:
                    self.cov_x = np.var(e)*np.linalg.inv(JtJ)
                    print(np.array_str(self.cov_x,max_line_width=120))
                else:
                    self.cov_x = None
                    print("Singular Jacobian, we cannot estimate the covariance")
            else:
                self.cov_x = None
                print("Singular Jacobian, we cannot estimate the covariance")
        self.model.setParameters(self.optimum)
        if self.verbose>0:
            print(self.model.getEquation())
            self.printFitting()
            print(" ")
        return self.optimum

    def setConfidenceInterval(self,confidenceInterval):
        if self.cov_x is not None:
            from scipy.stats import norm
            nstd = norm.ppf(1-(1-confidenceInterval/100)/2)
            perr = np.sqrt(np.clip(np.diag(self.cov_x),0.0,None))
            self.lowerBound = self.optimum-nstd*perr
            self.upperBound = self.optimum+nstd*perr

            self.significance = self.model.areParametersSignificant(self.lowerBound,self.upperBound)
            parameterNames = self.model.getParameterNames()
            print("Confidence intervals %f%% --------------------------"%confidenceInterval)
            print("ParameterName ParameterValue   ParameterConfidenceInterval  IsStatisticallySignificant")
            for n in range(0,len(self.optimum)):
                print("%s %f [%f,%f] %s"%(parameterNames[n],self.optimum[n],self.lowerBound[n],self.upperBound[n],self.significance[n]))

            self.model.setConfidenceInterval(self.lowerBound,self.upperBound)
        else:
            self.lowerBound=["NA"]*len(self.optimum)
            self.upperBound=["NA"]*len(self.optimum)
            self.significance=["NA"]*len(self.optimum)
            self.model.setConfidenceIntervalNA()


class PKPDSampleFit:
    READING_SAMPLEFITTINGS_NAME = 0
    READING_SAMPLEFITTINGS_MODELEQ = 1
    READING_SAMPLEFITTINGS_R2 = 2
    READING_SAMPLEFITTINGS_R2ADJ = 3
    READING_SAMPLEFITTINGS_AIC = 4
    READING_SAMPLEFITTINGS_AICc = 5
    READING_SAMPLEFITTINGS_BIC = 6
    READING_SAMPLEFITTINGS_PARAMETER_BOUNDS = 7
    READING_SAMPLEFITTINGS_SAMPLE_VALUES = 8

    def __init__(self):
        self.sampleName = ""

        # Lists with the sample fits values
        self.x = None
        self.y = None
        self.yp = None
        self.yl = None
        self.yu = None

        self.modelEquation = ""
        self.R2 = 0
        self.R2adj = 0
        self.AIC = 0
        self.AICc = 0
        self.BIC = 0
        self.parameters = None
        self.lowerBound = None
        self.upperBound = None
        self.significance = None

        self.multiOutputSeries = False

    def printForPopulation(self,fh,observations):
        outputStr = ""
        for parameter in self.parameters:
            outputStr += "%f "%parameter
        observations = np.vstack([observations, self.parameters])
        outputStr += " # %f %f %f %f %f"%(self.R2,self.R2adj,self.AIC,self.AICc,self.BIC)
        fh.write(outputStr+"\n")
        return observations

    def printForPopulationExcel(self,wb, row, observations):
        toPrint = []
        for parameter in self.parameters:
            toPrint.append("%f "%parameter)
        observations = np.vstack([observations, self.parameters])
        excelWriteRow(toPrint+[self.R2,self.R2adj,self.AIC,self.AICc,self.BIC], wb, row)
        return observations, row+1

    def getBasicInfo(self):
        """ Return a string with some basic information of the fitting. """
        info = ""
        info += "Sample name: %s\n" % self.sampleName
        info += "Model: %s\n" % self.modelEquation
        info += "R2: %f\n" % self.R2
        info += "R2adj: %f\n" % self.R2adj
        info += "AIC: %f\n" % self.AIC
        info += "AICc(Recommended): %f\n" % self.AICc
        info += "BIC: %f\n" % self.BIC

        return info

    def _printToStream(self, fh):
        if self.significance == None:
            self.significance = ["Undetermined"]*len(self.parameters)
        fh.write(self.getBasicInfo())
        fh.write("Parameter lowerBound upperBound IsStatisticallySignificant -------\n")
        for parameter, lower, upper, significance in izip(self.parameters,self.lowerBound,self.upperBound,\
                                                          self.significance):
            fh.write("%f [%s,%s] %s\n"%(parameter,str(lower),str(upper),significance))
        fh.write("X   Y   Ypredicted [Ylower,Yupper] -------\n")
        for j in range(len(self.x)):
            fh.write("Series %d -----\n"%j)
            xj = self.x[j]
            yj = self.y[j]
            ypj = self.yp[j]
            ylj = self.yl[j]
            yuj = self.yu[j]
            for x,y,yp,yl,yu in izip(xj,yj,ypj,ylj,yuj):
                fh.write("%f %s %s [%s,%s]\n"%(x,str(y),str(yp),str(yl),str(yu)))
        fh.write("\n")

    def _printToExcel(self, wb, row, parameterNames):
        if self.significance == None:
            self.significance = ["Undetermined"]*len(self.parameters)
        basicInfo = self.getBasicInfo()
        i=0
        for line in basicInfo.split('\n'):
            if i == 0:
                excelWriteRow(line.split(':'), wb, row, bold=True)
                excelFillCells(wb, row)
            else:
                excelWriteRow(line.split(':'), wb, row)
            row += 1
            i+=1
        excelWriteRow(["Parameter","Value","lowerBound","upperBound","IsStatisticallySignificant"],wb,row); row+=1
        for parameterName, parameter, lower, upper, significance in izip(parameterNames,self.parameters,
                                                          self.lowerBound,self.upperBound,self.significance):
            excelWriteRow([parameterName, parameter,str(lower),str(upper),significance], wb, row); row+=1
        row+=1
        excelWriteRow(["X","Y","Ypredicted", "[Ylower","Yupper]"],wb,row); row+=1

        for j in range(len(self.x)):
            excelWriteRow("Series %d"%j, wb, row); row+=1
            xj = self.x[j]
            yj = self.y[j]
            ypj = self.yp[j]
            ylj = self.yl[j]
            yuj = self.yu[j]
            for x,y,yp,yl,yu in izip(xj,yj,ypj,ylj,yuj):
                try:
                    ypToPrint=float(yp)
                except:
                    ypToPrint=str(yp)
                try:
                    ylToPrint=float(yl)
                except:
                    ylToPrint=str(yl)
                try:
                    yuToPrint=float(yu)
                except:
                    yuToPrint=str(yu)
                excelWriteRow([float(x), float(y), ypToPrint, ylToPrint, yuToPrint], wb, row); row+=1
        return row+1

    def restartReadingState(self):
        self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_NAME

    def readFromLine(self, line):
        if self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_NAME:
            tokens = line.split(':')
            self.sampleName = tokens[1].strip()
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_MODELEQ

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_MODELEQ:
            tokens = line.split(':')
            self.modelEquation = tokens[1].strip()
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_R2

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_R2:
            tokens = line.split(':')
            self.R2 = float(tokens[1])
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_R2ADJ

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_R2ADJ:
            tokens = line.split(':')
            self.R2adj = float(tokens[1])
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_AIC

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_AIC:
            tokens = line.split(':')
            self.AIC = float(tokens[1])
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_AICc

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_AICc:
            tokens = line.split(':')
            self.AICc = float(tokens[1])
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_BIC

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_BIC:
            tokens = line.split(':')
            self.BIC = float(tokens[1])
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_PARAMETER_BOUNDS

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_PARAMETER_BOUNDS:
            if line.startswith("Parameter lowerBound upperBound"):
                self.parameters=[]
                self.lowerBound=[]
                self.upperBound=[]
                self.significance=[]
            elif line.startswith("X   Y   Ypredicted"):
                self.state=PKPDSampleFit.READING_SAMPLEFITTINGS_SAMPLE_VALUES
                self.x=[]
                self.y=[]
                self.yp=[]
                self.yl=[]
                self.yu=[]
            else:
                tokens=line.split()
                self.parameters.append(float(tokens[0]))
                self.significance.append(tokens[2])
                tokens=(tokens[1])[1:-1].split(',')
                self.lowerBound.append(tokens[0])
                self.upperBound.append(tokens[1])

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_SAMPLE_VALUES:
            tokens=line.split()
            if tokens[0].strip()=="Series":
                self.x.append([])
                self.y.append([])
                self.yp.append([])
                self.yl.append([])
                self.yu.append([])
                self.multiOutputSeries = True
            else:
                if self.multiOutputSeries:
                    self.x[-1].append(float(tokens[0]))
                    self.y[-1].append(float(tokens[1]))
                    self.yp[-1].append(float(tokens[2]))
                    tokens=(tokens[3])[1:-1].split(',')
                    self.yl[-1].append(tokens[0])
                    self.yu[-1].append(tokens[1])
                else:
                    self.x.append(float(tokens[0]))
                    self.y.append(float(tokens[1]))
                    self.yp.append(float(tokens[2]))
                    tokens=(tokens[3])[1:-1].split(',')
                    self.yl.append(tokens[0])
                    self.yu.append(tokens[1])

    def copyFromOptimizer(self,optimizer):
        self.R2 = optimizer.R2
        self.R2adj = optimizer.R2adj
        self.AIC = optimizer.AIC
        self.AICc = optimizer.AICc
        self.BIC = optimizer.BIC
        self.significance = optimizer.significance
        self.lowerBound = optimizer.lowerBound
        self.upperBound = optimizer.upperBound


class PKPDSampleFitBootstrap:
    READING_SAMPLEFITTINGS_NAME = 0
    READING_SAMPLEFITTINGS_XB = 1
    READING_SAMPLEFITTINGS_YB = 2
    READING_SAMPLEFITTINGS_PARAMETERS = 3

    def __init__(self):
        self.sampleName = ""
        self.R2 = []
        self.R2adj = []
        self.AIC = []
        self.AICc = []
        self.BIC = []
        self.parameters = None
        self.xB = []
        self.yB = []

    def _printSample(self,fh,n, n0=0):
        outputStr = "%d: "%(n+n0)
        for parameter in self.parameters[n,:]:
            outputStr += "%f "%parameter
        outputStr += " # %f %f %f %f %f"%(self.R2[n],self.R2adj[n],self.AIC[n],self.AICc[n],self.BIC[n])
        fh.write(outputStr+"\n")

    def _printSampleExcel(self,wb, row, n, n0=0):
        toPrint = ["Sample %d"%(n+n0)]
        for parameter in self.parameters[n,:]:
            toPrint.append(parameter)
        excelWriteRow(toPrint+[self.R2[n],self.R2adj[n],self.AIC[n],self.AICc[n],self.BIC[n]],wb,row)
        return row+1

    def printForPopulation(self,fh,observations):
        for n in range(0,self.parameters.shape[0]):
            self._printSample(fh,n,observations.shape[0])
        observations = np.vstack([observations, self.parameters])
        return observations

    def printForPopulationExcel(self,wb,row,observations):
        for n in range(0,self.parameters.shape[0]):
            row=self._printSampleExcel(wb,row,n,observations.shape[0])
        observations = np.vstack([observations, self.parameters])
        return observations,row+1

    def _printToStream(self,fh):
        fh.write("Sample name: %s\n"%self.sampleName)
        for n in range(0,self.parameters.shape[0]):
            fh.write("xB: %s\n"%self.xB[n])
            fh.write("yB: %s\n"%self.yB[n])

            outputStr = ""
            for parameter in self.parameters[n,:]:
                outputStr += "%f "%parameter
            outputStr += " # %f %f %f %f %f"%(self.R2[n],self.R2adj[n],self.AIC[n],self.AICc[n],self.BIC[n])
            fh.write(outputStr+"\n")
        fh.write("\n")

    def _printToExcel(self,wb,row,parameterNames):
        excelWriteRow(["Sample name:",self.sampleName],wb,row,bold=True); row+=1
        for n in range(0,self.parameters.shape[0]):
            excelWriteRow("Sample %d"%n, wb, row); row+=1
            excelWriteRow(["xB:"]+self.xB[n][1:-1].split(), wb, row); row+=1
            excelWriteRow(["yB:"]+self.yB[n][1:-1].split(), wb, row); row+=1
        return row+1

    def restartReadingState(self):
        self.state = PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_NAME

    def readFromLine(self, line):
        if self.state==PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_NAME:
            tokens = line.split(':')
            self.sampleName = tokens[1].strip()
            self.strRead = ""
            self.state = PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_XB

        elif self.state==PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_XB:
            if ":" in line:
                tokens = line.split(':')
                self.strRead += tokens[1].strip()
            else:
                self.strRead += line
            if "]" in self.strRead:
                self.xB.append(self.strRead.strip())
                self.strRead=""
                self.state = PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_YB

        elif self.state==PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_YB:
            if ":" in line:
                tokens = line.split(':')
                self.strRead += tokens[1].strip()
            else:
                self.strRead += line
            if "]" in self.strRead:
                self.yB.append(self.strRead.strip())
                self.strRead=""
                self.state = PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_PARAMETERS

        elif self.state==PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_PARAMETERS:
            tokens = line.split('#')
            tokensParameters = tokens[0].strip().split(' ')
            tokensQuality = tokens[1].strip().split(' ')
            if self.parameters is None:
                self.parameters = np.empty((0,len(tokensParameters)),np.double)
            self.parameters = np.vstack([self.parameters, [float(prm) for prm in tokensParameters]])

            self.R2.append(float(tokensQuality[0]))
            self.R2adj.append(float(tokensQuality[1]))
            self.AIC.append(float(tokensQuality[2]))
            self.AICc.append(float(tokensQuality[3]))
            self.BIC.append(float(tokensQuality[4]))

            self.state = PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_XB

    def copyFromOptimizer(self,optimizer):
        self.R2.append(optimizer.R2)
        self.R2adj.append(optimizer.R2adj)
        self.AIC.append(optimizer.AIC)
        self.AICc.append(optimizer.AICc)
        self.BIC.append(optimizer.BIC)


class PKPDFitting(EMObject):
    READING_FITTING_EXPERIMENT = 1
    READING_FITTING_PREDICTOR = 2
    READING_FITTING_PREDICTED = 3
    READING_FITTING_PREDICTED_LIST = 4
    READING_FITTING_MODEL = 5
    READING_POPULATION_HEADER = 6
    READING_POPULATION = 7
    READING_SAMPLEFITTINGS_BEGIN = 8
    READING_SAMPLEFITTINGS_CONTINUE = 9

    def __init__(self, cls="", **args):
        EMObject.__init__(self, **args)
        self.fnFitting = String()
        self.fnExperiment = String()
        self.predictor = None
        self.predicted = None
        self.modelDescription = ""
        self.modelParameters = []
        self.modelParameterUnits = []
        self.sampleFits = []
        self.summaryLines = []
        if cls=="":
            self.sampleFittingClass = "PKPDSampleFit"
        else:
            self.sampleFittingClass = cls

    def isPopulation(self):
        if self.fnFitting.get() is None:
            return False
        return self.fnFitting.get().endswith("bootstrapPopulation.pkpd")

    def write(self, fnFitting, writeToExcel=True):
        fh=open(fnFitting,'w')
        self._printToStream(fh)
        fh.close()
        self.fnFitting.set(fnFitting)
        writeMD5(fnFitting)

        if writeToExcel:
            self.writeToExcel(os.path.splitext(fnFitting)[0] + ".xlsx")

    def getAllParameters(self):
        allParameters = np.empty((0,len(self.modelParameters)),np.double)
        for sampleFitting in self.sampleFits:
            allParameters = np.vstack([allParameters, sampleFitting.parameters])
        return allParameters

    def getStats(self, observations=None):
        if observations is None:
            observations = self.getAllParameters()
        mu=np.mean(observations,axis=0)
        C=np.cov(np.transpose(observations))
        sigma = np.sqrt(np.diag(C))
        R=np.corrcoef(np.transpose(observations))
        percentiles = np.percentile(observations,[0, 2.5, 25, 50, 75, 97.5, 100],axis=0)
        return mu, sigma, R, percentiles

    def _printToStream(self,fh):
        fh.write("[FITTING] ===========================\n")
        fh.write("Experiment: %s\n"%self.fnExperiment.get())
        fh.write("Predictor (X): ")
        self.predictor._printToStream(fh)
        fh.write("Predicted (Y): ")
        if type(self.predicted)==list:
            fh.write("Predicted list=%d\n"%len(self.predicted))
            for y in self.predicted:
                y._printToStream(fh)
        else:
            self.predicted._printToStream(fh)
        fh.write("Model: %s\n"%self.modelDescription)
        fh.write("\n")

        fh.write("[POPULATION PARAMETERS] =============\n")
        auxUnit = PKPDUnit()
        for paramName, paramUnits in izip(self.modelParameters, self.modelParameterUnits):
            auxUnit.unit = paramUnits
            fh.write("%s [%s] "%(paramName,auxUnit._toString()))
        fh.write(" # R2 R2adj AIC AICc BIC\n")
        observations = np.empty((0,len(self.modelParameters)),np.double)
        for sampleFitting in self.sampleFits:
            observations = sampleFitting.printForPopulation(fh,observations)
        fh.write("\n")

        mu=np.mean(observations,axis=0)
        if observations.shape[0]>2:
            C=np.cov(np.transpose(observations))
            if not C.shape:
                R=np.asarray([1.0])
                sigma = C
            else:
                R = np.corrcoef(np.transpose(observations))
                sigma = np.sqrt(np.diag(C))
            fh.write("Mean   parameters  = %s\n"%np.array_str(mu))
            fh.write("Median parameters  = %s\n"%np.array_str(np.median(observations,axis=0)))
            limits = np.percentile(observations,[2.5,97.5],axis=0)
            fh.write("Lower bound (2.5%%) = %s\n"%np.array_str(limits[0]))
            fh.write("Upper bound (97.5%%) = %s\n"%np.array_str(limits[1]))
            fh.write("Covariance matrix  =\n%s\n"%np.array_str(C,max_line_width=120))
            fh.write("Correlation matrix  =\n%s\n"%np.array_str(R,max_line_width=120))
        fh.write("\n")

        fh.write("[SAMPLE FITTINGS] ===================\n")
        for sampleFitting in self.sampleFits:
            sampleFitting._printToStream(fh)

    def writeToExcel(self,fnXls):
        wb = openpyxl.Workbook()
        wb.active.title = "Experiment"

        currentRow = 1
        excelWriteRow("FITTING",wb,currentRow,bold=True); excelFillCells(wb,currentRow); currentRow+=1
        excelWriteRow(["Experiment:",self.fnExperiment.get()], wb, currentRow); currentRow+=1
        excelWriteRow("Predictor (X): ", wb, currentRow);
        currentRow=self.predictor._printToExcel(wb,currentRow,2)
        toPrint=["Predicted (Y): "]
        if type(self.predicted)==list:
            excelWriteRow(toPrint+["Predicted list=%d"%len(self.predicted)], wb, currentRow); currentRow+=1
            for y in self.predicted:
                currentRow=y._printToExcel(wb,currentRow)
        else:
            excelWriteRow(toPrint, wb, currentRow)
            currentRow=self.predicted._printToExcel(wb,currentRow,2); currentRow+=1
        excelWriteRow(["Model:",self.modelDescription], wb, currentRow); currentRow+=1
        currentRow+=1

        excelWriteRow("POPULATION PARAMETERS",wb,currentRow,bold=True); excelFillCells(wb,currentRow); currentRow+=1
        auxUnit = PKPDUnit()
        toPrint=[]
        for paramName, paramUnits in izip(self.modelParameters, self.modelParameterUnits):
            auxUnit.unit = paramUnits
            toPrint.append("%s [%s] "%(paramName,auxUnit._toString()))
        excelWriteRow(toPrint+["R2","R2adj","AIC","AICc","BIC"],wb,currentRow); currentRow+=1
        observations = np.empty((0,len(self.modelParameters)),np.double)
        for sampleFitting in self.sampleFits:
            observations, currentRow = sampleFitting.printForPopulationExcel(wb,currentRow,observations)
        currentRow+=1

        mu=np.mean(observations,axis=0)
        if observations.shape[0]>2:
            excelWriteRow(["Mean   parameters  =",np.array_str(mu)], wb, currentRow); currentRow+=1
            excelWriteRow(["Median parameters  =",np.array_str(np.median(observations,axis=0))], wb, currentRow); currentRow+=1
            limits = np.percentile(observations,[2.5,97.5],axis=0)
            excelWriteRow(["Lower bound (2.5%%) =",np.array_str(limits[0])], wb, currentRow); currentRow+=1
            excelWriteRow(["Upper bound (97.5%%) =",np.array_str(limits[1])], wb, currentRow); currentRow+=1
        currentRow+=1

        excelWriteRow("SAMPLE FITTINGS",wb,currentRow,bold=True); excelFillCells(wb,currentRow); currentRow+=1
        for sampleFitting in self.sampleFits:
            currentRow=sampleFitting._printToExcel(wb,currentRow,self.modelParameters)

        excelAdjustColumnWidths(wb)
        wb.save(fnXls)

    def load(self, fnFitting=None):
        fnFitting = str(fnFitting or self.fnFitting)
        fh = open(fnFitting)
        if not fh:
            raise Exception("Cannot open %s" % fnFitting)
        if not verifyMD5(fnFitting):
            raise Exception("The file %s has been modified since its creation" % fnFitting)
        self.fnFitting.set(fnFitting)

        auxUnit = PKPDUnit()
        for line in fh.readlines():
            line=line.strip()
            if line=="":
                if state==PKPDFitting.READING_SAMPLEFITTINGS_CONTINUE:
                     state=PKPDFitting.READING_SAMPLEFITTINGS_BEGIN
                continue
            if line.startswith('[') and line.endswith('='):
                section = line.split('=')[0].strip().lower()
                if section=="[fitting]":
                    state=PKPDFitting.READING_FITTING_EXPERIMENT
                    self.summaryLines.append(line)
                elif section=="[population parameters]":
                    state=PKPDFitting.READING_POPULATION_HEADER
                    self.summaryLines.append(line)
                elif section=="[sample fittings]":
                    state=PKPDFitting.READING_SAMPLEFITTINGS_BEGIN
                else:
                    print("Skipping: ",line)

            elif state==PKPDFitting.READING_FITTING_EXPERIMENT:
                tokens = line.split(':')
                self.fnExperiment.set(tokens[1].strip())
                state = PKPDFitting.READING_FITTING_PREDICTOR
                self.summaryLines.append(line)

            elif state==PKPDFitting.READING_FITTING_PREDICTOR:
                tokens = line.split(':')
                self.predictor = PKPDVariable()
                self.predictor.parseTokens(tokens[1].split(';'))
                state = PKPDFitting.READING_FITTING_PREDICTED
                self.summaryLines.append(line)

            elif state==PKPDFitting.READING_FITTING_PREDICTED:
                tokens = line.split(':')
                if (tokens[1].strip().startswith("Predicted list")):
                    state = PKPDFitting.READING_FITTING_PREDICTED_LIST
                    self.remainingPredicted = int(tokens[1].split('=')[1])
                    self.predicted = []
                    self.summaryLines.append(line)
                else:
                    self.predicted = PKPDVariable()
                    self.predicted.parseTokens(tokens[1].split(';'))
                    state = PKPDFitting.READING_FITTING_MODEL
                    self.summaryLines.append(line)

            elif state==PKPDFitting.READING_FITTING_PREDICTED_LIST:
                    self.summaryLines.append(line)
                    newVar = PKPDVariable()
                    newVar.parseTokens(line.split(';'))
                    self.predicted.append(newVar)
                    self.remainingPredicted -= 1
                    if self.remainingPredicted == 0:
                        state = PKPDFitting.READING_FITTING_MODEL

            elif state==PKPDFitting.READING_FITTING_MODEL:
                tokens = line.split(':')
                self.modelDescription = tokens[1].strip()
                state = PKPDFitting.READING_POPULATION
                self.summaryLines.append(line)
                self.summaryLines.append("\n")

            elif state==PKPDFitting.READING_POPULATION_HEADER:
                lineParts = line.split('#')
                tokens = lineParts[0].strip().split(' ')
                for i in range(0,len(tokens),2):
                    self.modelParameters.append(tokens[i])
                    self.modelParameterUnits.append(auxUnit._fromString(tokens[i+1][1:-1]))
                state = PKPDFitting.READING_POPULATION

            elif state==PKPDFitting.READING_POPULATION:
                self.summaryLines.append(line)

            elif state==PKPDFitting.READING_SAMPLEFITTINGS_BEGIN:
                newSampleFit = eval("%s()"%self.sampleFittingClass)
                self.sampleFits.append(newSampleFit)
                self.sampleFits[-1].restartReadingState()
                self.sampleFits[-1].readFromLine(line)
                state = PKPDFitting.READING_SAMPLEFITTINGS_CONTINUE

            elif state==PKPDFitting.READING_SAMPLEFITTINGS_CONTINUE:
                self.sampleFits[-1].readFromLine(line)

        fh.close()

    def getSampleFit(self, sampleName):
        for sampleFit in self.sampleFits:
            if sampleFit.sampleName == sampleName:
                return sampleFit
        return None

    def loadExperiment(self):
        experiment = PKPDExperiment()
        experiment.load(self.fnExperiment.get())
        return experiment


class PKPDSampleSignalAnalysis:
    def __init__(self):
        self.sampleName = ""
        self.x = None
        self.y = None
        self.analysisVariables = []
        self.parameters = None
        self.lowerBound = None
        self.upperBound = None
        self.significance = None

    def _printToStream(self,fh):
        if self.significance == None:
            self.significance = ["Undetermined"]*len(self.parameters)
        fh.write("Sample name: %s\n"%self.sampleName)
        for varName, parameter in izip(self.analysisVariables,self.parameters):
            fh.write("%s = %f\n"%(varName,parameter))
        fh.write("\n")


class PKPDSignalAnalysis(EMObject):
    READING_EXPERIMENT = 1
    READING_PREDICTOR = 2
    READING_PREDICTED = 3
    READING_MODEL = 4
    READING_ANALYSIS_VARIABLES = 5
    READING_POPULATION_HEADER = 6
    READING_POPULATION = 7
    READING_SAMPLEANALYSIS_NAME = 8
    READING_SAMPLEANALYSIS_PARAMETERS = 9

    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.fnAnalysis = String()
        self.fnExperiment = String()
        self.predictor = None
        self.predicted = None
        self.analysisDescription = ""
        self.analysisParameters = []
        self.sampleAnalyses = []
        self.summaryLines = []

    def write(self, fnAnalysis):
        fh=open(fnAnalysis,'w')
        self._printToStream(fh)
        fh.close()
        self.fnAnalysis.set(fnAnalysis)
        writeMD5(fnAnalysis)

    def _printToStream(self,fh):
        fh.write("[SIGNAL ANALYSIS] ===========================\n")
        fh.write("Experiment: %s\n"%self.fnExperiment.get())
        fh.write("Predictor (X): ")
        self.predictor._printToStream(fh)
        fh.write("Predicted (Y): ")
        self.predicted._printToStream(fh)
        fh.write("Analysis: %s\n"%self.analysisDescription)
        fh.write("\n")

        fh.write("[POPULATION PARAMETERS] =============\n")
        fh.write(' '.join(self.analysisParameters)+"\n")
        i=0
        for sampleAnalysis in self.sampleAnalyses:
            outputStr = ""
            j=0
            for parameter in sampleAnalysis.parameters:
                outputStr += "%f "%parameter
                if i==0 and j==0:
                    observations = np.zeros([len(self.sampleAnalyses),len(sampleAnalysis.parameters)])
                observations[i,j]=parameter
                j+=1
            fh.write(outputStr+"\n")
            i+=1
        fh.write("\n")

        mu=np.mean(observations,axis=0)
        C=np.cov(np.transpose(observations))
        sigma = np.sqrt(np.diag(C))
        R=np.corrcoef(np.transpose(observations))
        fh.write("Mean   parameters  = %s\n"%np.array_str(mu))
        fh.write("Median parameters  = %s\n"%np.array_str(np.median(observations,axis=0)))
        fh.write("Lower bound (95%%, independent Gaussians) = %s\n"%np.array_str(mu-1.96*sigma))
        fh.write("Upper bound (95%%, independent Gaussians) = %s\n"%np.array_str(mu+1.96*sigma))
        fh.write("Covariance matrix  =\n%s\n"%np.array_str(C,max_line_width=120))
        fh.write("Correlation matrix  =\n%s\n"%np.array_str(R,max_line_width=120))
        fh.write("\n")

        fh.write("[SAMPLE ANALYSES] ===================\n")
        for sampleAnalysis in self.sampleAnalyses:
            sampleAnalysis._printToStream(fh)

    def load(self,fnAnalysis):
        fh=open(fnAnalysis)
        if not fh:
            raise Exception("Cannot open %s"%fnAnalysis)
        if not verifyMD5(fnAnalysis):
            raise Exception("The file %s has been modified since its creation"%fnAnalysis)
        self.fnAnalysis.set(fnAnalysis)

        for line in fh.readlines():
            line=line.strip()
            if line=="":
                if state==PKPDSignalAnalysis.READING_SAMPLEANALYSIS_PARAMETERS:
                     state=PKPDSignalAnalysis.READING_SAMPLEANALYSIS_NAME
                continue
            if line.startswith('[') and line.endswith('='):
                section = line.split('=')[0].strip().lower()
                if section=="[signal analysis]":
                    state=PKPDSignalAnalysis.READING_EXPERIMENT
                    self.summaryLines.append(line)
                elif section=="[population parameters]":
                    state=PKPDSignalAnalysis.READING_POPULATION_HEADER
                    self.summaryLines.append(line)
                elif section=="[sample analyses]":
                    state=PKPDSignalAnalysis.READING_SAMPLEANALYSIS_NAME
                else:
                    print("Skipping: ",line)

            elif state==PKPDSignalAnalysis.READING_EXPERIMENT:
                tokens = line.split(':')
                self.fnExperiment.set(tokens[1].strip())
                state = PKPDSignalAnalysis.READING_PREDICTOR
                self.summaryLines.append(line)

            elif state==PKPDSignalAnalysis.READING_PREDICTOR:
                tokens = line.split(':')
                self.predictor = PKPDVariable()
                self.predictor.parseTokens(tokens[1].split(';'))
                state = PKPDSignalAnalysis.READING_PREDICTED
                self.summaryLines.append(line)

            elif state==PKPDSignalAnalysis.READING_PREDICTED:
                tokens = line.split(':')
                self.predicted = PKPDVariable()
                self.predicted.parseTokens(tokens[1].split(';'))
                state = PKPDSignalAnalysis.READING_MODEL
                self.summaryLines.append(line)

            elif state==PKPDSignalAnalysis.READING_MODEL:
                tokens = line.split(':')
                self.analysisDescription = tokens[1].strip()
                state = PKPDSignalAnalysis.READING_POPULATION
                self.summaryLines.append(line)
                self.summaryLines.append("\n")

            elif state==PKPDSignalAnalysis.READING_POPULATION_HEADER:
                self.analysisParameters=line.split(' ')
                state = PKPDSignalAnalysis.READING_POPULATION

            elif state==PKPDSignalAnalysis.READING_POPULATION:
                self.summaryLines.append(line)

            elif state==PKPDSignalAnalysis.READING_SAMPLEANALYSIS_NAME:
                tokens = line.split(':')
                self.sampleAnalyses.append(PKPDSampleSignalAnalysis())
                self.sampleAnalyses[-1].sampleName = tokens[1].strip()
                self.sampleAnalyses[-1].parameters=[]
                state = PKPDSignalAnalysis.READING_SAMPLEANALYSIS_PARAMETERS

            elif state==PKPDSignalAnalysis.READING_SAMPLEANALYSIS_PARAMETERS:
                tokens=line.split('=')
                self.sampleAnalyses[-1].analysisVariables.append(tokens[0].strip())
                self.sampleAnalyses[-1].parameters.append(float(tokens[1].strip()))

        fh.close()

    def getSampleAnalysis(self,sampleName):
        for sampleAnalysis in self.sampleAnalyses:
            if sampleAnalysis.sampleName == sampleName:
                return sampleAnalysis
        return None


class PKPDAllometricScale(EMObject):
    READING_PREDICTOR = 1
    READING_MODELS = 2
    READING_X = 3
    READING_Y = 4

    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.fnScale = String()
        self.predictor = ""
        self.predictorUnits = ""
        self.scaled_vars = []
        self.averaged_vars = []
        self.models = {}
        self.qualifiers = {}
        self.confidence = 0
        self.X = None
        self.Y = None

    def write(self, fnScale):
        fh=open(fnScale,'w')
        self._printToStream(fh)
        fh.close()
        self.fnScale.set(fnScale)
        writeMD5(fnScale)

    def _printToStream(self,fh):
        fh.write("[ALLOMETRIC SCALING] ===========================\n")
        fh.write("Predictor (X): %s %s\n"%(self.predictor,self.predictorUnits))
        for varName, varUnits in self.scaled_vars:
            model = self.models[varName]
            qualifiers = self.qualifiers[varName]
            # fh.write("Scale model: %s=(%f)*%s^(%f) %f%% Confidence interval=[%f,%f]"%\
            #          (varName,model[0],model[1],self.confidence,qualifiers[0],qualifiers[1]))
            fh.write("Scale model: %s=(%f)*%s^(%f); R2=%f; %f%% Confidence intervals (y=k*x^a) k=[%f,%f] a=[%f,%f]; %s\n" % \
                     (varName, model[0], self.predictor, model[1], qualifiers[0], self.confidence,
                      qualifiers[1], qualifiers[2], qualifiers[3], qualifiers[4], varUnits))
        for varName, varUnits in self.averaged_vars:
            model = self.models[varName]
            qualifiers = self.qualifiers[varName]
            fh.write("Average model: %s=%f Std=%f; %s\n"%(varName,model[0],qualifiers[0],varUnits))

        fh.write("\n")
        fh.write("[DATA] =========================================\n")
        fh.write("%s= %s\n"%(self.predictor,str(self.X)))
        for varName, y in self.Y.iteritems():
            fh.write("%s= %s\n"%(varName,str(y)))

    def load(self,fnScale):
        fh=open(fnScale)
        if not fh:
            raise Exception("Cannot open %s"%fnScale)
        if not verifyMD5(fnScale):
            raise Exception("The file %s has been modified since its creation"%fnScale)
        self.fnScale.set(fnScale)

        for line in fh.readlines():
            line=line.strip()
            if line=="":
                continue
            if line.startswith('[') and line.endswith('='):
                section = line.split('=')[0].strip().lower()
                if section=="[allometric scaling]":
                    state=PKPDAllometricScale.READING_PREDICTOR
                elif section=="[data]":
                    state=PKPDAllometricScale.READING_X
                else:
                    print("Skipping: ",line)

            elif state==PKPDAllometricScale.READING_PREDICTOR:
                tokens = line.split(':')[1]
                tokens = tokens.split()
                self.predictor = tokens[0].strip()
                self.predictorUnits = tokens[1].strip()
                state = PKPDAllometricScale.READING_MODELS

            elif state==PKPDAllometricScale.READING_MODELS:
                tokens = line.split(':')
                if tokens[0]=="Scale model":
                    tokens = tokens[1].split(';')
                    # "%s=(%f)*%s^(%f); R2=%f; %f%% Confidence intervals (y=k*x^a) k=[%f,%f] a=[%f,%f]"
                    idxL=0
                    idxR=tokens[0].find('=')
                    varName=tokens[0][idxL:idxR].strip()
                    idxL=idxR+2
                    idxR=tokens[0].find(')',idxL)
                    k=float(tokens[0][idxL:idxR])
                    idxL=tokens[0].find('(',idxR)
                    idxR=tokens[0].find(')',idxL)
                    a=float(tokens[0][idxL+1:idxR])
                    self.models[varName]=[k, a]

                    idxL=tokens[1].find('=')
                    R2=float(tokens[1][idxL+1:])

                    idxR=tokens[2].find('%')
                    self.confidence=float(tokens[2][:idxR])

                    idxL=tokens[2].find('[',idxR+1)
                    idxR=tokens[2].find(',',idxL+1)
                    kl=float(tokens[2][idxL+1:idxR])
                    idxL=idxR
                    idxR=tokens[2].find(']',idxL+1)
                    ku=float(tokens[2][idxL+1:idxR])

                    idxL=tokens[2].find('[',idxR+1)
                    idxR=tokens[2].find(',',idxL+1)
                    al=float(tokens[2][idxL+1:idxR])
                    idxL=idxR
                    idxR=tokens[2].find(']',idxL+1)
                    au=float(tokens[2][idxL+1:idxR])

                    self.qualifiers[varName]=[R2,kl,ku,al,au]

                    varUnits = tokens[3].strip()
                    self.scaled_vars.append((varName,varUnits))

                elif tokens[0]=="Average model":
                    # Average model: %s=%f Std=%f
                    idxR=tokens[1].find('=')
                    varName=tokens[1][:idxR].strip()
                    idxL=idxR
                    idxR=tokens[1].find(' ',idxL)
                    mean=float(tokens[1][idxL+1:idxR])
                    idxL=tokens[1].find('=',idxR+1)
                    idxR=tokens[1].find(';',idxL+1)
                    std=float(tokens[1][idxL+1:idxR])
                    varUnits = tokens[1][idxR+1:]
                    self.models[varName]=[mean]
                    self.qualifiers[varName]=[std]
                    self.averaged_vars.append((varName,varUnits))
                else:
                    print("Skipping: ",line)

            elif state==PKPDAllometricScale.READING_X:
                tokens = line.split('=')
                self.X = []
                for x in tokens[1].strip()[1:-1].split(','):
                    self.X.append(float(x))
                self.Y = {}
                state = PKPDAllometricScale.READING_Y

            elif state==PKPDAllometricScale.READING_Y:
                tokens = line.split('=')
                varName = tokens[0].strip()
                self.Y[varName] = []
                for y in tokens[1].strip()[1:-1].split(','):
                    self.Y[varName].append(float(y))

        fh.close()


class PKPDDoseResponse(EMObject):
    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.doses=[]
        self.Nresponses=[]
        self.Npatients=[]
        self.responses=[]

    def read(self,inputStr,isFn=False):
        if isFn:
            fh=open(inputStr)
            inputStr=fh.readlines()
            fh.close()
        else:
            inputStr=inputStr.split("\n")
        for line in inputStr:
            line=line.strip()
            if line!="":
                tokens=line.split(":")
                dose=float(tokens[0].strip())
                self.doses.append(dose)
                self.responses.append(tokens[1].strip())

                k=0
                n=0
                for token in tokens[1].split():
                    n+=1
                    if token.strip()=="1":
                        k+=1
                self.Nresponses.append(k)
                self.Npatients.append(n)

    def write(self,fnOut):
        fh=open(fnOut,"w")
        for i in range(len(self.doses)):
            fh.write("%f: %s\n"%(self.doses[i],self.responses[i]))
        fh.close()

    def findDose(self,dose):
        # for n in range(len(self.doses)):
        #     if abs(self.doses[n]-dose)<1e-6:
        #         return n
        if len(self.doses)>0 and abs(self.doses[-1]-dose)<1e-6:
            return len(self.doses)-1
        else:
            return None

    def appendResponse(self,dose,response):
        # response must be Boolean
        n = self.findDose(dose)
        if n is None:
            n=len(self.doses)
            self.doses.append(dose)
            self.responses.append("")
            self.Nresponses.append(0)
            self.Npatients.append(0)

        responseStr=" 1" if response else " 0"
        self.responses[n]=(self.responses[n]+responseStr).strip()
        self.Npatients[n]+=1
        if response:
            self.Nresponses[n]+=1


def flattenArray(y):
    if type(y[0])!=list and type(y[0])!=np.ndarray:
        y = [np.array(y,dtype=np.float32)]
    else:
        y = [np.array(yi,dtype=np.float32) for yi in y]
    return y

def smartLog(y):
    return np.array([math.log10(yi) if np.isfinite(yi) and yi>0 else float("inf") for yi in y])


class PKPDDataSet:
    _datasetDict = {}  # store all created datasets

    def __init__(self, name, folder, files, url=None):
        """
        Params:

        #filesDict is dict with key, value pairs for each file
        """
        self._datasetDict[name] = self
        self.folder = folder
        import pkg_resources
        from pyworkflow.install.plugin_funcs import PluginInfo
        package = PluginInfo('scipion-pkpd', 'scipion-pkpd',
                             remote=False).pipName
        dist = pkg_resources.get_distribution(package).location
        self.path = join(dist, 'pkpd',
                         'data', 'test', folder)
        self.filesDict = files
        self.url = url

    def getFile(self, key):
        if key in self.filesDict:
            return join(self.path, self.filesDict[key])
        return join(self.path, key)

    def getPath(self):
        return self.path

    @classmethod
    def getDataSet(cls, name):
        """
        This method is called every time the dataset want to be retrieved
        """
        assert name in cls._datasetDict, "Dataset: %s dataset doesn't exist." % name

        ds = cls._datasetDict[name]
        folder = ds.folder
        url = '' if ds.url is None else ' -u ' + ds.url

        if not pwutils.envVarOn('SCIPION_TEST_NOSYNC'):
            command = ("%s %s testdata --download %s %s"
                       % (pw.PYTHON, pw.getScipionScript(), folder, url))
            print(">>>> %s" % command)
            os.system(command)
        return cls._datasetDict[name]
