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
from .protocol_pkpd_sa_base import ProtPKPDSABase
from pkpd.models.sa_models import NCAExpIVModel
from pkpd.pkpd_units import PKPDUnit

# TESTED in test_workflow_gabrielsson_pk07.py

class ProtPKPDNCAIVExp(ProtPKPDSABase):
    """ Non-compartmental analysis based on an exponential fitting.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'nca iv exponentials'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('protExponential', params.PointerParam, label="Input exponential fitting",
                      pointerClass='ProtPKPDExponentialFit',
                      help='The input experiment will be taken from the exponential fitting')
        form.addParam('protElimination', params.PointerParam, label="Elimination rate",
                      pointerClass='ProtPKPDEliminationRate',
                      help='Select an execution of a protocol estimating the elimination rate')
        form.addParam("absorptionF", params.FloatParam, label="Absorption fraction", default=1,
                      help="Between 0 (=no absorption) and 1 (=full absorption)")

    def getInputExperiment(self):
        return self.protExponential.get().inputExperiment.get()

    def getListOfFormDependencies(self):
        return [self.protElimination.get().getObjId(), self.protExponential.get().getObjId()]

    #--------------------------- STEPS functions --------------------------------------------
    def setupFromFormParameters(self):
        self.exponentialFitting = self.readFitting(self.protExponential.get().outputFitting.fnFitting.get())
        self.eliminationFitting = self.readFitting(self.protElimination.get().outputFitting.fnFitting.get())

    def getXYvars(self):
        self.varNameX = self.protElimination.get().predictor.get()
        self.varNameY = self.protElimination.get().predicted.get()

    def createAnalysis(self):
        self.analysis = NCAExpIVModel()
        self.analysis.setExperiment(self.experiment)
        self.analysis.setXVar(self.varNameX)
        self.analysis.setYVar(self.varNameY)
        self.analysis.F = self.absorptionF.get()
        self.analysis.CnUnits = PKPDUnit()
        self.analysis.CnUnits.unit = self.exponentialFitting.modelParameterUnits[0]
        self.analysis.lambdanUnits = PKPDUnit()
        self.analysis.lambdanUnits.unit = self.exponentialFitting.modelParameterUnits[1]

    def prepareForSampleAnalysis(self, sampleName):
        sample = self.experiment.samples[sampleName]

        sampleFit = self.eliminationFitting.getSampleFit(sampleName)
        sample.interpretDose()
        self.analysis.D = sample.getDoseAt(0.0)
        if sampleFit == None:
            print("  Cannot process %s because its elimination rate cannot be found\n\n"%sampleName)
            return False
        self.analysis.lambdaz = sampleFit.parameters[1]
        self.analysis.lambdazUnits = PKPDUnit()
        self.analysis.lambdazUnits.unit = self.eliminationFitting.modelParameterUnits[1]
        print("Elimination rate = %f [%s]"%(self.analysis.lambdaz,self.analysis.lambdazUnits._toString()))

        sampleFit = self.exponentialFitting.getSampleFit(sampleName)
        self.analysis.Cn = []
        self.analysis.lambdan = []
        for i in range(0,len(sampleFit.parameters),2):
            self.analysis.Cn.append(sampleFit.parameters[i])
            self.analysis.lambdan.append(sampleFit.parameters[i+1])
            print("C%d = %f [%s]"%(i/2,self.analysis.Cn[-1],self.analysis.CnUnits._toString()))
            print("lambda%d = %f [%s]"%(i/2,self.analysis.lambdan[-1],self.analysis.lambdanUnits._toString()))
        print(' ')

        return True

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        msg.append("Non-compartmental analysis for the observations of the variable %s"%self.protElimination.get().predicted.get())
        return msg

    def _validate(self):
        msg = []
        if self.protExponential.get().predictor.get()!=self.protElimination.get().predictor.get():
            msg.append("The predictor of the exponential fitting (%s) does not match the predictor of the elimination protocol (%s)"%\
                       (self.protExponential.get().predictor.get(),self.protElimination.get().predictor.get()))
        if self.protExponential.get().predicted.get()!=self.protElimination.get().predicted.get():
            msg.append("The predicted of the exponential fitting (%s) does not match the predicted of the elimination protocol (%s)"%\
                       (self.protExponential.get().predictor.get(),self.protElimination.get().predictor.get()))
        return msg

    def _warnings(self):
        msg = []
        experiment = self.readExperiment(self.getInputExperiment().fnPKPD,show=False)
        incorrectList = experiment.getNonBolusDoses()
        if len(incorrectList)!=0:
            msg.append("This protocol is meant only for intravenous bolus regimens. Check the doses for %s"%(','.join(incorrectList)))
        return msg

