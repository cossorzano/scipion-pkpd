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
from pkpd.models.sa_models import NCAEVModel
from pkpd.pkpd_units import PKPDUnit, strUnit
from pyworkflow.protocol.constants import LEVEL_ADVANCED

# TESTED in test_workflow_gabrielsson_pk06.py

class ProtPKPDNCAEV(ProtPKPDSABase):
    """ Non-compartmental analysis of a non-intravenous bolus.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'nca ev'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('protAbsorption', params.PointerParam, label="Non-IV Input experiment",
                      pointerClass='ProtPKPDAbsorptionRate', important=True,
                      help='Select an execution of a protocol estimating the absorption rate.')
        form.addParam('areaCalc', params.EnumParam, choices=["Trapezoidal","Mixed"],
                      label="Method for AUC, AUMC calculation", default=1, expertLevel=LEVEL_ADVANCED,
                      help='The mixed integration uses the trapezoidal is used in the raising side and the log-trapezoidal in the decay side.\n'
                            'See explanation at http://learnpkpd.com/2011/04/02/calculating-auc-linear-and-log-linear\n')
        form.addParam('protNCAIV', params.PointerParam, label="IV Input experiment (optional)", allowsNull=True,
                      pointerClass='ProtPKPDNCAIVExp,ProtPKPDNCAIVObs',
                      help='A companion experiment with intraveneous bolus is used to estimate the bioavailability and the Mean Absorption Time. '
                           'The sample names in both experiments must be exactly the same.')

    def getListOfFormDependencies(self):
        return [self.protAbsorption.get().getObjId()]

    #--------------------------- STEPS functions --------------------------------------------
    def getInputExperiment(self):
        return self.protAbsorption.get().outputExperiment

    def setupFromFormParameters(self):
        self.experimentIV = None
        if self.protNCAIV.get() != None:
            self.experimentIV=self.readExperiment(self.protNCAIV.get().outputExperiment.fnPKPD)

    def getXYvars(self):
        self.varNameX = self.protAbsorption.get().protElimination.get().predictor.get()
        self.varNameY = self.protAbsorption.get().protElimination.get().predicted.get()

    def createAnalysis(self):
        self.analysis = NCAEVModel()
        self.analysis.setExperiment(self.experiment)
        self.analysis.setXVar(self.varNameX)
        self.analysis.setYVar(self.varNameY)
        self.analysis.F = self.protAbsorption.get().absorptionF.get()
        if self.areaCalc == 0:
            self.analysis.areaCalc = "Trapezoidal"
        elif self.areaCalc == 1:
            self.analysis.areaCalc = "Mixed"

    def prepareForSampleAnalysis(self, sampleName):
        sample = self.experiment.samples[sampleName]
        sample.interpretDose()
        self.analysis.D = sample.getDoseAt(0.0)
        self.analysis.Ke = float(sample.descriptors["Ke"])
        return True

    def postAnalysis(self, sampleName):
        if self.experimentIV!=None:
            if sampleName in self.experimentIV.samples:
                sampleNIV = self.experiment.samples[sampleName]
                DNIV= self.analysis.D

                sampleIV = self.experimentIV.samples[sampleName]
                sampleIV.interpretDose()
                DIV = sampleIV.getDoseAt(0.0)

                F = float(sampleNIV.descriptors['AUC_0inf'])/float(sampleNIV.descriptors['AUC_0inf'])*DIV/DNIV
                MAT = float(sampleNIV.descriptors['MRT']) - float(sampleIV.descriptors['MRT'])
                self.experiment.addParameterToSample(sampleName, "F", PKPDUnit.UNIT_NONE, "Estimated bioavailability", F)
                self.experiment.addParameterToSample(sampleName, "MAT", self.analysis.parameterUnits[-1], "Estimated Mean Absorption Time", MAT)
                print("F = %f"%F)
                print("MAT = %f [%s]"%(MAT,strUnit(self.analysis.parameterUnits[-1])))

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        msg.append("Non-compartmental analysis for the observations of the variable %s"%self.protAbsorption.get().protElimination.get().predicted.get())
        return msg

    def _warnings(self):
        experiment = self.readExperiment(self.getInputExperiment().fnPKPD,show=False)
        incorrectList = experiment.getNonBolusDoses()
        if len(incorrectList)==0:
            return []
        else:
            return ["This protocol is meant only for intravenous bolus regimens. Check the doses for %s"%(','.join(incorrectList))]

