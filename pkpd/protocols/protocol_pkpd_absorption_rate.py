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
from .protocol_pkpd_fit_base import ProtPKPDFitBase
from pkpd.models.pk_models import PKPDSimpleEVModel
from pkpd.pkpd_units import strUnit
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.object import Integer
import math
import numpy as np

# TESTED in test_workflow_gabrielsson_pk02.py
# TESTED in test_workflow_gabrielsson_pk06.py

class ProtPKPDAbsorptionRate(ProtPKPDFitBase):
    """ Estimation of the absorption rate for a non-intravenous route. The estimation is performed after estimating
        the elimination rate. The experiment is determined by the\n
        Protocol created by http://www.kinestatpharma.com\n.
        See the theory at http://www.pharmpress.com/files/docs/Basic%20Pharmacokinetics%20sample.pdf"""
    _label = 'absorption rate'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('protElimination', params.PointerParam, label="Elimination rate",
                      pointerClass='ProtPKPDEliminationRate',
                      help='Select an execution of a protocol estimating the elimination rate')
        form.addParam("absorptionF", params.FloatParam, label="Absorption fraction", default=1,
                      help="Between 0 (=no absorption) and 1 (=full absorption)")

        form.addParam('bounds', params.StringParam, label="Ka, V, [tlag] bounds", default="", expertLevel=LEVEL_ADVANCED,
                      help='Bounds for Ka (absorption constant), V (distribution volume) and optionally tlag.\nExample 1: (0,1e-3);(30,50);(0.1,0.5) -> Ka in (0,1e-3), V in (30,50) and tlag in (0.1,0.5)\n')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')
        form.addParam('includeTlag', params.BooleanParam, label="Include tlag", default=True, expertLevel=LEVEL_ADVANCED,
                      help='Calculate the delay between administration and absorption')
        self.fitType=Integer() # Logarithmic fit
        self.fitType.set(1)

    def getListOfFormDependencies(self):
        return [self.protElimination.get().getObjId(), self.absorptionF.get(), self.bounds.get(), self.confidenceInterval.get()]

    #--------------------------- STEPS functions --------------------------------------------
    def setupFromFormParameters(self):
        self.eliminationFitting = self.readFitting(self.protElimination.get().outputFitting.fnFitting.get())
        self.model.F = self.absorptionF.get()

    def getXYvars(self):
        self.varNameX = self.protElimination.get().predictor.get()
        self.varNameY = self.protElimination.get().predicted.get()

    def createModel(self):
        return PKPDSimpleEVModel(self.includeTlag.get())

    def prepareForSampleAnalysis(self, sampleName):
        sample = self.experiment.samples[sampleName]
        sampleFit = self.eliminationFitting.getSampleFit(sampleName)
        sample.interpretDose()
        self.model.D = sample.getDoseAt(0.0)
        self.model.Dunits = sample.getDoseUnits()
        if sampleFit == None:
            print("  Cannot process %s because its elimination rate cannot be found\n\n"%sampleName)
            return False
        self.model.C0 = sampleFit.parameters[0]
        self.model.C0units = self.eliminationFitting.modelParameterUnits[0]
        self.model.Ke = sampleFit.parameters[1]
        self.model.KeUnits = self.eliminationFitting.modelParameterUnits[1]
        print("Concentration at t=0 = %f [%s]"%(self.model.C0,strUnit(self.model.C0units)))
        print("Elimination rate = %f [%s]"%(self.model.Ke,strUnit(self.model.KeUnits)))

        self.experiment.addParameterToSample(sampleName, "Ke", self.model.KeUnits, "Automatically estimated elimination rate", self.model.Ke)

        return True

    def postSampleAnalysis(self, sampleName):
        xunits = self.experiment.getVarUnits(self.varNameX)
        Cunits = self.experiment.getVarUnits(self.varNameY)
        Ka = self.model.parameters[0]
        Ke = self.model.Ke
        tmax = math.log(Ka/Ke)/(Ka-Ke)

        self.experiment.addParameterToSample(sampleName, "tmax", xunits, "Estimated time of the Maximum of the non-iv peak", tmax)
        Cmax = self.model.forwardModel(self.model.parameters,x=[np.atleast_1d(np.array(tmax))])[0][0]
        self.experiment.addParameterToSample(sampleName, "Cmax", Cunits, "Estimated concentration of the Maximum of the non-iv peak", Cmax)
        print("tmax = %f [%s]"%(tmax,strUnit(xunits)))
        print("Cmax = %f [%s]"%(Cmax,strUnit(Cunits)))

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        msg.append("Non-compartmental analysis for the observations of the variable %s"%self.protElimination.get().predicted.get())
        return msg

    def _warnings(self):
        msg = []
        experiment = self.readExperiment(self.getInputExperiment().fnPKPD,show=False)
        incorrectList = experiment.getNonBolusDoses()
        if len(incorrectList)!=0:
            msg.append("This protocol is meant only for bolus regimens. Check the doses for %s"%(','.join(incorrectList)))

        return msg

