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
from os.path import exists

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.pkpd_units import PKPDUnit

# TESTED in test_workflow_gabrielsson_pk06.py

class ProtPKPDNCAEstimateBioavailability(ProtPKPD):
    """ Estimate bioavailability as F=(AUCpo/Dpo) / (AUCiv/Div) [Gabrielsson 2010, p. 546], i.e., the ratio \n
        between the oral and intravenous AUC normalized by their respective doses.
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'bioavailability nca'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('protNCAEV', params.PointerParam, label="Extra-vascular NCA",
                      pointerClass='ProtPKPDNCAEV',
                      help='Non-compartmental analysis of the Extra-Vascular route, the dose must be a bolus\n'
                           'The sample names in both (EV and IV) analysis must be exactly the same.')
        form.addParam('protNCAIV', params.PointerParam, label="IV NCA",
                      pointerClass='ProtPKPDNCAIVExp,ProtPKPDNCAIVObs',
                      help='Non-compartmental analysis of the Extra-Vascular route, the dose must be a bolus\n'
                           'The sample names in both (EV and IV) analysis must be exactly the same.')

    #--------------------------- STEPS functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runEstimateF')
        self._insertFunctionStep('createOutput')

    def runEstimateF(self):
        self.experimentEV=self.readExperiment(self.protNCAEV.get().outputExperiment.fnPKPD)
        self.experimentIV=self.readExperiment(self.protNCAIV.get().outputExperiment.fnPKPD)

        Flist = []
        for sampleEVname, sampleEV in self.experimentEV.samples.iteritems():
            if sampleEVname in self.experimentIV.samples:
                sampleEV.interpretDose()
                DEV=sampleEV.getDoseAt(0)

                sampleIV=self.experimentIV.samples[sampleEVname]
                sampleIV.interpretDose()
                DIV=sampleIV.getDoseAt(0)

                AUCev=sampleEV.descriptors['AUC_0inf']
                AUCiv=sampleIV.descriptors['AUC_0inf']
                if AUCev!="NA" and AUCiv!="NA":
                    F=(float(AUCev)/DEV)/(float(AUCiv)/DIV)
                    Flist.append(F)
                    self.experimentEV.addParameterToSample(sampleEV.sampleName,"bioavailability",PKPDUnit.UNIT_NONE,
                                                           "bioavailability",F)
                else:
                    self.experimentEV.addParameterToSample(sampleEV.sampleName,"bioavailability",PKPDUnit.UNIT_NONE,
                                                           "bioavailability","NA")
            else:
                self.experimentEV.addParameterToSample(sampleEV.sampleName,"bioavailability",PKPDUnit.UNIT_NONE,
                                                       "bioavailability","NA")
        self.experimentEV.write(self._getPath("experiment.pkpd"))
        fhSummary=open(self._getPath("summary.txt"),"w")
        if len(Flist)>0:
            Fmean = np.mean(Flist)
            fhSummary.write("Average bioavailability (normalized between 0 and 1)=%f\n"%Fmean)
            if Fmean>1:
                fhSummary.write("Bioavailability larger than 1 is caused by PK non-linearities and dependences of AUC on the dose\n")
        fhSummary.close()


    def createOutput(self):
        self._defineOutputs(outputExperiment=self.experimentEV)
        self._defineSourceRelation(self.protNCAEV, self.experimentEV)
        self._defineSourceRelation(self.protNCAIV, self.experimentEV)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        self.addFileContentToMessage(msg,self._getPath("summary.txt"))
        return msg

    def _citations(self):
        return ['Gabrielsson2010']
