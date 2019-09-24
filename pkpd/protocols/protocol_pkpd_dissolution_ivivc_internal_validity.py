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
import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD

# Tested in test_workflow_levyplot.py
# Tested in test_workflow_deconvolution2.py

class ProtPKPDIVIVCInternalValidity(ProtPKPD):
    """ This protocol compares the AUC and Cmax predicted from an in vitro-in vivo experiment and
        the AUC and Cmax of the in vivo population from which the IVIV correlation was estimated.
        There should not be differences larger than 15%% between the two sets of variables.
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'ivivc internal validity'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="True in-vivo experiment",
                      pointerClass='PKPDExperiment',
                      help='Select the experiment with the measurements you want to analyze. It must have Cmax and AUC0t')
        form.addParam('inputSimulated', params.PointerParam, label="Simulated in-vivo experiment",
                      pointerClass='PKPDExperiment',
                      help='Select the experiment with the measurements you want to analyze. It must have Cmax and AUC0t'
                            'and it is the output of a IVIVC+PK simulation')

    #--------------------------- STEPS functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runAnalysis',self.inputExperiment.getObjId(),self.inputSimulated.getObjId())

    def runAnalysis(self,objId1,objId2):
        trueExp = self.readExperiment(self.inputExperiment.get().fnPKPD,False)
        simExp = self.readExperiment(self.inputSimulated.get().fnPKPD,False)

        if not "AUC0t" in trueExp.variables or not "Cmax" in trueExp.variables:
            raise Exception("Cannot find AUC0t or Cmax in the in vivo experiment")
        if not "AUC0t" in simExp.variables or not "Cmax" in simExp.variables or not "from" in simExp.variables:
            raise Exception("Cannot find AUC0t, Cmax or from in the simulated experiment")

        trueAUC = {}
        trueCmax = {}
        for sampleName, sample in trueExp.samples.iteritems():
            try:
                trueAUC[sampleName]=float(sample.descriptors["AUC0t"])
                trueCmax[sampleName]=float(sample.descriptors["Cmax"])
            except:
                pass

        simAUC = {}
        simCmax = {}
        simNames = []
        for sampleName, sample in simExp.samples.iteritems():
            fromNames=sample.descriptors["from"]
            sampleName=fromNames.split("---")[0]
            if not sampleName in simNames:
                simNames.append(sampleName)
            if not sampleName in simAUC:
                simAUC[sampleName]=[]
                simCmax[sampleName]=[]
            try:
                simAUC[sampleName].append(float(sample.descriptors["AUC0t"]))
                simCmax[sampleName].append(float(sample.descriptors["Cmax"]))
            except:
                pass

        errorAUC=[]
        errorCmax=[]
        for sampleName in trueAUC:
            tAUC = trueAUC[sampleName]
            tCmax = trueCmax[sampleName]
            if sampleName in simAUC:
                if tAUC>0:
                    for sAUC in simAUC[sampleName]:
                        errorAUC.append((tAUC-sAUC)/tAUC*100)
                if tCmax>0:
                    for sCmax in simCmax[sampleName]:
                        errorCmax.append((tCmax-sCmax)/tCmax*100)

        np.savetxt(self._getExtraPath("errorAUC.txt"),errorAUC)
        np.savetxt(self._getExtraPath("errorCmax.txt"),errorCmax)

        if len(errorAUC)==0:
            print("Cannot find any matching name between the true experiment (%s) and the simulated names (%s)"%\
                     (" ".join(trueExp.samples.keys())," ".join(simNames)))
        else:
            alpha_2 = (100-95)/2
            limits = np.percentile(errorAUC,[alpha_2,100-alpha_2])
            fhSummary=open(self._getPath("summary.txt"),"w")
            self.doublePrint(fhSummary,"error AUC (normalized to 100, (true-sim)/true) %f%% confidence interval=[%f,%f] mean=%f"%(95,limits[0],limits[1],np.mean(errorAUC)))
            limits = np.percentile(errorCmax,[alpha_2,100-alpha_2])
            self.doublePrint(fhSummary,"error Cmax (normalized to 100, (true-sim)/true) %f%% confidence interval=[%f,%f] mean=%f"%(95,limits[0],limits[1],np.mean(errorCmax)))
            fhSummary.close()

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        self.addFileContentToMessage(msg,self._getPath("summary.txt"))
        return msg

