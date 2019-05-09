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
import math
import random

import pyworkflow.protocol.params as params
from pkpd.objects import PKPDExperiment, PKPDSample, PKPDVariable
from .protocol_pkpd import ProtPKPD


# tested in ***

class ProtPKPDDissolutionF2(ProtPKPD):
    """ Calculate the f1 and f2 from two dissolution profiles."""

    _label = 'dissol f1 and f2'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputRef', params.PointerParam, label="Reference dissolution profiles",
                      pointerClass='PKPDExperiment', help='Select an experiment with dissolution profiles')
        form.addParam('inputTest', params.PointerParam, label="Test dissolution profiles",
                      pointerClass='PKPDExperiment', help='Select an experiment with dissolution profiles')
        form.addParam('dissolutionVar', params.StringParam, label="Dissolution variable", default="C",
                      help='Which variable contains the profile. It must be the same in both experiments.')
        form.addParam('Nboostrap', params.IntParam, label="Number of bootstrap samples", default=200,
                      help='Per pair of profiles, set to 0 for no bootstrapping')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllF',self.inputRef.get().getObjId(),self.inputTest.get().getObjId())

    #--------------------------- STEPS functions --------------------------------------------
    def getProfiles(self,prmExp,varName):
        experiment = self.readExperiment(prmExp.fnPKPD)
        allY = []
        for sampleName, sample in experiment.samples.iteritems():
            y=sample.getValues(varName)
            allY.append(np.asarray(y,dtype=np.float64))
        return allY

    def calculateF(self,pRef,pTest):
        diff = pRef-pTest
        D2 = (np.square(diff)).mean(axis=None)
        f2=50*math.log(100.0/math.sqrt(1+D2),10.0)
        f1= np.sum(np.abs(diff))/np.sum(pRef)*100
        return f1, f2

    def printStats(self,allF,Fstr,Fformula):
        mu=np.mean(allF)
        sigma = np.std(allF)
        percentiles = np.percentile(allF,[0, 2.5, 25, 50, 75, 97.5, 100])
        retval=""
        retval +="%s = %s\n"%(Fstr,Fformula)
        retval +="%s distribution with B=%d bootstrap samples (total of %d samples)\n"%(Fstr,self.Nboostrap.get(),len(allF))
        retval +="%s mean+-std: %f+-%f\n"%(Fstr,mu,sigma)
        retval +="%s minimum,maximum: [%f,%f]\n"%(Fstr,percentiles[0],percentiles[6])
        retval +="%s percentile [2.5,97.5]%%: [%f,%f]\n"%(Fstr,percentiles[1],percentiles[5])
        retval +="%s percentile [25,75]%%: [%f,%f]\n"%(Fstr,percentiles[2],percentiles[4])
        retval +="%s percentile 50%%: %f\n"%(Fstr,percentiles[3])
        return retval

    def calculateAllF(self, objId1, objId2):
        profilesRef=self.getProfiles(self.inputRef.get(),self.dissolutionVar.get())
        profilesTest=self.getProfiles(self.inputTest.get(),self.dissolutionVar.get())

        allF1=[]
        allF2=[]
        if len(profilesRef)>0:
            idx = [k for k in range(0, len(profilesRef[0]))]
            Nidx = len(idx)
        for profileRef in profilesRef:
            for profileTest in profilesTest:
                f1,f2=self.calculateF(profileRef,profileTest)
                allF1.append(f1)
                allF2.append(f2)
                if self.Nboostrap.get()>0:
                    for n in range(self.Nboostrap.get()):
                        idxB = sorted(np.random.choice(idx,Nidx))
                        profileRefB = np.asarray([profileRef[i] for i in idxB])
                        profileTestB = np.asarray([profileTest[i] for i in idxB])
                        f1, f2 = self.calculateF(profileRefB, profileTestB)
                        allF1.append(f1)
                        allF2.append(f2)
        strF1=self.printStats(allF1,"F1","sum(|pRef-pTest|)/sum(pRef)*100")
        strF2=self.printStats(allF2,"F2","50*log10(100/sqrt(1+mean(|pRef-pTest|^2)))")
        np.savetxt(self._getExtraPath("f1.txt"),allF1)
        np.savetxt(self._getExtraPath("f2.txt"),allF2)

        self.printSection("Results")
        fhSummary = open(self._getPath("summary.txt"),"w")
        self.doublePrint(fhSummary,strF2)
        self.doublePrint(fhSummary,"---------------------------")
        self.doublePrint(fhSummary,strF1)
        fhSummary.close()

    def _validate(self):
        return []

    def _summary(self):
        retval = []
        self.addFileContentToMessage(retval,self._getPath("summary.txt"))
        return retval

    def _citations(self):
        retval = ['Islam2018']
        return retval