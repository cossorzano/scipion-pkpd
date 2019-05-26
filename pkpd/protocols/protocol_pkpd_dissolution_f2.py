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
from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD


# tested in test_workflow_dissolution_f2.py

class ProtPKPDDissolutionF2(ProtPKPD):
    """ Calculate the f1 and f2 from two dissolution profiles."""

    _label = 'dissol f1 and f2'

    BYVECTOR = 0
    BYPOINT = 1

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputRef', params.PointerParam, label="Reference dissolution profiles",
                      pointerClass='PKPDExperiment', help='Select an experiment with dissolution profiles')
        form.addParam('inputTest', params.PointerParam, label="Test dissolution profiles",
                      pointerClass='PKPDExperiment', help='Select an experiment with dissolution profiles')
        form.addParam('timeVar', params.StringParam, label="Time variable", default="t",
                      help='Which variable contains the time points. It must be the same in both experiments.')
        form.addParam('dissolutionVar', params.StringParam, label="Dissolution variable", default="C",
                      help='Which variable contains the profile. It must be the same in both experiments.')
        form.addParam('Nbootstrap', params.IntParam, label="Number of bootstrap samples", default=200,
                      help='Per pair of profiles, set to 0 for no bootstrapping')
        form.addParam('bootstrapBy', params.EnumParam, label="Bootrstap by", choices=['Vessel','Time point'], default=0,
                      help='Bootstrapping per vessel will take all the samples from the same vessel (some time points may be repeated). '
                           'The total number of samples is Nref*Ntest*Nbootstrap where Nref is the number of reference vessels, Ntest the number of test vessels, '
                           'and Nbootstrap the number of bootstrap samples.'
                           'Bootstrapping per time point will mix all vessels so that at each time point a random sample from any of the vessels is taken.'
                           'The total number of samples is Nbootstrap samples')
        form.addParam('resampleT', params.FloatParam, label="Resample profiles (time step)", default=-1,
                      help='Resample the input profiles at this time step (make sure it is in the same units as the input). '
                           'Leave it to -1 for no resampling')
        form.addParam('confidence', params.FloatParam, label="Confidence (%%)", default=95,
                      help='Confidence level')


    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllF',self.inputRef.get().getObjId(),self.inputTest.get().getObjId())

    #--------------------------- STEPS functions --------------------------------------------
    def getProfiles(self,prmExp,varNameT,varNameC):
        experiment = self.readExperiment(prmExp.fnPKPD)
        allY = []
        for sampleName, sample in experiment.samples.iteritems():
            x=np.asarray(sample.getValues(varNameT),dtype=np.float64)
            y=np.asarray(sample.getValues(varNameC),dtype=np.float64)
            if self.resampleT.get()>0:
                B = InterpolatedUnivariateSpline(x, y, k=1)
                y = B(np.arange(np.min(x),np.max(x),self.resampleT.get()))
            allY.append(y)
        return allY

    def calculateF(self,pRef,pTest):
        idxm85=np.argwhere(np.logical_and(pRef<=85,pTest<=85)).tolist()
        idxm85 = [item for sublist in idxm85 for item in sublist]
        idxp85=np.argwhere(np.logical_or(pRef>85,pTest>85)).tolist()
        idxp85 = [item for sublist in idxp85 for item in sublist]
        idx=idxm85
        if len(idxp85)>0:
            idxp85=np.random.choice(idxp85,1)
            idx.append(idxp85[0])
        idx=sorted(idx)

        diff = pRef[idx]-pTest[idx]
        D2 = (np.square(diff)).mean(axis=None)
        f2=50*math.log(100.0/math.sqrt(1+D2),10.0)
        f1= np.sum(np.abs(diff))/np.sum(pRef[idx])*100

        print("Bootstrap sample %d" % self.b)
        print("Reference measures: %s"%np.array2string(pRef[idx],max_line_width=10000))
        print("Test measures: %s"%np.array2string(pTest[idx],max_line_width=10000))
        print("f1=%f f2=%f"%(f1,f2))
        print(" ")
        self.b = self.b + 1

        return f1, f2

    def printStats(self,allF,Fstr,Fformula):
        allF=[f for f in allF if not np.isnan(f)]
        mu=np.mean(allF)
        sigma = np.std(allF)
        alpha=1-self.confidence.get()/100.0
        percentiles = np.percentile(allF,[0, alpha/2*100, 25, 50, 75, (1-alpha/2)*100, 100])
        retval=""
        retval +="%s = %s\n"%(Fstr,Fformula)
        retval +="%s distribution with B=%d bootstrap samples (total of %d samples)\n"%(Fstr,self.Nbootstrap.get(),len(allF))
        retval +="%s mean+-std: %f+-%f\n"%(Fstr,mu,sigma)
        retval +="%s minimum,maximum: [%f,%f]\n"%(Fstr,percentiles[0],percentiles[6])
        retval +="%s percentile [%f,%f]%%: [%f,%f]\n"%(Fstr,alpha/2*100,(1-alpha/2)*100,percentiles[1],percentiles[5])
        retval +="%s percentile [25,75]%%: [%f,%f]\n"%(Fstr,percentiles[2],percentiles[4])
        retval +="%s percentile 50%%: %f\n"%(Fstr,percentiles[3])
        return retval

    def bootstrapByTimePoint(self,profilesList):
        Nsamples = len(profilesList)
        Ntimepoints = profilesList[0].size
        bootstrapSample = np.zeros((Ntimepoints))
        for j in range(Ntimepoints):
            i=np.random.randint(0,Nsamples)
            bootstrapSample[j]=profilesList[i][j]
        return bootstrapSample

    def calculateAllF(self, objId1, objId2):
        profilesRef=self.getProfiles(self.inputRef.get(), self.timeVar.get(), self.dissolutionVar.get())
        profilesTest=self.getProfiles(self.inputTest.get(),self.timeVar.get(), self.dissolutionVar.get())

        allF1=[]
        allF2=[]
        if len(profilesRef)>0:
            idx = [k for k in range(0, len(profilesRef[0]))]
            Nidx = len(idx)
        self.printSection("Calculations")
        self.b=1
        if self.bootstrapBy.get()==ProtPKPDDissolutionF2.BYVECTOR:
            for profileRef in profilesRef:
                print("Full reference profile: %s"%np.array2string(profileRef,max_line_width=10000))
                for profileTest in profilesTest:
                    print("Full test profile: %s" % np.array2string(profileTest, max_line_width=10000))
                    print(" ")
                    f1,f2=self.calculateF(profileRef,profileTest)
                    allF1.append(f1)
                    allF2.append(f2)
                    if self.Nbootstrap.get()>0:
                        for n in range(self.Nbootstrap.get()):
                            idxB = sorted(np.random.choice(idx,Nidx))
                            profileRefB = np.asarray([profileRef[i] for i in idxB])
                            profileTestB = np.asarray([profileTest[i] for i in idxB])
                            f1, f2 = self.calculateF(profileRefB, profileTestB)
                            allF1.append(f1)
                            allF2.append(f2)
        else:
            for n in range(self.Nbootstrap.get()):
                profileRefB = self.bootstrapByTimePoint(profilesRef)
                profileTestB = self.bootstrapByTimePoint(profilesTest)
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
