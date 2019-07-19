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
from scipy.stats import norm

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.utils import uniqueFloatValues


# tested in test_workflow_dissolution_f2.py

class ProtPKPDDissolutionF2(ProtPKPD):
    """ Calculate the f1 and f2 from two dissolution profiles. The bootstrap confidence interval is bias corrected and accelarated."""

    _label = 'dissol f1 and f2'

    BYVECTOR = 0
    BYPOINT = 1
    BYPOINTAVG = 2

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
        form.addParam('f2mode', params.EnumParam, label="F2 mode", choices=['EMA','FDA'], default=0,
                      help='The EMA allows only 1 sample when any of the reference or test goes above 85%%, '
                      'but the FDA allows 1 sample when both of them go above 85%%')
        form.addParam('Nbootstrap', params.IntParam, label="Number of bootstrap samples", default=200,
                      help='Per pair of profiles, set to 0 for no bootstrapping')
        form.addParam('bootstrapBy', params.EnumParam, label="Bootrstap by", choices=['Vessel','Time point','Time point and average'], default=2,
                      help='Bootstrapping per vessel will take all the samples from the same vessel (some time points may be repeated). '
                           'Bootstrapping per time point will mix the vessels to form a time profile with as many samples as the input ones. '
                           'These two techniques (bootstrapping per vessel or time point) calculates F1 and F2 at the level of single dissolution profile. '
                           'Opposed to this, bootstrapping by time point and average produces a whole experiment with as many vessels as the input, where '
                           'the profiles have been shuffled at the level of time. Then, the average profile of both bootstrap experiments is performed and '
                           'the F1 and F2 of the two averages are calculated. '
                           'The total number of samples is Nref*Ntest*Nbootstrap where Nref is the number of reference vessels, Ntest the number of test vessels, '
                           'and Nbootstrap the number of bootstrap samples.')
        form.addParam('resampleT', params.FloatParam, label="Resample profiles (time step)", default=-1,
                      help='Resample the input profiles at this time step (make sure it is in the same units as the input). '
                           'Leave it to -1 for no resampling')
        form.addParam('keepResample', params.BooleanParam, label="Generate output with resampled profiles", default=False,
                      condition='resampleT>0',
                      help='Create an output experiment with the resampled profiles')
        form.addParam('confidence', params.FloatParam, label="Confidence (%%)", default=95,
                      help='Confidence level')


    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateAllF',self.inputRef.get().getObjId(),self.inputTest.get().getObjId())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def getProfiles(self,prmExp,varNameT,varNameC):
        experiment = self.readExperiment(prmExp.fnPKPD)
        allY = []
        for sampleName, sample in experiment.samples.iteritems():
            x=np.asarray(sample.getValues(varNameT),dtype=np.float64)
            y=np.asarray(sample.getValues(varNameC),dtype=np.float64)
            if self.resampleT.get()>0:
                x, y = uniqueFloatValues(x, y)
                B = InterpolatedUnivariateSpline(x, y, k=1)
                xp = np.arange(np.min(x),np.max(x)+self.resampleT.get(),self.resampleT.get())
                y = B(xp)
                if self.keepResample.get():
                    sample.setValues(varNameT, [str(xi) for xi in xp.tolist()])
                    sample.setValues(varNameC, [str(yi) for yi in y.tolist()])
            allY.append(y)
        if self.resampleT.get()>0 and self.keepResample.get():
            self.experiment = experiment
            self.experiment.write(self._getPath("experiment.pkpd"))
        return allY

    def randomIdx(self,pRef,pTest):
        if self.f2mode.get()==0: # EMA, only 1 sample if any of the two is above 85 %
            idxm85=np.argwhere(np.logical_and(pRef<=85,pTest<=85)).tolist()
            idxp85=np.argwhere(np.logical_or(pRef>85,pTest>85)).tolist()
        else: #FDA, only 1 sample if both of them are above 85%
            idxp85=np.argwhere(np.logical_and(pRef>85,pTest>85)).tolist()
            idxm85=np.argwhere(np.logical_or(pRef<=85,pTest<=85)).tolist()
        idxm85 = [item for sublist in idxm85 for item in sublist]
        idxp85 = [item for sublist in idxp85 for item in sublist]
        idx = idxm85
        if len(idxp85) > 0:
            idxp85 = np.random.choice(idxp85, 1)
            idx.append(idxp85[0])
        return sorted(idx)

    def calculateF(self,pRef,pTest, randomizeIdx):
        if randomizeIdx:
            idx=self.randomIdx(pRef,pTest)
        else:
            idx=np.arange(0,len(pRef))
        counter=0
        while len(idx)<3 or len(set(idx))<2:
            if randomizeIdx:
                idx = self.randomIdx(pRef, pTest)
            else:
                idx = np.arange(0, len(pRef))
            counter+=1
            if counter>20:
                return np.nan,np.nan

        diff = pRef[idx]-pTest[idx]
        D2 = (np.square(diff)).mean(axis=None)
        f2=50*math.log(100.0/math.sqrt(1+D2),10.0)
        f1= np.sum(np.abs(diff))/np.sum(pRef[idx])*100

        # print("Reference measures: %s"%np.array2string(pRef[idx],max_line_width=10000))
        # print("Test measures: %s"%np.array2string(pTest[idx],max_line_width=10000))
        # print("f1=%f f2=%f"%(f1,f2))
        # print(" ")

        return f1, f2

    def printStats(self,allF,Fstr,Fformula,alphaL,alphaU):
        allF=[f for f in allF if not np.isnan(f)]
        mu=np.mean(allF)
        sigma = np.std(allF)
        percentiles = np.percentile(allF,[0, alphaL*100, 25, 50, 75, alphaU*100, 100])
        alpha=1-self.confidence.get()/100
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

    def bootstrapByTimePointAvg(self,profilesList):
        Nsamples = len(profilesList)
        Ntimepoints = profilesList[0].size
        randomExperiment = np.zeros((Nsamples,Ntimepoints))
        for i in range(Nsamples):
            for j in range(Ntimepoints):
                iFrom=np.random.randint(0,Nsamples)
                randomExperiment[i,j]=profilesList[iFrom][j]
        return np.mean(randomExperiment,axis=0)

    def calculateAllF(self, objId1, objId2):
        profilesRef=self.getProfiles(self.inputRef.get(), self.timeVar.get(), self.dissolutionVar.get())
        profilesTest=self.getProfiles(self.inputTest.get(),self.timeVar.get(), self.dissolutionVar.get())

        # Sample estimate ----------------
        self.printSection("Sample estimate")
        allF10=[]
        allF20=[]
        for profileRef in profilesRef:
            print("Full reference profile: %s" % np.array2string(profileRef, max_line_width=10000))
            for profileTest in profilesTest:
                print("Full test profile: %s" % np.array2string(profileTest, max_line_width=10000))
                print(" ")
                print("calculateAllF pRef",profileRef)
                print("calculateAllF pTest",profileTest)
                f1, f2 = self.calculateF(profileRef, profileTest, False)
                allF10.append(f1)
                allF20.append(f2)
#        f10,f20 = self.calculateF(np.mean(np.asarray(profilesRef),axis=0),np.mean(np.asarray(profilesTest),axis=0), False)
        f10=np.mean([f for f in allF10 if not np.isnan(f)])
        f20=np.mean([f for f in allF20 if not np.isnan(f)])

        # Jack knife ----------------
        self.printSection("Jackknife")
        allF1J=[]
        allF2J=[]
        for profileRef in profilesRef:
            print("Full reference profile: %s" % np.array2string(profileRef, max_line_width=10000))
            for profileTest in profilesTest:
                print("Full test profile: %s" % np.array2string(profileTest, max_line_width=10000))
                print(" ")
                for i in range(profileRef.size):
                    idx=[j for j in range(profileRef.size) if j != i]
                    f1, f2 = self.calculateF(profileRef[idx], profileTest[idx], False)
                    allF1J.append(f1)
                    allF2J.append(f2)
        f1J=np.mean([f for f in allF1J if not np.isnan(f)])
        f2J=np.mean([f for f in allF2J if not np.isnan(f)])
        allF1J=[f1 for f1 in allF1J if not np.isnan(f1)]
        allF2J=[f2 for f2 in allF2J if not np.isnan(f2)]

        # Bootstrapping ------------------
        self.printSection("Bootstrapping")
        allF1b=[]
        allF2b=[]
        if len(profilesRef)>0:
            idx = [k for k in range(0, len(profilesRef[0]))]
            Nidx = len(idx)
        self.b=1
        if self.bootstrapBy.get()==ProtPKPDDissolutionF2.BYVECTOR:
            for profileRef in profilesRef:
                for profileTest in profilesTest:
                    for n in range(self.Nbootstrap.get()):
                        idxB = sorted(np.random.choice(idx,Nidx))
                        while len(set(idxB))<3:
                            idxB = sorted(np.random.choice(idx, Nidx))
                        profileRefB = np.asarray([profileRef[i] for i in idxB])
                        profileTestB = np.asarray([profileTest[i] for i in idxB])
                        print("Bootstrap sample %d" % self.b)
                        print("calculateAllF pRef", profileRefB)
                        print("calculateAllF pTest", profileTestB)
                        f1, f2 = self.calculateF(profileRefB, profileTestB, True)
                        allF1b.append(f1)
                        allF2b.append(f2)
                        self.b = self.b + 1
        elif self.bootstrapBy.get() == ProtPKPDDissolutionF2.BYVECTOR:
            for n in range(self.Nbootstrap.get())*len(profilesRef)*len(profilesTest):
                profileRefB = self.bootstrapByTimePointAvg(profilesRef)
                profileTestB = self.bootstrapByTimePointAvg(profilesTest)
                f1, f2 = self.calculateF(profileRefB, profileTestB, True)
                allF1b.append(f1)
                allF2b.append(f2)
                self.b = self.b + 1
        else:
            for n in range(self.Nbootstrap.get())*len(profilesRef)*len(profilesTest):
                profileRefB = self.bootstrapByTimePoint(profilesRef)
                profileTestB = self.bootstrapByTimePoint(profilesTest)
                f1, f2 = self.calculateF(profileRefB, profileTestB, True)
                allF1b.append(f1)
                allF2b.append(f2)
                self.b = self.b + 1
        allF1b=[f1 for f1 in allF1b if not np.isnan(f1)]
        allF2b=[f2 for f2 in allF2b if not np.isnan(f2)]
        np.savetxt(self._getExtraPath("f1.txt"),allF1b)
        np.savetxt(self._getExtraPath("f2.txt"),allF2b)

        # Bias corrected and accelerated --------------------
        z0f1=np.clip(norm.ppf(float((np.asarray(allF1b)<f10).sum())/self.b),-10,10)
        z0f2=np.clip(norm.ppf(float((np.asarray(allF2b)<f20).sum())/self.b),-10,10)
        # print("f10",f10)
        # print("f20",f20)
        # print("self.b",self.b)
        # print("float(np.sum(allF1b<f10))",np.sum(allF1b<f10),(np.asarray(allF1b)<f10).sum(),float(np.sum(allF1b<f10)))
        # print("float(np.sum(allF2b<f20))",np.sum(allF2b<f20),(np.asarray(allF2b)<f20).sum(),float(np.sum(allF2b<f20)))
        # print("z0f1",z0f1)
        # print("z0f2",z0f2)
        af1=np.sum(np.power(allF1J-f1J,3.0))/(6*np.power(np.sum(np.power(allF1J-f1J,2.0)),1.5))
        af2=np.sum(np.power(allF2J-f2J,3.0))/(6*np.power(np.sum(np.power(allF2J-f2J,2.0)),1.5))
        # print("af1",af1)
        # print("af2",af2)
        alpha = 1-self.confidence.get()/100
        alphaLf1 = norm.cdf(z0f1+(z0f1+norm.ppf(alpha/2))/(1-af1*(z0f1+norm.ppf(alpha/2))))
        alphaLf2 = norm.cdf(z0f2+(z0f2+norm.ppf(alpha/2))/(1-af2*(z0f2+norm.ppf(alpha/2))))
        alphaUf1 = norm.cdf(z0f1+(z0f1+norm.ppf(1-alpha/2))/(1-af1*(z0f1+norm.ppf(1-alpha/2))))
        alphaUf2 = norm.cdf(z0f2+(z0f2+norm.ppf(1-alpha/2))/(1-af2*(z0f2+norm.ppf(1-alpha/2))))

        strF1Basic=self.printStats(allF1b,"Basic F1","sum(|pRef-pTest|)/sum(pRef)*100",alpha/2,1-alpha/2)
        strF2Basic=self.printStats(allF2b,"Basic F2","50*log10(100/sqrt(1+mean(|pRef-pTest|^2)))",alpha/2,1-alpha/2)

        strF1BCa=self.printStats(allF1b,"BCa F1","sum(|pRef-pTest|)/sum(pRef)*100",alphaLf1,alphaUf1)
        strF2BCa=self.printStats(allF2b,"Bca F2","50*log10(100/sqrt(1+mean(|pRef-pTest|^2)))",alphaLf2,alphaUf2)

        self.printSection("Results")
        fhSummary = open(self._getPath("summary.txt"),"w")

        self.doublePrint(fhSummary,"Sample F1 mean=%f"%f10)
        self.doublePrint(fhSummary,"Sample F2 mean=%f"%f20)
        self.doublePrint(fhSummary," ")

        self.doublePrint(fhSummary,"Jackknife F1 mean=%f"%f1J)
        self.doublePrint(fhSummary,"Jackknife F2 mean=%f"%f2J)
        self.doublePrint(fhSummary," ")

        print("Bias correction z0 F1=%f (bias)"%z0f1)
        print("Bias correction z0 F2=%f (bias)"%z0f2)
        print("Bias correction a F1=%f (acceleration)"%af1)
        print("Bias correction a F2=%f (acceleration)"%af2)
        print("Bias correction alphaL F1=%f"%alphaLf1)
        print("Bias correction alphaU F1=%f"%alphaUf1)
        print("Bias correction alphaL F2=%f"%alphaLf2)
        print("Bias correction alphaU F2=%f"%alphaUf2)
        print(" ")

        self.doublePrint(fhSummary,strF2Basic)
        self.doublePrint(fhSummary,strF2BCa)
        self.doublePrint(fhSummary,"---------------------------")
        self.doublePrint(fhSummary,strF1Basic)
        self.doublePrint(fhSummary,strF1BCa)
        fhSummary.close()

    def createOutputStep(self):
        if self.resampleT.get()>0 and self.keepResample.get():
            self._defineOutputs(outputExperiment=self.experiment)
            self._defineSourceRelation(self.inputRef, self.experiment)
            self._defineSourceRelation(self.inputTest, self.experiment)

    def _validate(self):
        return []

    def _summary(self):
        retval = []
        self.addFileContentToMessage(retval,self._getPath("summary.txt"))
        return retval

    def _citations(self):
        retval = ['Islam2018']
        return retval
