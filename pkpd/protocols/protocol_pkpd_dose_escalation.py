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
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDDoseResponse
import scipy.stats
from pyworkflow.protocol.constants import LEVEL_ADVANCED

def ClopperPearson(k,n,alpha=0.05):
    # https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    lo = scipy.stats.beta.ppf(alpha / 2, k, n - k + 1)
    hi = scipy.stats.beta.ppf(1 - alpha / 2, k + 1, n - k)
    if np.isnan(lo):
        lo = 0.0
    if np.isnan(hi):
        hi = 1.0
    return [lo,hi]

class ProtPKPDDoseEscalation(ProtPKPD):
    """ Given a set of binary responses (toxicity, response/not response, ...), estimate the next dose for a target response\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'dose escalation'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form, fullForm=True):
        form.addSection('Input')
        form.addParam('prot1ptr', params.PointerParam, label="Previous dose response",
                      pointerClass='ProtPKPDSimulateDoseEscalation,ProtPKPDDoseEscalation',allowsNull=True,
                      help='Optional. You may link this simulation to a previous dose simulation or analysis, so that you may know the sequence')
        form.addParam('measurements', params.TextParam, height=12, width=80, label="Measurements", default="",
                      help="If you give a previous dose response, these doses will be added to the previous ones. If you are coming from a escalation simulation, it is not necessary to fill this field.\n"\
                           "Each dose should occupy a line. The dose amount goes first, then semicolon and 0 or 1 depending on the response. Example\n" \
                           "0.05:  0 0 0\n" \
                           "0.10:  0 0 0\n" \
                           "0.167: 0 1 0\n" \
                      )
        form.addParam('scaleFactorsForm', params.StringParam, label="Dose scale factors", default="1, 2, 1.67, 1.4, 1.33", expertLevel=LEVEL_ADVANCED,
                      help="Dose escalation sequence. The default value is d1=1*d1, d2=2*d1, d3=1.67*d2, d4=1.4*d3, d5=1.33*d4, d6=1.33*d5, ..., dn=1.33*dn-" \
                      )

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runEstimate')

    #--------------------------- STEPS functions --------------------------------------------
    def runEstimate(self):
        self.parseInput()
        self.threePlusThree()
        self.bestOfFive()
        self.upAndDown()
        self.storerC()
        self.storerBC()
        self.dr.write(self._getPath("dose_response.txt"))

    def getPreviousDose(self,n0=1):
        n = self.findDose(self.dr.doses[-1])
        if n>=n0:
            return self.doseList[n-n0]
        else:
            return self.doseList[0]

    def getCurrentDose(self):
        return self.dr.doses[-1]

    def getNextDose(self):
        return self.doseList[self.findDose(self.dr.doses[-1])+1]

    def findDose(self,dose):
        for n in range(len(self.doseList)):
            if abs(self.doseList[n]-dose)<1e-6:
                return n
        return None

    def getLastResponse(self):
        return self.dr.responses[-1].split(" ")[-1]

    def getBeforeLastResponse(self):
        allResponses = " ".join(self.dr.responses)
        tokens = allResponses.split(" ")
        if len(tokens)>=2:
            return tokens[-2]
        else:
            return None

    def parseInput(self):
        self.dr=PKPDDoseResponse()
        if self.prot1ptr.get() is not None:
            self.dr.read(self.prot1ptr.get()._getPath("dose_response.txt"),True)
        self.dr.read(self.measurements.get())
        print("Current responses ===============================================================")
        for i in range(len(self.dr.doses)):
            print("Dose %f DLT Responses: %s"%(self.dr.doses[i],self.dr.responses[i]))

        print("\n\nAnalysis of confidence intervals ================================================")
        for i in range(len(self.dr.doses)):
            k=self.dr.Nresponses[i]
            n=self.dr.Npatients[i]
            lo,hi=ClopperPearson(k,n)
            print("Dose=%f Prob=%d/%d=%f. 95%% Confidence interval: [%f,%f]" % (self.dr.doses[i],k,n,float(k)/float(n),lo,hi))

        self.scaleFactors=[]
        for token in self.scaleFactorsForm.get().split(','):
            self.scaleFactors.append(float(token.strip()))

        self.doseList=[]
        # Append all doses that have been administered so far
        for dose in self.dr.doses:
            if self.findDose(dose) is None:
                self.doseList.append(dose)
        # Append one more dose
        self.doseList = sorted(self.doseList)
        scaleFactor = self.scaleFactors[min(len(self.scaleFactors)-1,len(self.doseList))]
        self.doseList.append(self.doseList[-1]*scaleFactor)

    def threePlusThree(self):
        print("\n\n3+3 Strategy ====================================================================")
        Nresponse=self.dr.Nresponses[-1]
        Npatients=self.dr.Npatients[-1]
        if Npatients==3:
            if Nresponse==0:
                print("0 DLTs in 3. Go to next Dose=%f"%self.getNextDose())
            elif Nresponse==1:
                print("1 DLT in 3. Try 3 more patients at this Dose=%f"%self.getCurrentDose())
            else:
                print("Discontinue escalation. Maximum Tolerable Dose=%f"%self.getPreviousDose())
        elif Npatients==6:
            if Nresponse==1:
                print("1 DLT in 6. Go to next Dose=%f"%self.getNextDose())
            else:
                print("Discontinue escalation. Maximum Tolerable Dose=%f"%self.getPreviousDose())
        else:
            print("This strategy does not know how to handle this situation")

    def bestOfFive(self):
        print("\n\nBest of 5 Strategy ==============================================================")
        Nresponse=self.dr.Nresponses[-1]
        Npatients=self.dr.Npatients[-1]
        if Npatients==3:
            if Nresponse==0:
                print("0 DLTs in 3. Go to next Dose=%f"%self.getNextDose())
            elif Nresponse==1:
                print("1 DLT in 3. Try 1 more patient at this Dose=%f"%self.getCurrentDose())
            else:
                print("Discontinue escalation. Maximum Tolerable Dose=%f"%self.getPreviousDose())
        elif Npatients==4:
            if Nresponse==1:
                print("1 DLT in 4. Go to next dose d=%f"%self.getNextDose())
            elif Nresponse == 2:
                print("2 DLTs in 3. Try 1 more patient at this Dose=%f"%self.getCurrentDose())
            else:
                print("Discontinue escalation. Maximum Tolerable Dose=%f"%self.getPreviousDose())
        elif Npatients==5:
            if Nresponse==2:
                print("2 DLTs in 5. Go to next Dose=%f"%self.getNextDose())
            else:
                print("Discontinue escalation. Maximum Tolerable Dose=%f"%self.getPreviousDose())
        else:
            print("This strategy does not know how to handle this situation")

    def upAndDown(self):
        print("\n\nUp And Down Strategy ============================================================")
        lastResponse = self.getLastResponse()
        NTotalPatients = np.sum(self.dr.Npatients)
        print("Current number of patients: %d"%NTotalPatients)
        if lastResponse=="0":
            doseNextIndividual = self.getNextDose()
            print("Last individual did not have DLT. Go to next Dose=%f"%doseNextIndividual)
        else:
            doseNextIndividual = self.getPreviousDose()
            print("Last individual had DLT. Go to previous Dose=%f"%doseNextIndividual)
        print("Current Maximum Tolerable Dose=%f"%doseNextIndividual)
        print("This strategy should be run up to a prespecified total number of patients")

    def storerC(self):
        print("\n\nStorer's C Strategy ============================================================")
        lastResponse = self.getLastResponse()
        NTotalPatients = np.sum(self.dr.Npatients)
        print("Current number of patients: %d"%NTotalPatients)
        if lastResponse=="0":
            if self.dr.Npatients[-1]==1:
                doseNextIndividual = self.getCurrentDose()
                print("Last patient with no DLT. Try 1 more patient at this Dose=%f" % doseNextIndividual)
            elif self.dr.Npatients[-1]>=2:
                beforeLastResponse = self.getBeforeLastResponse()
                if beforeLastResponse=="0":
                    doseNextIndividual = self.getNextDose()
                    print("Last 2 patients did not have DLT. Go to next Dose d=%f" % doseNextIndividual)
                else:
                    doseNextIndividual = -1
                    print("This strategy cannot handle this situation")
        else:
            doseNextIndividual = self.getPreviousDose()
            print("Last individual had DLT. Go to previous Dose=%f"%doseNextIndividual)
        print("Current Maximum Tolerable Dose=%f"%doseNextIndividual)
        print("This strategy should be run up to a prespecified total number of patients")

    def storerBC(self):
        print("\n\nStorer's Two-stage (BC) Strategy ================================================")
        responsesToThisDose = self.dr.responses[-1].split(" ")
        lastResponse = self.getLastResponse()
        NTotalPatients = np.sum(self.dr.Npatients)
        NTotalResponses = np.sum(self.dr.Nresponses)
        print("Current number of patients: %d"%NTotalPatients)
        if NTotalResponses==0: # Stage 1
            doseNextIndividual = self.getNextDose()
            print("Last patient did not have DLT. Go to next Dose d=%f" % doseNextIndividual)
        elif NTotalResponses==1 and lastResponse=="1": # End of stage 1
            doseNextIndividual = self.getPreviousDose()
            print("Last patient had DLT. Go to next Dose d=%f" % doseNextIndividual)
        else: # Stage 2
            if lastResponse=="0":
                if self.dr.Npatients[-1]==1:
                    doseNextIndividual = self.getCurrentDose()
                    print("Last patient with no DLT. Try 1 more patient at this Dose=%f" % doseNextIndividual)
                elif self.dr.Npatients[-1]>=2:
                    beforeLastResponse = self.getBeforeLastResponse()
                    if beforeLastResponse=="0":
                        doseNextIndividual = self.getNextDose()
                        print("Last 2 patients did not have DLT. Go to next Dose d=%f" % doseNextIndividual)
                    else:
                        doseNextIndividual = -1
                        print("This strategy cannot handle this situation")
            else:
                beforeLastResponse = self.getBeforeLastResponse()
                if beforeLastResponse=="1":
                    doseNextIndividual = self.getPreviousDose(2)
                    print("Last two patients had DLT. Go to previous Dose=%f" % doseNextIndividual)
                else:
                    doseNextIndividual = self.getPreviousDose()
                    print("Last individual had DLT. Go to previous Dose=%f"%doseNextIndividual)
        print("Current Maximum Tolerable Dose=%f"%doseNextIndividual)
        print("This strategy should be run up to a prespecified total number of patients")

        # bcrm, CRM, dfcrm, crmpack

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        return [self.measurements.get()]
