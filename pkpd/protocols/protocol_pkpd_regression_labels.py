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

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDVariable


class ProtPKPDRegressionLabel(ProtPKPD):
    """ Perform a regression between two labels\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'regression labels'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('labelY', params.StringParam, label="Y Label", default="", help='Name of a numerical label')
        form.addParam('labelX', params.StringParam, label="X Label", default="", help='Name of a numerical label')
        form.addParam('regressionType', params.EnumParam, choices=["Y=f(X)","Y=f(log(X))","log(Y)=f(X)","log(Y)=f(log(X))"],
                      label="Regression type", default=0,
                      help='f is always a polynomial of a given degree. log is the natural logarithm in base e\n')
        form.addParam('degree', params.IntParam, label="Degree", default=1, help='Degree of the polynomial f')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runRegression',self.inputExperiment.get().getObjId(), self.labelY.get(), self.labelX.get(),
                                 self.regressionType.get(), self.degree.get())

    #--------------------------- STEPS functions --------------------------------------------

    def runRegression(self, objId, labelY, labelX, regressionType, degree):
        X, Y, XtoUse, YtoUse = self.getXYValues()

        self.printSection("Performing regressions")

        p,V = np.polyfit(XtoUse,YtoUse,self.degree.get(),cov=True)
        mu = p
        Ypredicted = np.polyval(p, XtoUse)

        fhOut = open(self._getPath("results.txt"),'w')

        degree = self.degree.get()
        if self.regressionType.get()==0:
            modelToPrint = "%s="%labelY
            for i in range(0,degree+1):
                modelToPrint+="+(%f)*%s^(%d)"%(p[i],labelX,degree-i)
            self.doublePrint(fhOut,"Model: %s"%modelToPrint)
            self.doublePrint(fhOut," ")
            self.doublePrint(fhOut,"==========================================")
            self.doublePrint(fhOut,"%s    %s    %spredicted   Error=%s-%spredicted "%(labelX,labelY,labelY,labelY,labelY))
            self.doublePrint(fhOut,"==========================================")
            for i in range(0,len(X)):
                self.doublePrint(fhOut,"%f %f %f"%(X[i],Y[i],Ypredicted[i]))
        elif self.regressionType.get()==1:
            modelToPrint = "%s="%labelY
            for i in range(0,degree+1):
                modelToPrint+="+(%f)*(log(%s))^(%d)"%(p[i],labelX,degree-i)
            self.doublePrint(fhOut,"Model: %s"%modelToPrint)
            self.doublePrint(fhOut," ")
            self.doublePrint(fhOut,"=============================================")
            self.doublePrint(fhOut,"%s  log(%s)  %s    %spredicted Error=%s-%spredicted"%(labelX,labelX,labelY,labelY,labelY,labelY))
            self.doublePrint(fhOut,"=============================================")
            for i in range(0,len(X)):
                self.doublePrint(fhOut,"%f %f %f %f"%(X[i],math.log(X[i]),Y[i],Ypredicted[i]))
        elif self.regressionType.get()==2:
            modelToPrint = "log(%s)="%labelY
            for i in range(0,degree+1):
                modelToPrint+="+(%f)*%s^(%d)"%(p[i],labelX,degree-i)
            self.doublePrint(fhOut,"Model: %s"%modelToPrint)
            self.doublePrint(fhOut," ")
            self.doublePrint(fhOut,"============================================================================")
            self.doublePrint(fhOut,"%s     %s    log(%s)  log(%spredicted)  %spredicted  Error=log(%s)-log(%spredicted)"%\
                             (labelX,labelY,labelY,labelY,labelY,labelY,labelY,labelY))
            self.doublePrint(fhOut,"============================================================================")
            for i in range(0,len(X)):
                self.doublePrint(fhOut,"%f %f %f %f %f"%(X[i],Y[i],math.log(Y[i]), Ypredicted[i], math.exp(Ypredicted[i])))
        elif self.regressionType.get()==3:
            modelToPrint = "log(%s)="%labelY
            for i in range(0,degree+1):
                modelToPrint+="+(%f)*(log(%s))^(%d)"%(p[i],labelX,degree-i)
            self.doublePrint(fhOut,"Model: %s"%modelToPrint)
            self.doublePrint(fhOut," ")
            self.doublePrint(fhOut,"======================================================================================")
            self.doublePrint(fhOut,"%s     log(%s)     %s    log(%s)  log(%spredicted)  %spredicted Error=log(%s)-log(%spredicted)"%\
                             (labelX,labelX,labelY,labelY,labelY,labelY,labelY,labelY))
            self.doublePrint(fhOut,"======================================================================================")
            for i in range(0,len(X)):
                self.doublePrint(fhOut,"%f %f %f %f %f %f"%(X[i],math.log(X[i]),Y[i],math.log(Y[i]), Ypredicted[i], math.exp(Ypredicted[i])))

        e = YtoUse-Ypredicted
        R2 = (1-np.var(e)/np.var(YtoUse))
        n=XtoUse.shape[0] # Number of samples
        p=degree
        if n-p>0:
            R2adj = 1-R2*(n-1)/(n-p)*(1-R2)
        else:
            R2adj = -1
        logL = 0.5*(-n*(math.log(2*math.pi)+1-math.log(n)+math.log(np.sum(np.multiply(e,e)))))
        AIC = 2*p-2*logL
        if n-p-1>0:
            AICc = AIC+2*p*(p+1)/(n-p-1)
        else:
            AICc = -1
        BIC = p*math.log(n)-2*logL

        self.doublePrint(fhOut,"------------------------")
        self.doublePrint(fhOut,"Mean error = %f"%np.mean(e))
        self.doublePrint(fhOut,"Std error = %f"%np.std(e))
        self.doublePrint(fhOut,"R2 = %f"%R2)
        self.doublePrint(fhOut,"R2adj = %f"%R2adj)
        self.doublePrint(fhOut,"AIC = %f"%AIC)
        self.doublePrint(fhOut,"AICc(Recommended) = %f"%AICc)
        self.doublePrint(fhOut,"BIC = %f"%BIC)
        self.doublePrint(fhOut,"------------------------")

        from scipy.stats import norm
        nstd = norm.ppf(0.975)
        perr = np.sqrt(np.diag(V))
        lowerBound = mu-nstd*perr
        upperBound = mu+nstd*perr

        self.doublePrint(fhOut,"------------------------")
        for i in range(0,degree+1):
            self.doublePrint(fhOut,"Term of degree %d = %f Confidence interval 95%% [%f,%f]"%(degree-i,mu[i],lowerBound[i],upperBound[i]))
            if lowerBound[i]<0 and upperBound[i]>0:
                self.doublePrint(fhOut,"   This term may not be significant")
        self.doublePrint(fhOut,"------------------------")
        fhOut.close()

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=["%s=f(%s) where f is a polynomial of degree %d"%(self.labelY.get(),self.labelX.get(),self.degree.get())]
        self.addFileContentToMessage(msg,self._getPath("results.txt"))
        return msg

    def _validate(self):
        msg = []
        experiment = self.loadInputExperiment()

        label = self.labelY.get()
        if not label in experiment.variables:
            msg.append("%s is not a label of the experiment"%label)
        else:
            variable = experiment.variables[label]
            if variable.role != PKPDVariable.ROLE_LABEL:
                msg.append('%s is not a label'%label)
            else:
                if variable.varType!=PKPDVariable.TYPE_NUMERIC:
                    msg.append("%s is not numeric"%label)

        label = self.labelX.get()
        if not label in experiment.variables:
            msg.append("%s is not a label of the experiment"%label)
        else:
            variable = experiment.variables[label]
            if variable.role != PKPDVariable.ROLE_LABEL:
                msg.append('%s is not a label'%label)
            else:
                if variable.varType!=PKPDVariable.TYPE_NUMERIC:
                    msg.append("%s is not numeric"%label)
        return msg

    #--------------------------- UTILS functions --------------------------------------------
    def getXYValues(self, printExperiment=True):
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD,
                                         printExperiment)
        X = []
        Y = []
        labelX = self.labelX.get()
        labelY = self.labelY.get()
        for sampleName, sample in self.experiment.samples.iteritems():
            X.append(float(sample.descriptors[labelX]))
            Y.append(float(sample.descriptors[labelY]))
        X = np.asarray(X, np.double)
        Y = np.asarray(Y, np.double)
        logX = np.log(X)
        logY = np.log(Y)

        if self.regressionType.get() == 0 or self.regressionType.get() == 2:
            XtoUse = X
        else:
            XtoUse = logX

        if self.regressionType.get() == 0 or self.regressionType.get() == 1:
            YtoUse = Y
        else:
            YtoUse = logY

        return X, Y, XtoUse, YtoUse

    def getFunction(self):
        """ Parse the function string from the given result file.
        Replace the '^' char by the Python '**' operator. """
        result = None
        f = open(self._getPath('results.txt'))
        for line in f:
            if line.strip().startswith('Model:'):
                result = line.split('=')[1].replace('^', '**')
                break
        f.close()
        return result

    def evalFunction(self, xValues, func=None):
        """ Evaluate the function with the given values
        and return the list of resulting values. """
        func = func or self.getFunction()
        xLabel = self.labelX.get()
        evalDict = {'log': np.log}

        yValues = []
        for x in xValues:
            evalDict[xLabel] = x
            yValues.append(eval(func, globals(), evalDict))

        return yValues

    def filterVarForWizard(self, v):
        """ Define the type of variables required (used in wizard). """
        return v.isNumeric() and v.isLabel()
