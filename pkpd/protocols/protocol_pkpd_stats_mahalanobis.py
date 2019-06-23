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
from scipy import stats
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.spatial import distance


import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.utils import upper_tri_masking, uniqueFloatValues

# Tested in test_workflow_levyplot.py

class ProtPKPDStatsMahalanobis(ProtPKPD):
    """ Experiment 1 defines the mean and covariance for the Mahalanobis distance.
        Then, the Mahalanobis distance of all elements in Experiment 1 with respect to the mean is calculated\n
        If a second experiment is given, then all distances from the second to the mean of the first experiment\n
        are also calculated.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'Mahalanobis'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment1', params.PointerParam, label="Experiment 1",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('labels', params.StringParam, label="Labels", default="",
                      help='Name of the labels to construct the multivariate vectors, e.g. Cl V\n'
                           'The set of labels is either a set of labels or a single measurement.\n'
                           'If you use a measurement, make sure that it does not have a measurement\n'
                           'that is always the same, typically (0,0)')
        form.addParam('expression1', params.StringParam, label="Subgroup 1 (optional)", default="",
                      help='For example, $(weight)<100 and $(sex)=="male"')
        form.addParam('inputExperiment2', params.PointerParam, label="Experiment 2 (optional)",
                      pointerClass='PKPDExperiment', allowsNull=True,
                      help='Select an experiment with samples')
        form.addParam('expression2', params.StringParam, label="Subgroup 2 (optional)", default="",
                      help='For example, $(weight)>=100 and $(sex)=="male". If it is empty, the same Expression 1 will be used for grouping in Experiment 2')
        form.addParam('resampleT', params.FloatParam, label="Resample profiles (time step)", default=-1,
                      help='Resample the input profiles at this time step (make sure it is in the same units as the input). '
                           'Leave it to -1 for no resampling. This is only valid when the label to compare is a measurement.')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runCompare',self.inputExperiment1.get().getObjId(),
                                 self.labels.get(), self.expression1.get(), self.expression2.get())

    #--------------------------- STEPS functions --------------------------------------------
    def getLabels(self):
        labels=[]
        for token in self.labels.get().split():
            labels.append(token.strip())
        return labels

    def analysisType(self, experiment):
        Nlabels=0
        Nmeasurements=0
        N=len(self.getLabels())
        for label in self.getLabels():
            if label in experiment.variables:
                if experiment.variables[label].isLabel():
                    Nlabels+=1
                elif experiment.variables[label].isMeasurement():
                    Nmeasurements+=1
            else:
                raise Exception("Cannot find the label %s in one of the experiments"%label)
        if Nlabels==N:
            return 0
        if Nmeasurements==1 and N==1:
            return 1
        raise Exception("The set of labels is not correctly formed. It should be all labels or one measurement")

    def getValues(self,experiment,expression):
        allX=[]
        if self.analysisType(experiment)==0: # All labels
            for label in self.getLabels():
                x1 = [float(x) for x in experiment.getSubGroupLabels(expression, label)]
                allX.append(x1)
            X = np.asarray(allX, dtype=np.double)
            return np.transpose(X)
        else: # A measurement
            Ylabel=self.getLabels()[0]
            Xlabel=experiment.getTimeVariable()
            temp=[]
            minX=1e38
            maxX=-1e38
            for sampleName, sample in experiment.samples.iteritems():
                x,y=sample.getXYValues(Xlabel,Ylabel)
                minX=min(minX,np.min(x))
                maxX=max(maxX,np.max(x))
                temp.append((x[0],y[0]))

            for x,y in temp:
                if self.resampleT.get()>0:
                    x,y = uniqueFloatValues(x, y)
                    B = InterpolatedUnivariateSpline(x, y, k=1)
                    xp = np.arange(minX, maxX + self.resampleT.get(), self.resampleT.get())
                    y = B(xp)
                allX.append(y)
            X = np.asarray(allX, dtype=np.double)

            v = np.var(X,0)
            X=X[:,v>0] # Remove columns with no variance
            return X

    def printStats(self,allF,Fstr,explanation):
        allF=[f for f in allF if not np.isnan(f)]
        mu=np.mean(allF)
        sigma = np.std(allF)
        alpha=1-95/100.0
        percentiles = np.percentile(allF,[0, alpha/2*100, 25, 50, 75, (1-alpha/2)*100, 100])
        retval=""
        retval +="%s (%s)\n"%(Fstr,explanation)
        retval +="%s mean+-std: %f+-%f\n"%(Fstr,mu,sigma)
        retval +="%s minimum,maximum: [%f,%f]\n"%(Fstr,percentiles[0],percentiles[6])
        retval +="%s percentile [%f,%f]%%: [%f,%f]\n"%(Fstr,alpha/2*100,(1-alpha/2)*100,percentiles[1],percentiles[5])
        retval +="%s percentile [25,75]%%: [%f,%f]\n"%(Fstr,percentiles[2],percentiles[4])
        retval +="%s percentile 50%%: %f\n"%(Fstr,percentiles[3])
        return retval

    def runCompare(self, objId1, labels, expression1, expression2):
        fh = open(self._getPath("report.txt"),'w')

        self.experiment1 = self.readExperiment(self.inputExperiment1.get().fnPKPD, False)
        X1 = self.getValues(self.experiment1,self.expression1.get())
        X2 = None
        if self.inputExperiment2.get() is not None:
            self.experiment2 = self.readExperiment(self.inputExperiment2.get().fnPKPD, False)
            expression2ToUse = self.expression1.get() if self.expression2.get()=="" else self.expression2.get()
            X2 = self.getValues(self.experiment2, expression2ToUse)

        print("Values in SubGroup 1:\n%s"%np.array2string(X1))
        if X2 is not None:
            print("Values in SubGroup 2:\n%s" % np.array2string(X2))

        try:
            self.printSection("Results")
            C1 = np.cov(np.transpose(X1))
            if np.abs(np.linalg.det(C1))<1e-10:
                distanceStr='Euclidean'
                self.doublePrint(fh,"The covariance matrix is singular (either there is a column of the data that is always the same or not enough data to define the covariance)")
                self.doublePrint(fh,"Using Euclidean distance instead")
                D11 = upper_tri_masking(distance.cdist(X1, X1, 'euclidean'))
                D11/= X1.shape[1]
            else:
                distanceStr='Mahalanobis'
                C1inv = np.linalg.inv(C1)
                D11 = upper_tri_masking(distance.cdist(X1,X1,'mahalanobis',VI=C1inv))

            np.savetxt(self._getExtraPath("D11.txt"),D11)

            str11 = self.printStats(D11, "D11", "%s distance Set 1 vs Set1"%distanceStr)
            self.doublePrint(fh, str11)
            if X2 is not None:
                self.doublePrint(fh, "---------------------------")
                if distanceStr=='Mahalanobis':
                    D12 = distance.cdist(X1,X2,'mahalanobis',VI=C1inv).flatten()
                    D12 /= X1.shape[1]
                else:
                    D12 = distance.cdist(X1, X2, 'euclidean').flatten()
                np.savetxt(self._getExtraPath("D12.txt"),D12)

                str12 = self.printStats(D12, "D12", "%s distance Set 1 vs Set2"%distanceStr)
                self.doublePrint(fh, str12)

                [D,pval] = stats.ks_2samp(D11, D12)
                self.doublePrint(fh, "---------------------------")
                self.doublePrint(fh,"Kolmogorov-Smirnov test for the compatibility of D11 and D12: D-statistic=%f p-value=%f"%(D,pval))
            fh.close()
        except Exception as e:
            print(e)

        fh.close()

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        expression2ToUse = self.expression1.get() if self.expression2.get()=="" else self.expression2.get()
        msg=["Comparison between %s in Subgroup1 (%s) and %s in Subgroup2 (%s), independent samples"%(self.labels.get(),self.expression1.get(),self.labels.get(),expression2ToUse)]
        msg.append(' ')
        self.addFileContentToMessage(msg,self._getPath("report.txt"))
        return msg

    def _validate(self):
        msg=[]
        # self.experiment1 = self.readExperiment(self.inputExperiment1.get().fnPKPD,False)
        # if self.inputExperiment2.get() is not None:
        #     self.experiment2 = self.readExperiment(self.inputExperiment2.get().fnPKPD)
        # for label in self.getLabels():
        #     if not label in self.experiment1.variables:
        #         msg.append("Cannot find %s amongst the Experiment 1 variables"%label)
        #     if self.inputExperiment2.get() is not None:
        #         if not label in self.experiment2.variables:
        #             msg.append("Cannot find %s amongst the Experiment 2 variables" % label)
        return msg
