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
from pkpd.objects import PKPDExperiment, PKPDVariable
from scipy import stats
import numpy as np

# Tested in test_workflow_levyplot.py

class ProtPKPDStatsExp2Subgroups2Mean(ProtPKPD):
    """ Compare two means from two subgroups from the same experiment .\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'Exp2 SubGr2 Mean'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment1', params.PointerParam, label="Experiment 1", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('label1', params.StringParam, label="Label 1", default="",
                      help='Name of the label in the first experiment to compare between the two subroups')
        form.addParam('expression1', params.StringParam, label="Subgroup 1 (optional)", default="",
                      help='For example, $(weight)<100 and $(sex)=="male"')
        form.addParam('inputExperiment2', params.PointerParam, label="Experiment 2", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('label2', params.StringParam, label="Label 2", default="",
                      help='Name of the label in the second experiment to compare between the two subroups. '
                      'If it is empty, Label 1 is also used as Label 2')
        form.addParam('paired', params.BooleanParam, label="Paired samples?", default=False,
                      help='If samples are paired, then the sample names found in Experiment 1 will be sought in Experiment 2')
        form.addParam('expression2', params.StringParam, label="Subgroup 2 (optional)", default="", condition="not paired",
                      help='For example, $(weight)>=100 and $(sex)=="male". If it is empty, the same Expression 1 will be used for grouping in Experiment 2')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runCompare',self.inputExperiment1.get().getObjId(), self.inputExperiment2.get().getObjId(),
                                 self.label1.get(), self.label2.get(), self.paired.get(),
                                 self.expression1.get(), self.expression2.get())

    #--------------------------- STEPS functions --------------------------------------------
    def runCompare(self, objId1, objId2, label1, label2, paired, expression1, expression2):
        fh = open(self._getPath("report.txt"),'w')

        self.experiment1 = self.readExperiment(self.inputExperiment1.get().fnPKPD)
        self.experiment2 = self.readExperiment(self.inputExperiment2.get().fnPKPD)
        label2ToUse = self.label1.get() if self.label2.get()=="" else self.label2.get()
        if self.paired:
            x1=[]
            x2=[]
            for sampleName, sample in self.experiment1.getSubGroup(self.expression1.get()).iteritems():
                x1.append(float(sample.descriptors[self.label1.get()]))
                if sampleName in self.experiment2.samples:
                    x2.append(float(self.experiment2.samples[sampleName].descriptors[label2ToUse]))
                else:
                    raise "Cannot find sample %s in Experiment 2"%sample.sampleName
        else:
            expression2ToUse = self.expression1.get() if self.expression2.get()=="" else self.expression2.get()
            x1 = [float(x) for x in self.experiment1.getSubGroupLabels(self.expression1.get(),self.label1.get())]
            x2 = [float(x) for x in self.experiment2.getSubGroupLabels(expression2ToUse,label2ToUse)]
        self.doublePrint(fh,"Values in SubGroup 1: %s"%str(x1))
        self.doublePrint(fh,"Values in SubGroup 2: %s"%str(x2))
        self.doublePrint(fh,"Testing H0: mu1=mu2")
        self.doublePrint(fh," ")

        try:
            if self.paired:
                [t,pval] = stats.ttest_rel(np.asarray(x1,np.double),np.asarray(x2,np.double),True)
                self.doublePrint(fh,"T-test two paired samples: t-statistic=%f p-value=%f"%(t,pval))
            else:
                [t,pval] = stats.ttest_ind_from_stats(np.mean(x1),np.std(x1),len(x1),
                                                      np.mean(x2),np.std(x2),len(x2),True)
                self.doublePrint(fh,"T-test two independent samples (same variance): t-statistic=%f p-value=%f"%(t,pval))
        except Exception as e:
            print(e)

        if not self.paired:
            try:
                [t,pval] = stats.ttest_ind_from_stats(np.mean(x1),np.std(x1),len(x1),
                                                      np.mean(x2),np.std(x2),len(x2),False)
                self.doublePrint(fh,"T-test two independent samples (different variance, Welch's test): t-statistic=%f p-value=%f"%(t,pval))
            except Exception as e:
                print(e)

        try:
            if self.paired:
                [w,pval] = stats.wilcoxon(x1, x2, correction=True)
                self.doublePrint(fh,"Wilcoxon signed rank test for two paired samples: w-statistic=%f p-value=%f"%(w,pval))
            else:
                [u,pval] = stats.mannwhitneyu(x1, x2, True)
                self.doublePrint(fh,"Mann-Whitney U test for two independent samples: u-statistic=%f p-value=%f"%(u,pval))
        except:
            pass

        fh.close()

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        import os
        label2ToUse = self.label1.get() if self.label2.get()=="" else self.label2.get()
        if self.paired:
            msg=["Comparison between %s in Subgroup1 (%s) and %s in Subgroup2, paired samples"%(self.label1.get(),self.expression1.get(),label2ToUse)]
        else:
            expression2ToUse = self.expression1.get() if self.expression2.get()=="" else self.expression2.get()
            msg=["Comparison between %s in Subgroup1 (%s) and %s in Subgroup2 (%s), independent samples"%(self.label1.get(),self.expression1.get(),label2ToUse,expression2ToUse)]
        msg.append(' ')
        self.addFileContentToMessage(msg,self._getPath("report.txt"))
        return msg

    def _validate(self):
        msg=[]
        self.experiment1 = self.readExperiment(self.inputExperiment1.get().fnPKPD,False)
        if not self.label1.get() in self.experiment1.variables:
            msg.append("Cannot find %s amongst the Experiment 1 variables"%self.label1.get())
        else:
            variable = self.experiment1.variables[self.label1.get()]
            if variable.role != PKPDVariable.ROLE_LABEL:
                msg.append("Variable %s is not a label in Experiment 1"%self.label1.get())
            if variable.varType != PKPDVariable.TYPE_NUMERIC:
                msg.append("Variable %s is not a number in Experiment 1"%self.label1.get())

        label2ToUse = self.label1.get() if self.label2.get()=="" else self.label2.get()
        self.experiment2 = self.readExperiment(self.inputExperiment2.get().fnPKPD,False)
        if not label2ToUse in self.experiment2.variables:
            msg.append("Cannot find %s amongst the Experiment 2 variables"%label2ToUse)
        else:
            variable = self.experiment1.variables[label2ToUse]
            if variable.role != PKPDVariable.ROLE_LABEL:
                msg.append("Variable %s is not a label in Experiment 2"%label2ToUse)
            if variable.varType != PKPDVariable.TYPE_NUMERIC:
                msg.append("Variable %s is not a number in Experiment 2"%label2ToUse)
        return msg

    def filterVarForWizard(self, v):
        """ Define the type of variables required (used in wizard). """
        return v.isNumeric() and v.isLabel()