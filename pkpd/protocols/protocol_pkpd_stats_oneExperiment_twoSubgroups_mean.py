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

class ProtPKPDStatsExp1Subgroups2Mean(ProtPKPD):
    """ Compare two means from two subgroups from the same experiment .\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'Exp1 SubGr2 Mean'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment", important=True,
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('labelToCompare', params.StringParam, label="Label to compare", default="",
                      help='Name of the label to compare between the two subroups')
        form.addParam('expression1', params.StringParam, label="Subgroup 1", default="",
                      help='For example, $(weight)<100 and $(sex)=="male"')
        form.addParam('expression2', params.StringParam, label="Subgroup 2", default="",
                      help='For example, $(weight)>=100 and $(sex)=="male"')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runCompare',self.inputExperiment.get().getObjId(), self.labelToCompare.get(),
                                 self.expression1.get(), self.expression2.get())

    #--------------------------- STEPS functions --------------------------------------------
    def runCompare(self, objId, labelToAdd, expression1, expression2):
        fh = open(self._getPath("report.txt"),'w')

        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        x1 = [float(x) for x in self.experiment.getSubGroupLabels(self.expression1.get(),self.labelToCompare.get())]
        x2 = [float(x) for x in self.experiment.getSubGroupLabels(self.expression2.get(),self.labelToCompare.get())]
        self.doublePrint(fh,"Values in SubGroup 1: %s"%str(x1))
        self.doublePrint(fh,"Values in SubGroup 2: %s"%str(x2))
        self.doublePrint(fh,"Testing H0: mu1=mu2")
        self.doublePrint(fh," ")

        try:
            [t,pval] = stats.ttest_ind(np.asarray(x1,np.double),np.asarray(x2,np.double),True)
            self.doublePrint(fh,"T-test two independent samples (same variance): t-statistic=%f p-value=%f"%(t,pval))
        except:
            pass

        try:
            [t,pval] = stats.ttest_ind(x1,x2, False)
            self.doublePrint(fh,"T-test two independent samples (different variance, Welch's test): t-statistic=%f p-value=%f"%(t,pval))
        except:
            pass

        try:
            [u,pval] = stats.mannwhitneyu(x1, x2, True)
            self.doublePrint(fh,"Mann-Whitney U test for two independent samples: u-statistic=%f p-value=%f"%(u,pval))
        except:
            pass

        fh.close()

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        import os
        msg=["Comparison between Subgroup1 (%s) and Subgroup2 (%s)"%(self.expression1.get(),self.expression2.get())]
        msg.append(' ')
        self.addFileContentToMessage(msg,self._getPath("report.txt"))
        return msg

    def _validate(self):
        msg=[]
        self.experiment = self.readExperiment(self.inputExperiment.get().fnPKPD,False)
        if not self.labelToCompare.get() in self.experiment.variables:
            msg.append("Cannot find %s amongst the experiment variables"%self.labelToCompare.get())
        else:
            variable = self.experiment.variables[self.labelToCompare.get()]
            if variable.role != PKPDVariable.ROLE_LABEL:
                msg.append("Variable %s is not a label"%self.labelToCompare.get())
            if variable.varType != PKPDVariable.TYPE_NUMERIC:
                msg.append("Variable %s is not a number"%self.labelToCompare.get())
        return msg

    def filterVarForWizard(self, v):
        """ Define the type of variables required (used in wizard). """
        return v.isNumeric() and v.isLabel()