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
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.object import String, Integer
from pkpd.models.pk_models import PKPDExponentialModel
from .protocol_pkpd_fit_base import ProtPKPDFitBase


# TESTED in test_workflow_gabrielsson_pk01.py
# TESTED in test_workflow_gabrielsson_pk02.py
# TESTED in test_workflow_gabrielsson_pk07.py

class ProtPKPDExponentialFit(ProtPKPDFitBase):
    """ Fit a set of exponentials. The observed measurement is modelled as Y=sum_{i=1}^N c_i exp(-lambda_i * X).\n
Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'fit exponentials'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form, fullForm=True):
        self._defineParams1(form,"t","Cp")
        if fullForm:
            form.addParam('fitType', params.EnumParam, choices=["Linear","Logarithmic","Relative"], label="Fit mode", default=1,
                          help='Linear: sum (Cobserved-Cpredicted)^2\nLogarithmic: sum(log10(Cobserved)-log10(Cpredicted))^2\n'\
                               "Relative: sum ((Cobserved-Cpredicted)/Cobserved)^2")
            form.addParam('Nexp', params.IntParam, label="Number of exponentials", default=1,
                          help='Number of exponentials to fit')
        else:
            self.fitType=Integer()
            self.fitType.set(1)
            self.Nexp=Integer()
            self.Nexp.set(1)
        form.addParam('bounds', params.StringParam, label="Amplitude and time constant bounds", default="", expertLevel=LEVEL_ADVANCED,
                      help='Bounds for the c_i amplitudes and lambdas.\nExample 1: (0,10);(0,1e-2) -> c1 in (0,10), lambda1 in (0,1e-2)\n'\
                           'Example 2: (0,10);(0,1e-2);(0,1);(0,1e-1) -> c1 in (0,10), lambda1 in (0,1e-2), c2 in (0,1), lambda2 in (0,1e-1)')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval=", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')
        if fullForm:
            form.addParam('reportX', params.StringParam, label="Evaluate at X=", default="", expertLevel=LEVEL_ADVANCED,
                          help='Evaluate the model at these X values\nExample 1: [0,5,10,20,40,100]\nExample 2: 0:0.55:10, from 0 to 10 in steps of 0.5')
        else:
            self.reportX=String()
            self.reportX.set("")

    def getListOfFormDependencies(self):
        return [self.predictor.get(), self.predicted.get(), self.fitType.get(), self.bounds.get()]

    def createModel(self):
        return PKPDExponentialModel()

    def setupFromFormParameters(self):
        self.model.Nexp=self.Nexp.get()

    #--------------------------- INFO functions --------------------------------------------
    def _warnings(self):
        warnings = []
        experiment = self.readExperiment(self.getInputExperiment().fnPKPD,show=False)
        incorrectList = experiment.getNonBolusDoses()
        if len(incorrectList)>0:
            warnings.append("This protocol is meant only for intravenous bolus regimens. Check the doses for %s"%(','.join(incorrectList)))
        return warnings

    def _validate(self):
        errors=ProtPKPDFitBase._validate(self)
        if self.Nexp.get()<1:
            errors.append("The number of exponentials has to be larger than 0")
        return errors
