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
from .protocol_pkpd_fit_base import ProtPKPDFitBase
from pkpd.models.dissolution_models import *


class ProtPKPDDissolutionFit(ProtPKPDFitBase):
    """ Fit a dissolution model. The observed measurement is modelled as Y=f(t).\n
Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'fit dissolution'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineParams1(form,"t","C")
        form.addParam('allowTlag', params.BooleanParam,
                      label="Allow lag", default=False,
                      help='Allow lag time before starting dissolution (t-tlag)')
        form.addParam('modelType', params.EnumParam, choices=["Zero order","First order","Fractional","Weibull","Higuchi",
                                                              "Korsmeyer-Peppas","Hixson-Crowell","Hopfenberg"],
                      label="Dissolution model", default=3,
                      help='Zero order: Y=K*(t-[tlag])\n'\
                           'First order: Y=Ymax*(1-exp(-beta*(t-[tlag])))\n'\
                           'Fractional order: Y=Ymax-pow(Amax^alpha-alpha*beta*t,1/alpha))\n'\
                           'Weibull: Y=Ymax*(1-exp(-lambda*t^b))\n'\
                           'Higuchi: Y=Ymax*t^0.5\n'\
                           'Korsmeyer-Peppas: Y=Ymax*t^m\n'\
                           'Hixson-Crowell: Y=Ymax*(1-(1-K*t)^3)\n'
                           'Hopfenberg: Y=Ymax*(1-(1-K*t)^m)')
        form.addParam('fitType', params.EnumParam, choices=["Linear","Logarithmic","Relative"], label="Fit mode", default=0,
                      expertLevel=LEVEL_ADVANCED,
                      help='Linear: sum (Cobserved-Cpredicted)^2\nLogarithmic: sum(log10(Cobserved)-log10(Cpredicted))^2\n'\
                           "Relative: sum ((Cobserved-Cpredicted)/Cobserved)^2")
        form.addParam('bounds', params.StringParam, label="Bounds (optional)", default="",
                      help='Parameter values for the simulation.\nExample: (1,10);(0,0.05) is (1,10) for the first parameter, (0,0.05) for the second parameter\n'
                           'Zero order: [tlag];K\n'
                           'First order: [tlag];Ymax;beta\n'
                           'Fractional order: [tlag]; Ymax;beta;alpha\n'
                           'Weibull: [tlag]; Ymax;lambda;b\n'
                           'Higuchi: [tlag]; Ymax\n'
                           'Korsmeyer-Peppas: [tlag]; Ymax; m\n'
                           'Hixson-Crowell: [tlag]; Ymax; K\n'
                           'Hopfenberg: [tlag]; Ymax; K; m')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval=", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')

    def getListOfFormDependencies(self):
        return [self.allowTlag.get(), self.modelType.get(), self.bounds.get(), self.confidenceInterval.get()]

    #--------------------------- STEPS functions --------------------------------------------
    def createModel(self):
        if self.modelType.get() == 0:
            return Dissolution0()
        elif self.modelType.get() == 1:
            return Dissolution1()
        elif self.modelType.get() == 2:
            return DissolutionAlpha()
        elif self.modelType.get() == 3:
            return DissolutionWeibull()
        elif self.modelType.get() == 4:
            return DissolutionHiguchi()
        elif self.modelType.get() == 5:
            return DissolutionKorsmeyer()
        elif self.modelType.get() == 6:
            return DissolutionHixson()
        elif self.modelType.get() == 7:
            return DissolutionHopfenberg()

    def setupFromFormParameters(self):
        self.model.allowTlag = self.allowTlag.get()
