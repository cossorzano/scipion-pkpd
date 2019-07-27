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
from pkpd.objects import PKPDExperiment


class ProtPKPDCompareExperiments(ProtPKPD):
    """ This protocol compares two experiments by plotting a summary of both.
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'compare experiments'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment1', params.PointerParam, label="Experiment 1",
                      pointerClass='PKPDExperiment',
                      help='Select the experiment with the measurements you want to analyze.')
        form.addParam('X1', params.StringParam, label='X1', default='t', help='X for the plot in the 1st experiment.')
        form.addParam('Y1', params.StringParam, label='Y1', default='Cp', help='Y for the plot in the 1st experiment.')
        form.addParam('inputExperiment2', params.PointerParam, label="Experiment 2",
                      pointerClass='PKPDExperiment',
                      help='Select the experiment with the measurements you want to analyze.')
        form.addParam('X2', params.StringParam, label='X2', default='', help='X for the plot in the 2nd experiment. If empty, the same as X1.')
        form.addParam('Y2', params.StringParam, label='Y2', default='', help='Y for the plot in the 2nd experiment. If empty, the same as Y1.')

    #--------------------------- STEPS functions --------------------------------------------
    def _insertAllSteps(self):
        pass

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        retval=[]
        experiment1 = PKPDExperiment()
        experiment1.load(self.inputExperiment1.get().fnPKPD,fullRead=False)
        if not self.X1.get() in experiment1.variables:
            retval.append("Cannot find %s in Experiment1"%self.X1.get())
        if not self.Y1.get() in experiment1.variables:
            retval.append("Cannot find %s in Experiment1"%self.Y1.get())

        X2=self.X1.get() if self.X2.get()=="" else self.X2.get()
        Y2=self.Y1.get() if self.Y2.get()=="" else self.Y2.get()
        experiment2 = PKPDExperiment()
        experiment2.load(self.inputExperiment2.get().fnPKPD,fullRead=False)
        if not X2 in experiment2.variables:
            retval.append("Cannot find %s in Experiment2"%X2)
        if not Y2 in experiment2.variables:
            retval.append("Cannot find %s in Experiment2"%Y2)
        return retval