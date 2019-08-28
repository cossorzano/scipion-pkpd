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

import copy
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

import pyworkflow.protocol.params as params
from .protocol_pkpd import ProtPKPD
from pkpd.objects import PKPDExperiment, PKPDSample
from pkpd.utils import uniqueFloatValues

# Tested in test_workflow_dissolution_f2.py

class ProtPKPDAverageSample(ProtPKPD):
    """ Produce an experiment with a single sample whose value is the average of all the input samples.\n
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'average sample'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputExperiment', params.PointerParam, label="Input experiment",
                      pointerClass='PKPDExperiment',
                      help='Select an experiment with samples')
        form.addParam('resampleT', params.FloatParam, label="Resample profiles (time step)", default=-1,
                      help='Resample the input profiles at this time step (make sure it is in the same units as the input). '
                           'Leave it to -1 for no resampling. This is only valid when the label to compare is a measurement.')
        form.addParam('condition', params.StringParam, default="", label="Condition",
                      help='You must use available labels. Example: $(Oral_F0)<0.25 and $(sex)=="Female"')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runJoin',self.inputExperiment.get().getObjId(), self.resampleT.get())
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runJoin(self, objId, resampleT):
        experiment = self.readExperiment(self.inputExperiment.get().fnPKPD)
        tvarName = experiment.getTimeVariable()
        mvarNames = experiment.getMeasurementVariables()
        self.experiment = PKPDExperiment()

        # General
        self.experiment.general["title"]="Average of "+experiment.general["title"]
        if self.condition.get()!="":
            self.experiment.general["title"]+=" Condition: %s"%self.condition.get()
        self.experiment.general["comment"]=copy.copy(experiment.general["comment"])
        for key, value in experiment.general.iteritems():
            if not (key in self.experiment.general):
                self.experiment.general[key] = copy.copy(value)

        # Variables
        for key, value in experiment.variables.iteritems():
            if not (key in self.experiment.variables):
                self.experiment.variables[key] = copy.copy(value)

        # Vias
        for key, value in experiment.vias.iteritems():
            if not (key in self.experiment.vias):
                self.experiment.vias[key] = copy.copy(value)

        # Doses
        doseName = None
        for key, value in experiment.doses.iteritems():
            dose = copy.copy(value)
            self.experiment.doses[dose.doseName] = dose
            doseName = dose.doseName

        # Samples
        self.printSection("Averaging")
        self.experiment.samples["avg"] = PKPDSample()
        tokens=["AverageSample"]
        if doseName is not None:
            tokens.append("dose=%s"%doseName)
        self.experiment.samples["avg"].parseTokens(tokens, self.experiment.variables, self.experiment.doses,
                                                   self.experiment.groups)

        for mvarName in mvarNames:
            observations={}
            for sampleName, sample in experiment.getSubGroup(self.condition.get()).iteritems():
                print("%s participates in the average"%sampleName)
                t, y = sample.getXYValues(tvarName,mvarName)
                t=t[0] # [array]
                y=y[0] # [array]
                for i in range(len(t)):
                    ti = float(t[i])
                    if not ti in observations.keys():
                        observations[ti]=[]
                    observations[ti].append(y[i])
            for ti in observations:
                if len(observations[ti])>0:
                    observations[ti]=np.mean(observations[ti])
                else:
                    observations[ti]=np.nan

            t=[]
            yavg=[]
            for ti in sorted(observations):
                t.append(ti)
                yavg.append(observations[ti])

            if self.resampleT.get()>0:
                t,yavg=uniqueFloatValues(t,yavg)
                B = InterpolatedUnivariateSpline(t, yavg, k=1)
                t = np.arange(np.min(t), np.max(t) + self.resampleT.get(), self.resampleT.get())
                yavg = B(t)
            self.experiment.samples["avg"].addMeasurementColumn(tvarName,t)
            self.experiment.samples["avg"].addMeasurementColumn(mvarName,yavg)
        print(" ")

        # Print and save
        self.writeExperiment(self.experiment,self._getPath("experiment.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputExperiment=self.experiment)
        self._defineSourceRelation(self.inputExperiment, self.experiment)

