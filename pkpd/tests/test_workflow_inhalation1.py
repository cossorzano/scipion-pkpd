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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import csv
import numpy as np
import os
import pandas

from pyworkflow.tests import *
from pkpd.protocols import *
from pkpd.objects import PKPDDataSet
from .test_workflow import TestWorkflow

class TestInhalation1Workflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.dataset = PKPDDataSet.getDataSet('Inhalation')

    def testDissolutionWorkflow(self):
        # Lung parameters
        print("Define lung parameters")
        protLung = self.newProtocol(ProtPKPDInhLungPhysiology,
                                      objLabel='pkpd - lung parameters')
        self.launchProtocol(protLung)
        self.assertIsNotNone(protLung.outputLungParameters.fnPhys, "There was a problem with the lung definition")

        # Substance parameters
        print("Substance (gold) ...")
        protGold= self.newProtocol(ProtPKPDInhSubstanceProperties,
                                   objLabel='pkpd - gold properties',
                                   name='gold',
                                   rho=19.3, MW=197,
                                   kdiss_alv=0, kdiss_br=0, kp_alv=0, kp_br=0, Cs_alv=1, Cs_br=1,
                                   Kpl_alv=1, Kpl_br=1, fu=1, R=0
                                   )
        self.launchProtocol(protGold)
        self.assertIsNotNone(protGold.outputSubstanceParameters.fnSubst, "There was a problem with the gold definition")

        # Deposition parameters
        print("Deposition ...")
        protDepo = self.newProtocol(ProtPKPDInhImportDepositionProperties,
                                    objLabel='pkpd - deposition',
                                    depositionFile=self.dataset.getFile('deposition1'))
        protDepo.substance.set(protGold.outputSubstanceParameters)
        protDepo.lungModel.set(protLung.outputLungParameters)
        self.launchProtocol(protDepo)
        self.assertIsNotNone(protDepo.outputDeposition.fnDeposition, "There was a problem with the deposition")

        # PK parameters
        print("PK parameters ...")
        protPK = self.newProtocol(ProtPKPDCreateExperiment,
                                  objLabel='pkpd - pk parameters',
                                  newTitle='Gold PK parameters',
                                  newVariables='Cl ; mL/min ; numeric[%f] ; label ; Two compartments, central clearance\n'
                                               'V ; mL ; numeric[%f] ; label ; Two compartments, central volume\n'
                                               'Vp ; mL ; numeric[%f] ; label ; Two compartments, peripheral volume\n'
                                               'Q ; mL/min ; numeric[%f] ; label ; Two compartments, passage rate from central to peripheral and viceversa\n'
                                               'F ; none ; numeric[%f] ; label ; Fraction that is absorbed orally\n'
                                               'k01 ; 1/min ; numeric[%f] ; label ; 1st order absorption rate of the oral fraction\n',
                                  newSamples='Individual1; Cl=0; V=1000; Vp=1000; Q=0; F=0; k01=0')
        self.launchProtocol(protPK)
        self.assertIsNotNone(protPK.outputExperiment.fnPKPD, "There was a problem with the PK parameters")

        # Simulate inhalation
        print("Inhalation simulation ...")
        protSimulate = self.newProtocol(ProtPKPDInhSimulate,
                                        objLabel='pkpd - simulate inhalation',
                                        deltaT=1.8)
        protSimulate.ptrDeposition.set(protDepo.outputDeposition)
        protSimulate.ptrPK.set(protPK.outputExperiment)
        self.launchProtocol(protSimulate)
        self.assertIsNotNone(protSimulate.outputExperiment.fnPKPD, "There was a problem with the simulation")
        experiment = PKPDExperiment()
        experiment.load(protSimulate.outputExperiment.fnPKPD)
        t = np.asarray([float (x) for x in experiment.samples['simulation'].getValues('t')])
        retention = np.asarray([float(x) for x in experiment.samples['simulation'].getValues('Retention')])
        self.assertTrue(abs(retention[0]-100.0)<0.001)
        self.assertTrue(abs(retention[8000]-1.3809)<0.001)

        # Plot short term
        dataSmith=pandas.read_csv(self.dataset.getFile('SmithPSLGold6'))
        IDs=pandas.unique(dataSmith['Subject'])

        import matplotlib.pyplot as plt
        plt.figure(figsize=(12,9))
        plt.plot(t/60,retention)
        plt.title('Short-term retention')
        plt.xlabel('Time (hours)')
        plt.ylabel('Lung retention (% lung dose)')
        plt.ylim([-2, 102])
        plt.xlim([0, 24])
        plt.xticks(np.arange(0,24,4))

        tmax_splot1 = 24 * 60  # 24 h
        colors = ['r','g','b','c','y','m','k']
        i=0
        legends=['PDE model prediction']
        for subject in IDs:
            data_i = dataSmith[dataSmith.Subject==subject]
            PSL_i = data_i[data_i.Particle=='PSL']
            Gold_i = data_i[data_i.Particle=='Gold']

            indPSL_splot1 = PSL_i.Time_h <= tmax_splot1 / 60
            indGld_splot1 = Gold_i.Time_h <= tmax_splot1 / 60

            plt.plot(PSL_i.Time_h[indPSL_splot1], PSL_i.lungRetention[indPSL_splot1], "%sv"%colors[i])
            plt.plot(Gold_i.Time_h[indGld_splot1], Gold_i.lungRetention[indGld_splot1], "%sx"%colors[i])
            legends.append('PSL %s'%subject)
            legends.append('Gold %s'%subject)
            i+=1
        plt.legend(legends)
        plt.savefig('shortTerm.png')

        # Plot long term
        plt.figure(figsize=(12, 9))
        plt.plot(t / (60*24), retention)
        plt.title('Long-term retention')
        plt.xlabel('Time (days)')
        plt.ylabel('Lung retention (% lung dose)')
        plt.ylim([-2, 102])
        plt.xlim([0, 10.5])

        i = 0
        legends = ['PDE model prediction']
        for subject in IDs:
            data_i = dataSmith[dataSmith.Subject == subject]
            PSL_i = data_i[data_i.Particle == 'PSL']
            Gold_i = data_i[data_i.Particle == 'Gold']

            plt.plot(PSL_i.Time_h/24, PSL_i.lungRetention, "%sv" % colors[i])
            plt.plot(Gold_i.Time_h/24, Gold_i.lungRetention, "%sx" % colors[i])
            legends.append('PSL %s' % subject)
            legends.append('Gold %s' % subject)
            i += 1
        plt.legend(legends)
        plt.savefig('longTerm.png')

if __name__ == "__main__":
    unittest.main()
