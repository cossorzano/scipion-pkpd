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

class TestInhalation3Workflow(TestWorkflow):
    # Hartung2020_MATLAB/scripts/simulation_PKstudiesHealthy.m

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
        print("Substance (fluticasone propionate) ...")
        protSubstFP= self.newProtocol(ProtPKPDInhSubstanceProperties,
                                     objLabel='pkpd - fluticasone propionate properties',
                                     name='fluticasone propionate',
                                     rho=1.43*500.57/1e3, MW=500.57,
                                     kdiss_alv=6.17e-5, kdiss_br=6.17e-5*0.2, kp_alv=92.6e-6*60, kp_br=92.6e-6*60,
                                     Cs_alv=11985*1e-3, Cs_br=11985*1e-3,
                                     Kpl_alv=2.47, Kpl_br=2.47, fu=0.0116, R=0.95
                                     )
        self.launchProtocol(protSubstFP)
        self.assertIsNotNone(protSubstFP.outputSubstanceParameters.fnSubst, "There was a problem with the substance definition")

        # Substance parameters
        print("Substance (budesonide) ...")
        protSubstBud=self.newProtocol(ProtPKPDInhSubstanceProperties,
                                     objLabel='pkpd - budesonide properties',
                                     name='budesonide',
                                     rho=1.3e9/430.53, MW=430.53,
                                     kdiss_alv=3.3e-4, kdiss_br=3.3e-4*0.2, kp_alv=5.33e-6*60, kp_br=5.33e-6*60,
                                     Cs_alv=69.797, Cs_br=69.797,
                                     Kpl_alv=8, Kpl_br=8, fu=0.161, R=0.8
                                     )
        self.launchProtocol(protSubstBud)
        self.assertIsNotNone(protSubstBud.outputSubstanceParameters.fnSubst, "There was a problem with the substance definition")

        def launchDepo(label, fn, protSubst):
            print("Deposition %s ..."%label)
            protDepo = self.newProtocol(ProtPKPDInhImportDepositionProperties,
                                        objLabel='pkpd - deposition %s'%label,
                                        depositionFile=self.dataset.getFile(fn))
            protDepo.substance.set(protSubst.outputSubstanceParameters)
            protDepo.lungModel.set(protLung.outputLungParameters)
            self.launchProtocol(protDepo)
            self.assertIsNotNone(protDepo.outputDeposition.fnDeposition, "There was a problem with the deposition 15")
            return protDepo

        # Deposition parameters
        protDepoFPH200 = launchDepo('FP Healthy 200 ug', 'FP_Diskus_healthy_200ug.txt', protSubstFP)
        protDepoFPH500 = launchDepo('FP Healthy 500 ug', 'FP_Diskus_healthy_500ug.txt', protSubstFP)
        protDepoFPH1000 = launchDepo('FP Healthy 1000 ug', 'FP_Diskus_healthy_1000ug.txt', protSubstFP)
        protDepoFPA1000 = launchDepo('FP Asthma 1000 ug', 'FP_Diskus_asthmatic_1000ug.txt', protSubstFP)
        protDepoBudH400 = launchDepo('Bud Healthy 400 ug', 'Bud_Turbohaler_healthy_400ug.txt', protSubstBud)
        protDepoBudH1000 = launchDepo('Bud Healthy 1000 ug', 'Bud_Turbohaler_healthy_1000ug.txt', protSubstBud)
        protDepoBudH1200 = launchDepo('Bud Healthy 1200 ug', 'Bud_Turbohaler_healthy_1200ug.txt', protSubstBud)
        protDepoBudA1200 = launchDepo('Bud Asthma 1200 ug', 'Bud_Turbohaler_asthmatic_1200ug.txt', protSubstBud)

        # PK parameters
        print("FP PK parameters ...")
        CL_L_h  = 73;          #[L/h]
        Vc_L    = 31;          #[L]
        k12_1_h = 1.78;        #[1/h]
        k21_1_h = 0.09;        #[1/h]
        Q_L_h   = k12_1_h * Vc_L;
        Vp_L    = Q_L_h/k21_1_h;

        CL_mL_min = 1000 / 60 * CL_L_h; # [L / h] --> [mL / min]
        Q_mL_min = 1000 / 60 * Q_L_h; # [L / h] --> [mL / min]

        Vc_mL = 1000 * Vc_L; # [L] --> [mL]
        Vp_mL = 1000 * Vp_L; # [L] --> [mL]

        protPKFP = self.newProtocol(ProtPKPDCreateExperiment,
                                    objLabel='pkpd - FP pk parameters',
                                    newTitle='Fluticasone propionate PK parameters',
                                    newVariables='Cl ; mL/min ; numeric[%f] ; label ; Two compartments, central clearance\n'
                                                 'V ; mL ; numeric[%f] ; label ; Two compartments, central volume\n'
                                                 'Vp ; mL ; numeric[%f] ; label ; Two compartments, peripheral volume\n'
                                                 'Q ; mL/min ; numeric[%f] ; label ; Two compartments, passage rate from central to peripheral and viceversa\n'
                                                 'F ; none ; numeric[%f] ; label ; Fraction that is absorbed orally\n'
                                                 'k01 ; 1/min ; numeric[%f] ; label ; 1st order absorption rate of the oral fraction\n',
                                    newSamples='Individual1; Cl=%f; V=%f; Vp=%f; Q=%f; F=0; k01=0'%(CL_mL_min,Vc_mL,
                                                                                                  Vp_mL, Q_mL_min))
        self.launchProtocol(protPKFP)
        self.assertIsNotNone(protPKFP.outputExperiment.fnPKPD, "There was a problem with the FP PK parameters")

        print("Bud PK parameters ...")
        CL_L_h = 85; # [L / h]
        Vc_L = 100; # [L]
        k12_1_h = 20.01; # [1 / h]
        k21_1_h = 11.06; # [1 / h]
        Q_L_h = k12_1_h * Vc_L;
        Vp_L = Q_L_h / k21_1_h;
        Foral = 0.11; # []
        ka = 0.45; # [1 / h]

        CL_mL_min = 1000 / 60 * CL_L_h; # [L / h] --> [mL / min]
        Q_mL_min = 1000 / 60 * Q_L_h; # [L / h] --> [mL / min]

        Vc_mL = 1000 * Vc_L; # [L] --> [mL]
        Vp_mL = 1000 * Vp_L; # [L] --> [mL]

        protPKBud = self.newProtocol(ProtPKPDCreateExperiment,
                                    objLabel='pkpd - Bud pk parameters',
                                    newTitle='Budesonide PK parameters',
                                    newVariables='Cl ; mL/min ; numeric[%f] ; label ; Two compartments, central clearance\n'
                                                 'V ; mL ; numeric[%f] ; label ; Two compartments, central volume\n'
                                                 'Vp ; mL ; numeric[%f] ; label ; Two compartments, peripheral volume\n'
                                                 'Q ; mL/min ; numeric[%f] ; label ; Two compartments, passage rate from central to peripheral and viceversa\n'
                                                 'F ; none ; numeric[%f] ; label ; Fraction that is absorbed orally\n'
                                                 'k01 ; 1/min ; numeric[%f] ; label ; 1st order absorption rate of the oral fraction\n',
                                    newSamples='Individual1; Cl=%f; V=%f; Vp=%f; Q=%f; F=%f; k01=%f'%(CL_mL_min,Vc_mL,
                                                                                                      Vp_mL, Q_mL_min,
                                                                                                      Foral, ka/60))
        self.launchProtocol(protPKBud)
        self.assertIsNotNone(protPKBud.outputExperiment.fnPKPD, "There was a problem with the Bud PK parameters")

        def simulate(protDepo, subst, label):
            if subst=="FP":
                simulationTime = 24*60
                deltaT = simulationTime/2000
                protPK = protPKFP
            else:
                simulationTime = 12*60
                deltaT = simulationTime/20000
                protPK = protPKBud

            print("Inhalation simulation %s ..."%label)
            protSimulate = self.newProtocol(ProtPKPDInhSimulate,
                                            objLabel='pkpd - simulate inhalation %s'%label,
                                            diameters="0.1,1.1,0.1; 1.2,24.2,0.2",
                                            simulationTime=simulationTime,
                                            deltaT=deltaT)
            protSimulate.ptrDeposition.set(protDepo.outputDeposition)
            protSimulate.ptrPK.set(protPK.outputExperiment)
            self.launchProtocol(protSimulate)
            self.assertIsNotNone(protSimulate.outputExperiment.fnPKPD, "There was a problem with the simulation")
            experiment = PKPDExperiment()
            experiment.load(protSimulate.outputExperiment.fnPKPD)
            t = np.asarray([float(x) for x in experiment.samples['simulation'].getValues('t')])
            Cnmol = np.asarray([float(x) for x in experiment.samples['simulation'].getValues('Cnmol')])
        simulate(protDepoFPH200, 'FP', 'FP healthy 200')


if __name__ == "__main__":
    unittest.main()
