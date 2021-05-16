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

def NCA(ti,Cci):
    AUC = np.trapz(Cci, ti)
    indmax = np.argmax(Cci)
    tmax = ti[indmax]
    Cmax = Cci[indmax]
    return [tmax, Cmax, AUC]

class TestInhalation4Workflow(TestWorkflow):
    # Hartung2020_MATLAB/scripts/simulation_PKindexTable.m

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
        protSubstFP1= self.newProtocol(ProtPKPDInhSubstanceProperties,
                                     objLabel='pkpd - fluticasone propionate properties 1',
                                     name='fluticasone propionate',
                                     rho=1.43*500.57/1e3, MW=500.57,
                                     kdiss_alv=6.17e-5, kdiss_br=6.17e-5, kp_alv=92.6e-6*60, kp_br=92.6e-6*60,
                                     Cs_alv=11985*1e-3, Cs_br=11985*1e-3,
                                     Kpl_alv=2.47, Kpl_br=2.47, fu=0.0116, R=0.95
                                     )
        self.launchProtocol(protSubstFP1)
        self.assertIsNotNone(protSubstFP1.outputSubstanceParameters.fnSubst, "There was a problem with the substance definition")

        print("Substance (fluticasone propionate) ...")
        protSubstFP2= self.newProtocol(ProtPKPDInhSubstanceProperties,
                                     objLabel='pkpd - fluticasone propionate properties 2',
                                     name='fluticasone propionate',
                                     rho=1.43*500.57/1e3, MW=500.57,
                                     kdiss_alv=6.17e-5, kdiss_br=6.17e-5*0.2, kp_alv=92.6e-6*60, kp_br=92.6e-6*60,
                                     Cs_alv=11985*1e-3, Cs_br=11985*1e-3,
                                     Kpl_alv=2.47, Kpl_br=2.47, fu=0.0116, R=0.95
                                     )
        self.launchProtocol(protSubstFP2)
        self.assertIsNotNone(protSubstFP2.outputSubstanceParameters.fnSubst, "There was a problem with the substance definition")

        # Substance parameters
        print("Substance (budesonide) ...")
        protSubstBud1=self.newProtocol(ProtPKPDInhSubstanceProperties,
                                     objLabel='pkpd - budesonide properties 1',
                                     name='budesonide',
                                     rho=1.3, MW=430.53,
                                     kdiss_alv=3.3e-4, kdiss_br=3.3e-4, kp_alv=5.33e-6*60, kp_br=5.33e-6*60,
                                     Cs_alv=69.797, Cs_br=69.797,
                                     Kpl_alv=8, Kpl_br=8, fu=0.161, R=0.8
                                     )
        self.launchProtocol(protSubstBud1)
        self.assertIsNotNone(protSubstBud1.outputSubstanceParameters.fnSubst, "There was a problem with the substance definition")

        # Substance parameters
        print("Substance (budesonide) ...")
        protSubstBud2=self.newProtocol(ProtPKPDInhSubstanceProperties,
                                     objLabel='pkpd - budesonide properties 2',
                                     name='budesonide',
                                     rho=1.3, MW=430.53,
                                     kdiss_alv=3.3e-4, kdiss_br=3.3e-4*0.2, kp_alv=5.33e-6*60, kp_br=5.33e-6*60,
                                     Cs_alv=69.797, Cs_br=69.797,
                                     Kpl_alv=8, Kpl_br=8, fu=0.161, R=0.8
                                     )
        self.launchProtocol(protSubstBud2)
        self.assertIsNotNone(protSubstBud2.outputSubstanceParameters.fnSubst, "There was a problem with the substance definition")

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
        protDepoFPCA1550_1 = launchDepo('FP chamber asth 1,5um 50ug', 'FP_chamber_asthmatic_1,5um_50ug.txt', protSubstFP1)
        protDepoFPCA3050_1 = launchDepo('FP chamber asth 3,0um 50ug', 'FP_chamber_asthmatic_3,0um_50ug.txt', protSubstFP1)
        protDepoFPCA6050_1 = launchDepo('FP chamber asth 6,0um 50ug', 'FP_chamber_asthmatic_6,0um_50ug.txt', protSubstFP1)
        protDepoFPDH200_1  = launchDepo('FP diskus Healthy 200 ug', 'FP_Diskus_healthy_200ug.txt', protSubstFP1)
        protDepoFPDH500_1  = launchDepo('FP diskus Healthy 500 ug', 'FP_Diskus_healthy_500ug.txt', protSubstFP1)
        protDepoFPDH1000_1 = launchDepo('FP diskus Healthy 1000 ug', 'FP_Diskus_healthy_1000ug.txt', protSubstFP1)
        protDepoFPDA1000_1 = launchDepo('FP diskus Asthma 1000 ug', 'FP_Diskus_asthmatic_1000ug.txt', protSubstFP1)
        protDepoBudH400_1  = launchDepo('Bud Healthy 400 ug', 'Bud_Turbohaler_healthy_400ug.txt', protSubstBud1)
        protDepoBudH1000_1 = launchDepo('Bud Healthy 1000 ug', 'Bud_Turbohaler_healthy_1000ug.txt', protSubstBud1)
        protDepoBudH1200_1 = launchDepo('Bud Healthy 1200 ug', 'Bud_Turbohaler_healthy_1200ug.txt', protSubstBud1)
        protDepoBudA1200_1 = launchDepo('Bud Asthma 1200 ug', 'Bud_Turbohaler_asthmatic_1200ug.txt', protSubstBud1)

        protDepoFPCA1550_2 = launchDepo('FP chamber asth 1,5um 50ug', 'FP_chamber_asthmatic_1,5um_50ug.txt', protSubstFP2)
        protDepoFPCA3050_2 = launchDepo('FP chamber asth 3,0um 50ug', 'FP_chamber_asthmatic_3,0um_50ug.txt', protSubstFP2)
        protDepoFPCA6050_2 = launchDepo('FP chamber asth 6,0um 50ug', 'FP_chamber_asthmatic_6,0um_50ug.txt', protSubstFP2)
        protDepoFPDH200_2  = launchDepo('FP diskus Healthy 200 ug', 'FP_Diskus_healthy_200ug.txt', protSubstFP2)
        protDepoFPDH500_2  = launchDepo('FP diskus Healthy 500 ug', 'FP_Diskus_healthy_500ug.txt', protSubstFP2)
        protDepoFPDH1000_2 = launchDepo('FP diskus Healthy 1000 ug', 'FP_Diskus_healthy_1000ug.txt', protSubstFP2)
        protDepoFPDA1000_2 = launchDepo('FP diskus Asthma 1000 ug', 'FP_Diskus_asthmatic_1000ug.txt', protSubstFP2)
        protDepoBudH400_2  = launchDepo('Bud Healthy 400 ug', 'Bud_Turbohaler_healthy_400ug.txt', protSubstBud2)
        protDepoBudH1000_2 = launchDepo('Bud Healthy 1000 ug', 'Bud_Turbohaler_healthy_1000ug.txt', protSubstBud2)
        protDepoBudH1200_2 = launchDepo('Bud Healthy 1200 ug', 'Bud_Turbohaler_healthy_1200ug.txt', protSubstBud2)
        protDepoBudA1200_2 = launchDepo('Bud Asthma 1200 ug', 'Bud_Turbohaler_asthmatic_1200ug.txt', protSubstBud2)

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
        ka = 0.45; # [1 / min]

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
                                                                                                      Foral, ka))
        self.launchProtocol(protPKBud)
        self.assertIsNotNone(protPKBud.outputExperiment.fnPKPD, "There was a problem with the Bud PK parameters")

        def simulate(protDepo, subst, label, tmax0, Cmax0, AUC0, MCC0):
            if subst=="FP":
                simulationTime = 12*60
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
            MCC = float(experiment.samples['simulation'].getDescriptorValue("mcc_cleared_lung_dose_fraction"))
            [tmax, Cmax, AUC] = NCA(t,Cnmol)
            print([tmax, Cmax, AUC, MCC])
            self.assertTrue(abs(tmax - tmax0) < 0.01*tmax0)
            self.assertTrue(abs(Cmax - Cmax0) < 0.01*Cmax0)
            self.assertTrue(abs(AUC - AUC0) < 0.01*AUC0)
            self.assertTrue(abs(MCC - MCC0) < 0.01*MCC0)

        simulate(protDepoFPCA1550_1,  'FP', 'FP chamber asth  1,5um 50ug (1)', 35.64, 1.9472e-4, 0.0285, 0.10489)
        simulate(protDepoFPCA3050_1,  'FP', 'FP chamber asth  3,0um 50ug (1)', 65.88, 5.2204e-05, 0.0216, 0.1796)
        simulate(protDepoFPCA6050_1,  'FP', 'FP chamber asth  6,0um 50ug (1)', 278.64, 0.1281e-4, 0.0088, 0.3477)
        simulate(protDepoFPDH200_1,   'FP', 'FP diskus Healthy 200 ug (1)', 41.04, 0.8315e-4, 0.0233, 0.2445)
        simulate(protDepoFPDH500_1,   'FP', 'FP diskus Healthy 500 ug (1)', 41.04, 0.2074e-3, 0.0581, 0.2449)
        simulate(protDepoFPDH1000_1,  'FP', 'FP diskus Healthy 1000 ug (1)', 41.04, 0.4133e-3, 0.1155, 0.2456)
        simulate(protDepoFPDA1000_1,  'FP', 'FP diskus Asthma 1000 ug (1)', 37.44, 0.3359e-3, 0.0823, 0.4699)
        simulate(protDepoBudH400_1,  'Bud', 'Bud Healthy 400 ug (1)', 49.392, 0.9162e-3, 0.2546, 0.0269)
        simulate(protDepoBudH1000_1, 'Bud', 'Bud Healthy 1000 ug (1)', 49.176, 0.00227909, 0.6338, 0.0316)
        simulate(protDepoBudH1200_1, 'Bud', 'Bud Healthy 1200 ug (1)', 49.104, 2.73152e-3, 0.7594, 0.0329)
        simulate(protDepoBudA1200_1, 'Bud', 'Bud Asthma 1200 ug (1)', 49.68, 0.001859812, 0.5662158, 0.2111286)

        simulate(protDepoFPCA1550_2,  'FP', 'FP chamber asth  1,5um 50ug (2)', 36.72, 1.4470e-04, 0.0237, 0.2008)
        simulate(protDepoFPCA3050_2,  'FP', 'FP chamber asth  3,0um 50ug (2)', 69.84, 3.4221e-05, 0.0143, 0.2834)
        simulate(protDepoFPCA6050_2,  'FP', 'FP chamber asth  6,0um 50ug (2)', 404.28, 0.0698e-4, 0.0048, 0.4028)
        simulate(protDepoFPDH200_2,   'FP', 'FP diskus Healthy 200 ug (2)', 41.40, 0.7590e-4, 0.0214, 0.2751)
        simulate(protDepoFPDH500_2,   'FP', 'FP diskus Healthy 500 ug (2)', 41.40, 0.1897e-3, 0.0533, 0.2751)
        simulate(protDepoFPDH1000_2,  'FP', 'FP diskus Healthy 1000 ug (2)', 41.04, 0.3794e-3, 0.1065, 0.2752)
        simulate(protDepoFPDA1000_2,  'FP', 'FP diskus Asthma 1000 ug (2)', 39.60, 0.2557e-3, 0.0651, 0.5343)
        simulate(protDepoBudH400_2,  'Bud', 'Bud Healthy 400 ug (2)', 48.348, 0.8997e-3, 0.2487, 0.0546)
        simulate(protDepoBudH1000_2, 'Bud', 'Bud Healthy 1000 ug (2)', 48.276, 0.00224676, 0.6199, 0.0572)
        simulate(protDepoBudH1200_2, 'Bud', 'Bud Healthy 1200 ug (2)', 48.24, 2.695265e-3, 0.7432, 0.0580)
        simulate(protDepoBudA1200_2, 'Bud', 'Bud Asthma 1200 ug (2)', 45.90, 0.00178895, 0.520985, 0.288722)

if __name__ == "__main__":
    unittest.main()
