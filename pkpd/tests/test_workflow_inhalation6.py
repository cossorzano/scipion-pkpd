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

class TestInhalation6Workflow(TestWorkflow):
    # Hartung2020_MATLAB/scripts/simulation_sensitivityAnalysis_FP.m

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
                                     objLabel='pkpd - fluticasone propionate properties 2',
                                     name='fluticasone propionate',
                                     rho=1.43*500.57/1e3, MW=500.57,
                                     kdiss_alv=6.17e-5, kdiss_br=6.17e-5*0.2, kp_alv=92.6e-6*60, kp_br=92.6e-6*60,
                                     Cs_alv=11985*1e-3, Cs_br=11985*1e-3,
                                     Kpl_alv=2.47, Kpl_br=2.47, fu=0.0116, R=0.95
                                     )
        self.launchProtocol(protSubstFP)
        self.assertIsNotNone(protSubstFP.outputSubstanceParameters.fnSubst, "There was a problem with the substance definition")


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
        protDepoFPDA250_x1 = launchDepo('FP diskus Asthma 250 ug', 'FP_Diskus_asthmatic_250ug.txt', protSubstFP)
        protDepoFPDA250_x2 = launchDepo('FP diskus Asthma 250 ug x2', 'FP_Diskus_asthmatic_250ug_2FoldLarger.txt', protSubstFP)
        protDepoFPDA250_x05 = launchDepo('FP diskus Asthma 250 ug x0.5', 'FP_Diskus_asthmatic_250ug_2FoldSmaller.txt', protSubstFP)

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

        def simulate(protDepo, subst, label, tmax0, Cmax0, AUC0):
            if subst=="FP":
                simulationTime = 24*60
                deltaT = simulationTime/2000
                protPK = protPKFP

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
            print(np.mean(Cnmol))
            # [tmax, Cmax, AUC] = NCA(t,Cnmol)
            # print([tmax, Cmax, AUC])
            # self.assertTrue(abs(tmax - tmax0) < 0.001*tmax0)
            # self.assertTrue(abs(Cmax - Cmax0) < 0.001*Cmax0)
            # self.assertTrue(abs(AUC - AUC0) < 0.001*AUC0)

        simulate(protDepoFPDA250_x1,  'FP', 'Reference', 39.60, 0.25556811e-3, 0.08235074)
        # simulate(protDepoFPDH1000_2,  'FP', 'FP diskus Healthy 1000 ug (2)', 41.76, 0.379130786e-3, 0.1347830376)

if __name__ == "__main__":
    unittest.main()
