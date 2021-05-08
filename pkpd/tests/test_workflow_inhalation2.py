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

from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

from pyworkflow.tests import *
from pkpd.protocols import *
from pkpd.objects import PKPDDataSet
from .test_workflow import TestWorkflow

def NCA(ti,Cci,MW):
    AUC_vec = cumtrapz(Cci, ti, initial=0)
    interpolator = interp1d(ti, AUC_vec);
    AUC_12h = interpolator(12*60)
    AUC_Sys_12h_scaled = AUC_12h * MW * 1000 / 60; # [nmol * min / mL] * [g / mol] * 1000 / 60 = [pg * h / mL]
    indmax = np.argmax(Cci)
    tmax = ti[indmax]
    Cmax = Cci[indmax] * MW * 1000; # [nmol/mL]*[g/mol]*1000 = [pg/mL]
    return [tmax, Cmax, AUC_Sys_12h_scaled]

class TestInhalation2Workflow(TestWorkflow):
    # Hartung2020_MATLAB/scripts/simulation_particleSizeFP.m

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
        MW = 500.57
        protSubst= self.newProtocol(ProtPKPDInhSubstanceProperties,
                                   objLabel='pkpd - fluticasone propionate properties',
                                   name='fluticasone propionate',
                                   rho=1.43*500.57/1e3, MW=MW,
                                   kdiss_alv=6.17e-5, kdiss_br=6.17e-5*0.2, kp_alv=92.6e-6*60, kp_br=92.6e-6*60,
                                   Cs_alv=11985*1e-3, Cs_br=11985*1e-3,
                                   Kpl_alv=2.47, Kpl_br=2.47, fu=0.0116, R=0.95
                                   )
        self.launchProtocol(protSubst)
        self.assertIsNotNone(protSubst.outputSubstanceParameters.fnSubst, "There was a problem with the substance definition")

        # Deposition1 parameters
        print("Deposition15 ...")
        protDepo15 = self.newProtocol(ProtPKPDInhImportDepositionProperties,
                                      objLabel='pkpd - deposition 1.5',
                                      depositionFile=self.dataset.getFile('FP_chamber_asthmatic_1.5um_50ug'))
        protDepo15.substance.set(protSubst.outputSubstanceParameters)
        protDepo15.lungModel.set(protLung.outputLungParameters)
        self.launchProtocol(protDepo15)
        self.assertIsNotNone(protDepo15.outputDeposition.fnDeposition, "There was a problem with the deposition 15")

        print("Deposition30 ...")
        protDepo30 = self.newProtocol(ProtPKPDInhImportDepositionProperties,
                                      objLabel='pkpd - deposition 3.0',
                                      depositionFile=self.dataset.getFile('FP_chamber_asthmatic_3.0um_50ug'))
        protDepo30.substance.set(protSubst.outputSubstanceParameters)
        protDepo30.lungModel.set(protLung.outputLungParameters)
        self.launchProtocol(protDepo30)
        self.assertIsNotNone(protDepo30.outputDeposition.fnDeposition, "There was a problem with the deposition 30")

        print("Deposition60 ...")
        protDepo60 = self.newProtocol(ProtPKPDInhImportDepositionProperties,
                                      objLabel='pkpd - deposition 6.0',
                                      depositionFile=self.dataset.getFile('FP_chamber_asthmatic_6.0um_50ug'))
        protDepo60.substance.set(protSubst.outputSubstanceParameters)
        protDepo60.lungModel.set(protLung.outputLungParameters)
        self.launchProtocol(protDepo60)
        self.assertIsNotNone(protDepo60.outputDeposition.fnDeposition, "There was a problem with the deposition 60")

        # PK parameters
        print("PK parameters ...")
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

        protPK = self.newProtocol(ProtPKPDCreateExperiment,
                                  objLabel='pkpd - pk parameters',
                                  newTitle='Fluticasone propionate PK parameters',
                                  newVariables='Cl ; mL/min ; numeric[%f] ; label ; Two compartments, central clearance\n'
                                               'V ; mL ; numeric[%f] ; label ; Two compartments, central volume\n'
                                               'Vp ; mL ; numeric[%f] ; label ; Two compartments, peripheral volume\n'
                                               'Q ; mL/min ; numeric[%f] ; label ; Two compartments, passage rate from central to peripheral and viceversa\n'
                                               'F ; none ; numeric[%f] ; label ; Fraction that is absorbed orally\n'
                                               'k01 ; 1/min ; numeric[%f] ; label ; 1st order absorption rate of the oral fraction\n',
                                  newSamples='Individual1; Cl=%f; V=%f; Vp=%f; Q=%f; F=0; k01=0'%(CL_mL_min,Vc_mL,
                                                                                                  Vp_mL, Q_mL_min))
        self.launchProtocol(protPK)
        self.assertIsNotNone(protPK.outputExperiment.fnPKPD, "There was a problem with the PK parameters")

        def simulate(protDepo, label, tmax0, Cmax0, AUC0):
            print("Inhalation simulation %s ..."%label)
            protSimulate = self.newProtocol(ProtPKPDInhSimulate,
                                                objLabel='pkpd - simulate inhalation %s'%label,
                                                simulationTime=12 * 60,
                                                deltaT=0.18)
            protSimulate.ptrDeposition.set(protDepo.outputDeposition)
            protSimulate.ptrPK.set(protPK.outputExperiment)
            self.launchProtocol(protSimulate)
            self.assertIsNotNone(protSimulate.outputExperiment.fnPKPD, "There was a problem with the simulation")
            experiment = PKPDExperiment()
            experiment.load(protSimulate.outputExperiment.fnPKPD)
            t = np.asarray([float(x) for x in experiment.samples['simulation'].getValues('t')])
            Cnmol = np.asarray([float(x) for x in experiment.samples['simulation'].getValues('Cnmol')])
            [tmax, Cmax, AUC_Sys_12h_scaled] = NCA(t, Cnmol, MW)
            print([tmax, Cmax, AUC_Sys_12h_scaled])
            self.assertTrue(abs(tmax - tmax0) < 0.0001)
            self.assertTrue(abs(Cmax - Cmax0) < 0.0001)
            self.assertTrue(abs(AUC_Sys_12h_scaled - AUC0) < 0.0001)
            return [tmax, Cmax, AUC_Sys_12h_scaled]

        # Simulate inhalations
        [tmax15x1, Cmax15x1, AUC_Sys_12h_scaled15x1] = simulate(protDepo15, "1.5x1",  36.36, 72.4503, 197.5372)
        [tmax30x1, Cmax30x1, AUC_Sys_12h_scaled30x1] = simulate(protDepo30, "3.0x1",  69.30, 17.1307, 119.5475)
        [tmax60x1, Cmax60x1, AUC_Sys_12h_scaled60x1] = simulate(protDepo60, "6.0x1", 404.28,  3.4926,  40.1811)

        print([tmax15x1, Cmax15x1, AUC_Sys_12h_scaled15x1])
        print([tmax30x1, Cmax30x1, AUC_Sys_12h_scaled30x1])
        print([tmax60x1, Cmax60x1, AUC_Sys_12h_scaled60x1])

if __name__ == "__main__":
    unittest.main()
