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

        AUC_brtis = []
        Ct_brtis = []
        AUC_sys = []
        Cmax_sys = []

        def simulate(protDepo, subst, label, Cmax0, AUC0, Ct_brtisend0, AUCbr0,
                     substanceMultiplier="1 1 1 1 1 1 1 1",
                     physiologyMultiplier="1 1 1 1 1 1 1 1 1",
                     pkMultiplier="1 1 1 1 1 1",
                     volMultiplier=1,
                     ciliarySpeedType=2,
                     diameters="0.1,1.1,0.1; 1.2,24.2,0.2"):
            if subst=="FP":
                simulationTime = 24*60
                deltaT = simulationTime/2000
                protPK = protPKFP

            print("Inhalation simulation %s ..."%label)
            protSimulate = self.newProtocol(ProtPKPDInhSimulate,
                                            objLabel='pkpd - simulate inhalation %s'%label,
                                            diameters=diameters,
                                            simulationTime=simulationTime,
                                            substanceMultiplier=substanceMultiplier,
                                            physiologyMultiplier=physiologyMultiplier,
                                            pkMultiplier=pkMultiplier,
                                            ciliarySpeedType=ciliarySpeedType,
                                            volMultiplier=volMultiplier,
                                            deltaT=deltaT)
            protSimulate.ptrDeposition.set(protDepo.outputDeposition)
            protSimulate.ptrPK.set(protPK.outputExperiment)
            self.launchProtocol(protSimulate)
            self.assertIsNotNone(protSimulate.outputExperiment.fnPKPD, "There was a problem with the simulation")
            experiment = PKPDExperiment()
            experiment.load(protSimulate.outputExperiment.fnPKPD)
            t = np.asarray([float(x) for x in experiment.samples['simulation'].getValues('t')])
            Cnmol = np.asarray([float(x) for x in experiment.samples['simulation'].getValues('Cnmol')])
            CbrTis = np.asarray([float(x) for x in experiment.samples['simulation'].getValues('CbrTis')])
            [_, _, AUCbr] = NCA(t,CbrTis)
            Ct_brtisend = CbrTis[-1]
            [_, Cmax, AUC] = NCA(t,Cnmol)
            print([Cmax, AUC, Ct_brtisend, AUCbr])
            self.assertTrue(abs(Cmax - Cmax0) < 0.001*Cmax0)
            self.assertTrue(abs(AUC - AUC0) < 0.001*AUC0)
            self.assertTrue(abs(Ct_brtisend - Ct_brtisend0) < 0.001*Ct_brtisend0)
            self.assertTrue(abs(AUCbr - AUCbr0) < 0.001*AUCbr0)

            AUC_brtis.append(AUCbr)
            Ct_brtis.append(Ct_brtisend)
            AUC_sys.append(AUC)
            Cmax_sys.append(Cmax)

        simulate(protDepoFPDA250_x1,  'FP', 'Reference', 0.000064061396015, 0.020939149041018, 0.000047763114995, 0.201099418846000)
        simulate(protDepoFPDA250_x1,  'FP', 'Solubility br x2', 0.000064090773860, 0.021001912431073, 0.000047340079904, 0.203125096139258,
                 substanceMultiplier="2 1 1 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Solubility br x0.5', 0.000064003680934, 0.020816536963454, 0.000048382688314, 0.197114852202011,
                 substanceMultiplier="0.5 1 1 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Solubility alv x2', 0.000064061433533, 0.020939149234721, 0.000047763114528, 0.201099419350972,
                 substanceMultiplier="1 2 1 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Solubility alv x0.5', 0.000064061320978, 0.020939148653611, 0.000047763115930, 0.201099417836058,
                 substanceMultiplier="1 0.5 1 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Dissolution rate br x2', 0.000070148629113, 0.023161888978793, 0.000042478479307, 0.270857982092349,
                 substanceMultiplier="1 1 2 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Dissolution rate br x0.5', 0.000060926178857, 0.019112370601936, 0.000042405011823, 0.141402354488804,
                 substanceMultiplier="1 1 0.5 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Dissolution rate alv x2', 0.000104435140614, 0.021167186535426, 0.000047202018843, 0.201693928634031,
                 substanceMultiplier="1 1 1 2 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Dissolution rate alv x0.5', 0.000038324095550, 0.020205292617607, 0.000049918823834, 0.199185196481824,
                 substanceMultiplier="1 1 1 0.5 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Permeability br x2', 0.000064714794121, 0.021237282325590, 0.000043666859564, 0.210101389923662,
                 substanceMultiplier="1 1 1 1 2 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Permeability br x0.5', 0.000063423760931, 0.020510505894042, 0.000049668053713, 0.187456517185514,
                 substanceMultiplier="1 1 1 1 0.5 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Permeability alv x2', 0.000064061465104, 0.020939155418676, 0.000047763100214, 0.201099435470356,
                 substanceMultiplier="1 1 1 1 1 2 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Permeability alv x0.5', 0.000064061256795, 0.020939136285631, 0.000047763144557, 0.201099385597105,
                 substanceMultiplier="1 1 1 1 1 0.5 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Partition coefficient br x2', 0.000063757685846, 0.020931764321828, 0.000095845591570, 0.401884069594298,
                 substanceMultiplier="1 1 1 1 1 1 2 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Partition coefficient br x0.5', 0.000064159478182, 0.020942830269960, 0.000023841842712, 0.100588868560786,
                 substanceMultiplier="1 1 1 1 1 1 0.5 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Partition coefficient alv x2', 0.000063683703509, 0.020933726571695, 0.000047771314560, 0.201085296873244,
                 substanceMultiplier="1 1 1 1 1 1 1 2")
        simulate(protDepoFPDA250_x1,  'FP', 'Partition coefficient alv x0.5', 0.000064247463985, 0.020941858364839, 0.000047759016417, 0.201106474860665,
                 substanceMultiplier="1 1 1 1 1 1 1 0.5")
        simulate(protDepoFPDA250_x1,  'FP', 'Systemic clearance x2', 0.000043673632008, 0.011538065829228, 0.000040622785497, 0.176677121918005,
                 pkMultiplier="2 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Systemic clearance x0.5', 0.000084136170595, 0.034603184700807, 0.000063271700881, 0.236581349709199,
                 pkMultiplier="0.5 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Mucociliary clearance x2', 0.000063421544571, 0.019851111548667, 0.000038971635122, 0.164585396158581,
                 physiologyMultiplier="2 1 1 1 1 1 1 1 1", ciliarySpeedType=1)
        simulate(protDepoFPDA250_x1,  'FP', 'Mucociliary clearance x0.5', 0.000065861537692, 0.022147849161159, 0.000060279674836, 0.241230162147393,
                 physiologyMultiplier="0.5 1 1 1 1 1 1 1 1", ciliarySpeedType=1)
        simulate(protDepoFPDA250_x1,  'FP', 'Fluid volume br x2', 0.000063463707336, 0.020595106538030, 0.000049657196576, 0.190262165108595,
                 physiologyMultiplier="1 2 1 1 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Fluid volume br x0.5', 0.000064669808631, 0.021157512298910, 0.000044666634208, 0.207609943505915,
                 physiologyMultiplier="1 0.5 1 1 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Fluid volume alv x2', 0.000064061162448, 0.020939134289553, 0.000047763147251, 0.201099380399566,
                 physiologyMultiplier="1 1 2 1 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Fluid volume alv x0.5', 0.000064061512194, 0.020939156416699, 0.000047763098867, 0.201099438069086,
                 physiologyMultiplier="1 1 0.5 1 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Tissue volume br x2', 0.000063757685846, 0.020931764321828, 0.000095845591570, 0.401884069594298,
                 physiologyMultiplier="1 1 1 2 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Tissue volume br x0.5', 0.000064159478182, 0.020942830269960, 0.000023841842712, 0.100588868560786,
                 physiologyMultiplier="1 1 1 0.5 1 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Tissue volume alv x2', 0.000063683703509, 0.020933726571695, 0.000047771314560, 0.201085296873244,
                 physiologyMultiplier="1 1 1 1 2 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Tissue volume alv x0.5', 0.000064247463985, 0.020941858364839, 0.000047759016417, 0.201106474860665,
                 physiologyMultiplier="1 1 1 1 0.5 1 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Surface area br x2', 0.000064061396015, 0.020939149041018, 0.000047763114995, 0.201099418846000,
                 physiologyMultiplier="1 1 1 1 1 2 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Surface area br x0.5', 0.000064061396015, 0.020939149041018, 0.000047763114995, 0.201099418846000,
                 physiologyMultiplier="1 1 1 1 1 0.5 1 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Surface area alv x2', 0.000064061465104, 0.020939155418676, 0.000047763100214, 0.201099435470356,
                 physiologyMultiplier="1 1 1 1 1 1 2 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Surface area alv x0.5', 0.000064061256795, 0.020939136285631, 0.000047763144557, 0.201099385597105,
                 physiologyMultiplier="1 1 1 1 1 1 0.5 1 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Perfusion br x2', 0.000064086164157, 0.020941917285760, 0.000029578336958, 0.127803936169596,
                 physiologyMultiplier="1 1 1 1 1 1 1 2 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Perfusion br x0.5', 0.000063936357042, 0.020933596057896, 0.000084313842427, 0.347532487928486,
                 physiologyMultiplier="1 1 1 1 1 1 1 0.5 1")
        simulate(protDepoFPDA250_x1,  'FP', 'Perfusion alv x2', 0.000064065174363, 0.020939390820754, 0.000047762555549, 0.201100049079817,
                 physiologyMultiplier="1 1 1 1 1 1 1 1 2")
        simulate(protDepoFPDA250_x1,  'FP', 'Perfusion alv x0.5', 0.000064052167212, 0.020938665356367, 0.000047764234239, 0.201098158051882,
                 physiologyMultiplier="1 1 1 1 1 1 1 1 0.5")
        simulate(protDepoFPDA250_x2,  'FP', 'Particle size x2', 0.000014337075251, 0.008736680705496, 0.000020425178372, 0.065368750134521,
                 diameters="0.1,1.1,0.1; 1.2,24.2,0.2; 25,33,1", volMultiplier=2)
        simulate(protDepoFPDA250_x05,  'FP', 'Particle size x0.5', 0.000173018319322, 0.031918987519192, 0.000056919164003, 0.402057950230398,
                 diameters="0.1,1.1,0.1; 1.2,24.2,0.2; 25,33,1", volMultiplier=0.5)

        # Plots
        AUCratio = np.divide(AUC_brtis,AUC_sys)
        AUCincr_brtis = 100 * np.divide(AUC_brtis[1::2],AUC_brtis[0]);
        AUCdecr_brtis = 100 * np.divide(AUC_brtis[2::2],AUC_brtis[0]);
        Ctincr_brtis = 100 * np.divide(Ct_brtis[1::2],Ct_brtis[0]);
        Ctdecr_brtis = 100 * np.divide(Ct_brtis[2::2],Ct_brtis[0]);
        AUCincr_sys = 100 * np.divide(AUC_sys[1::2],AUC_sys[0]);
        AUCdecr_sys = 100 * np.divide(AUC_sys[2::2],AUC_sys[0]);
        Cmaxincr_sys = 100 * np.divide(Cmax_sys[1::2],Cmax_sys[0]);
        Cmaxdecr_sys = 100 * np.divide(Cmax_sys[2::2],Cmax_sys[0]);
        AUCratioincr = 100 * np.divide(AUCratio[1::2],AUCratio[0]);
        AUCratiodecr = 100 * np.divide(AUCratio[2::2],AUCratio[0]);

        plot_fact = 10 / 3;
        plot_range = 100 * np.asarray([1/plot_fact, plot_fact]);
        XTicks = [30, 50, 100, 150, 200, 300];
        ordering = [x-1 for x in [2, 6, 8, 1, 4, 5, 3, 9, 7, 15, 12, 16, 14, 18, 11, 10, 17, 13, 19]]
        pos = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21] # [1:9 11:19 21];

        par = [
           # Substance related
           'Solubility (br)',
           'Solubility (alv)',
           'Dissolution rate (br)',
           'Dissolution rate (alv)',
           'Permeability (br)',
           'Permeability (alv)',
           'Partition coefficient (br)',
           'Partition coefficient (alv)',
           'Systemic clearance',
           # physiological
           'Mucociliary clearance',
           'Fluid volume (br)',
           'Fluid volume (alv)',
           'Tissue volume (br)',
           'Tissue volume (alv)',
           'Surface area (br)',
           'Surface area (alv)',
           'Perfusion (br)',
           'Perfusion (alv)',
           # deposition-dependent
           'Particle size'
        ]

        import matplotlib.pyplot as plt
        import matplotlib.ticker as mtick
        plt.figure(figsize=(15, 9))
        plt.title('AUC (0-24h) in lung tissue')
        plt.xlabel('% of reference value')
        plt.grid()
        plt.yticks(pos, labels=[par[x] for x in ordering])
        baseline = 100
        plt.barh(pos, [x-baseline for x in (AUCincr_brtis[ordering]).tolist()], color='b')
        plt.barh(pos, [x-baseline for x in (AUCdecr_brtis[ordering]).tolist()], color='r')
        plt.legend(['2-fold parameter increase','2-fold parameter decrease'])
        plt.gca().xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: x + baseline))
        plt.savefig('AUC0_24_lungtissue.png')

        plt.figure(figsize=(15, 9))
        plt.title('Average lung tissue concentration after 24h')
        plt.xlabel('% of reference value')
        plt.grid()
        plt.yticks(pos, labels=[par[x] for x in ordering])
        baseline = 100
        plt.barh(pos, [x - baseline for x in (Ctincr_brtis[ordering]).tolist()], color='b')
        plt.barh(pos, [x - baseline for x in (Ctdecr_brtis[ordering]).tolist()], color='r')
        plt.legend(['2-fold parameter increase', '2-fold parameter decrease'])
        plt.gca().xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: x + baseline))
        plt.savefig('lungtissue24.png')

        plt.figure(figsize=(15, 9))
        plt.title('AUC (0-24h) ratio (lung selectivity)')
        plt.xlabel('% of reference value')
        plt.grid()
        plt.yticks(pos, labels=[par[x] for x in ordering])
        baseline = 100
        plt.barh(pos, [x - baseline for x in (AUCratioincr[ordering]).tolist()], color='b')
        plt.barh(pos, [x - baseline for x in (AUCratiodecr[ordering]).tolist()], color='r')
        plt.legend(['2-fold parameter increase', '2-fold parameter decrease'])
        plt.gca().xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: x + baseline))
        plt.savefig('AUC0_24ratio.png')

        plt.figure(figsize=(15, 9))
        plt.title('Systemic AUC (0-24h)')
        plt.xlabel('% of reference value')
        plt.grid()
        plt.yticks(pos, labels=[par[x] for x in ordering])
        baseline = 100
        plt.barh(pos, [x-baseline for x in (AUCincr_sys[ordering]).tolist()], color='b')
        plt.barh(pos, [x-baseline for x in (AUCdecr_sys[ordering]).tolist()], color='r')
        plt.legend(['2-fold parameter increase','2-fold parameter decrease'])
        plt.gca().xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: x + baseline))
        plt.savefig('AUC0_24_systemic.png')

        plt.figure(figsize=(15, 9))
        plt.title('Systemic Cmax')
        plt.xlabel('% of reference value')
        plt.grid()
        plt.yticks(pos, labels=[par[x] for x in ordering])
        baseline = 100
        plt.barh(pos, [x - baseline for x in (Cmaxincr_sys[ordering]).tolist()], color='b')
        plt.barh(pos, [x - baseline for x in (Cmaxdecr_sys[ordering]).tolist()], color='r')
        plt.legend(['2-fold parameter increase', '2-fold parameter decrease'])
        plt.gca().xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: x + baseline))
        plt.savefig('Cmax_systemic.png')

if __name__ == "__main__":
    unittest.main()
