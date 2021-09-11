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

class TestInhalation7Workflow(TestWorkflow):
    # Hartung2020_MATLAB/scripts/simulation_sensitivityAnalysis_Bud.m

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
        print("Substance (budesonide) ...")
        protSubstBud=self.newProtocol(ProtPKPDInhSubstanceProperties,
                                     objLabel='pkpd - budesonide properties 2',
                                     name='budesonide',
                                     rho=1.3, MW=430.53,
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
        protDepoBudA800_x1 = launchDepo('Bud Turbohaler Asthma 800ug', 'Bud_Turbohaler_asthmatic_800ug.txt', protSubstBud)
        protDepoBudA800_x2 = launchDepo('Bud Turbohaler Asthma 800ug x2', 'Bud_Turbohaler_asthmatic_800ug_2FoldLarger.txt', protSubstBud)
        protDepoBudA800_x05 = launchDepo('Bud Turbohaler Asthma 800ug x0.5', 'Bud_Turbohaler_asthmatic_800ug_2FoldSmaller.txt', protSubstBud)

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
            if subst=="Bud":
                simulationTime = 12*60
                deltaT = simulationTime/20000
                protPK = protPKBud

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

        simulate(protDepoBudA800_x1,  'Bud', 'Reference', 0.001201526822791, 0.354301881503904, 0.005404220428749, 11.153499130558401)
        simulate(protDepoBudA800_x1,  'Bud', 'Solubility br x2', 0.001215249991232, 0.364751343404497, 0.005557880051374, 12.623390558321629,
                 substanceMultiplier="2 1 1 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Solubility br x0.5', 0.001186305288001, 0.342424052023733, 0.004819626351434, 9.418317380724629,
                 substanceMultiplier="0.5 1 1 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Solubility alv x2', 0.001201536074, 0.354301966246, 0.005404216073, 11.153500025985,
                 substanceMultiplier="1 2 1 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Solubility alv x0.5', 0.001201508343, 0.354301712013, 0.005404229141, 11.153497339648,
                 substanceMultiplier="1 0.5 1 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Dissolution rate br x2', 0.001227952943, 0.368997900815, 0.005801659439, 13.271698208948,
                 substanceMultiplier="1 1 2 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Dissolution rate br x0.5', 0.001179676390, 0.340687336041, 0.004762102000, 9.150525483341,
                 substanceMultiplier="1 1 0.5 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Dissolution rate alv x2', 0.001279708140, 0.354667170527, 0.005385444582, 11.157358986897,
                 substanceMultiplier="1 1 1 2 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Dissolution rate alv x0.5', 0.001084542590, 0.353383210864, 0.005451400917, 11.143791697266,
                 substanceMultiplier="1 1 1 0.5 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Permeability br x2', 0.001222087033, 0.369563298392, 0.005946524395, 13.270632796688,
                 substanceMultiplier="1 1 1 1 2 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Permeability br x0.5', 0.001183387572, 0.340832933134, 0.004581676018, 9.214131272824,
                 substanceMultiplier="1 1 1 1 0.5 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Permeability alv x2', 0.001201546308, 0.354303355338, 0.005404144677, 11.153514702969,
                 substanceMultiplier="1 1 1 1 1 2 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Permeability alv x0.5', 0.001201486505, 0.354298932965, 0.005404371977, 11.153467976528,
                 substanceMultiplier="1 1 1 1 1 0.5 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Partition coefficient br x2', 0.001182898425, 0.353413245682, 0.011050937432, 22.164490945241,
                 substanceMultiplier="1 1 1 1 1 1 2 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Partition coefficient br x0.5', 0.001211654566, 0.354732091309, 0.002673296897, 5.593978829352,
                 substanceMultiplier="1 1 1 1 1 1 0.5 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Partition coefficient alv x2', 0.001188247165, 0.353630921033, 0.005423243646, 11.146579433854,
                 substanceMultiplier="1 1 1 1 1 1 1 2")
        simulate(protDepoBudA800_x1,  'Bud', 'Partition coefficient alv x0.5', 0.001208200016, 0.354632504525, 0.005394732660, 11.156910141397,
                 substanceMultiplier="1 1 1 1 1 1 1 0.5")
        simulate(protDepoBudA800_x1,  'Bud', 'Systemic clearance x2', 0.001023319331, 0.185236386141, 0.004583124259, 9.471920613496,
                 pkMultiplier="2 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Systemic clearance x0.5', 0.001344706991, 0.594367656929, 0.008264150735, 13.522551410075,
                 pkMultiplier="0.5 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Bioavailability x2', 0.001598714792, 0.457373704413, 0.005586996435, 12.182194162445,
                 pkMultiplier="1 1 1 1 1 2")
        simulate(protDepoBudA800_x1,  'Bud', 'Bioavailability x0.5', 0.001006166885, 0.302765972014, 0.005312832435, 10.639151634150,
                 pkMultiplier="1 1 1 1 1 0.5")
        simulate(protDepoBudA800_x1,  'Bud', 'Mucociliary clearance x2', 0.001213312664, 0.346463121087, 0.004038271823, 9.624049361203,
                 physiologyMultiplier="2 1 1 1 1 1 1 1 1", ciliarySpeedType=1)
        simulate(protDepoBudA800_x1,  'Bud', 'Mucociliary clearance x0.5', 0.001180163892, 0.375489909837, 0.005415777203, 14.267137616955,
                 physiologyMultiplier="0.5 1 1 1 1 1 1 1 1", ciliarySpeedType=1)
        simulate(protDepoBudA800_x1,  'Bud', 'Fluid volume br x2', 0.001195881254, 0.350488722832, 0.004950826813, 10.638338474227,
                 physiologyMultiplier="1 2 1 1 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Fluid volume br x0.5', 0.001205027051, 0.356370214148, 0.005648126925, 11.429902033095,
                 physiologyMultiplier="1 0.5 1 1 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Fluid volume alv x2', 0.001201479412, 0.354297856668, 0.005404398993, 11.153456915002,
                 physiologyMultiplier="1 1 2 1 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Fluid volume alv x0.5', 0.001201549895, 0.354303892660, 0.005404131202, 11.153520225101,
                 physiologyMultiplier="1 1 0.5 1 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Tissue volume br x2', 0.001182898425, 0.353413245682, 0.011050937432, 22.164490945241,
                 physiologyMultiplier="1 1 1 2 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Tissue volume br x0.5', 0.001211654566, 0.354732091309, 0.002673296897, 5.593978829352,
                 physiologyMultiplier="1 1 1 0.5 1 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Tissue volume alv x2', 0.001188247165, 0.353630921033, 0.005423243646, 11.146579433854,
                 physiologyMultiplier="1 1 1 1 2 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Tissue volume alv x0.5', 0.001208200016, 0.354632504525, 0.005394732660, 11.156910141397,
                 physiologyMultiplier="1 1 1 1 0.5 1 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Surface area br x2', 0.001201526823, 0.354301881504, 0.005404220429, 11.153499130558,
                 physiologyMultiplier="1 1 1 1 1 2 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Surface area br x0.5', 0.001201526823, 0.354301881504, 0.005404220429, 11.153499130558,
                 physiologyMultiplier="1 1 1 1 1 0.5 1 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Surface area alv x2', 0.001201546308, 0.354303355338, 0.005404144677, 11.153514702969,
                 physiologyMultiplier="1 1 1 1 1 1 2 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Surface area alv x0.5', 0.001201486505, 0.354298932965, 0.005404371977, 11.153467976528,
                 physiologyMultiplier="1 1 1 1 1 1 0.5 1 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Perfusion br x2', 0.001209175911, 0.354613419397, 0.003204896790, 7.363561656543,
                 physiologyMultiplier="1 1 1 1 1 1 1 2 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Perfusion br x0.5', 0.001187232599, 0.353664155102, 0.009904871423, 18.657914258309,
                 physiologyMultiplier="1 1 1 1 1 1 1 0.5 1")
        simulate(protDepoBudA800_x1,  'Bud', 'Perfusion alv x2', 0.001201677614, 0.354314572750, 0.005403585248, 11.153633037620,
                 physiologyMultiplier="1 1 1 1 1 1 1 1 2")
        simulate(protDepoBudA800_x1,  'Bud', 'Perfusion alv x0.5', 0.001201128766, 0.354276360663, 0.005405497649, 11.153229856200,
                 physiologyMultiplier="1 1 1 1 1 1 1 1 0.5")
        simulate(protDepoBudA800_x2,  'Bud', 'Particle size x2', 0.000877633922, 0.286477232680, 0.003452813688, 6.923518661367,
                 diameters="0.1,1.1,0.1; 1.2,24.2,0.2; 25,33,1", volMultiplier=2)
        simulate(protDepoBudA800_x05,  'Bud', 'Particle size x0.5', 0.001353341968, 0.389232974573, 0.005399471025, 15.308840049263,
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
        ordering = [x-1 for x in [10, 2, 6, 8, 1, 4, 5, 3, 9, 7, 16, 13, 17, 15, 19, 12, 11, 18, 14, 20]]
        pos = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22] # [1:10 12:20 22]

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
           'Oral bioavailability',
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
        plt.title('AUC (0-12h) in lung tissue')
        plt.xlabel('% of reference value')
        plt.grid()
        plt.yticks(pos, labels=[par[x] for x in ordering])
        baseline = 100
        plt.barh(pos, [x-baseline for x in (AUCincr_brtis[ordering]).tolist()], color='b')
        plt.barh(pos, [x-baseline for x in (AUCdecr_brtis[ordering]).tolist()], color='r')
        plt.legend(['2-fold parameter increase','2-fold parameter decrease'])
        plt.gca().xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: x + baseline))
        plt.savefig('AUC0_12_lungtissue.png')

        plt.figure(figsize=(15, 9))
        plt.title('Average lung tissue concentration after 12')
        plt.xlabel('% of reference value')
        plt.grid()
        plt.yticks(pos, labels=[par[x] for x in ordering])
        baseline = 100
        plt.barh(pos, [x - baseline for x in (Ctincr_brtis[ordering]).tolist()], color='b')
        plt.barh(pos, [x - baseline for x in (Ctdecr_brtis[ordering]).tolist()], color='r')
        plt.legend(['2-fold parameter increase', '2-fold parameter decrease'])
        plt.gca().xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: x + baseline))
        plt.savefig('lungtissue12.png')

        plt.figure(figsize=(15, 9))
        plt.title('AUC (0-12h) ratio (lung selectivity)')
        plt.xlabel('% of reference value')
        plt.grid()
        plt.yticks(pos, labels=[par[x] for x in ordering])
        baseline = 100
        plt.barh(pos, [x - baseline for x in (AUCratioincr[ordering]).tolist()], color='b')
        plt.barh(pos, [x - baseline for x in (AUCratiodecr[ordering]).tolist()], color='r')
        plt.legend(['2-fold parameter increase', '2-fold parameter decrease'])
        plt.gca().xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: x + baseline))
        plt.savefig('AUC0_12ratio.png')

        plt.figure(figsize=(15, 9))
        plt.title('Systemic AUC (0-12h)')
        plt.xlabel('% of reference value')
        plt.grid()
        plt.yticks(pos, labels=[par[x] for x in ordering])
        baseline = 100
        plt.barh(pos, [x-baseline for x in (AUCincr_sys[ordering]).tolist()], color='b')
        plt.barh(pos, [x-baseline for x in (AUCdecr_sys[ordering]).tolist()], color='r')
        plt.legend(['2-fold parameter increase','2-fold parameter decrease'])
        plt.gca().xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: x + baseline))
        plt.savefig('AUC0_12_systemic.png')

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
