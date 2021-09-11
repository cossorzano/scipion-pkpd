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

from pkpd.objects import PKPhysiologyLungParameters

# Tested in test_workflow_inhalation1

class ProtPKPDInhLungPhysiology(ProtPKPD):
    """ Produce a description of the lung physiological parameters\n
        See Hartung2020.
        Protocol created by http://www.kinestatpharma.com\n """
    _label = 'lung parameters'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        # Default parameters from
        # Hartung2020_MATLAB/physiol_subst/input_physiology.m

        form.addSection('Input')
        groupS = form.addGroup('Systemic parameters')
        groupS.addParam('Qco', params.FloatParam, label="Cardiac output [ml/min]", default=312/60*1000,
                      help='6000 from  Boger 2016, Table S3. 312/60*1000 from Brown 1997')
        groupS.addParam('OWlun', params.FloatParam, label="Lung tissue weight, without blood [g]", default=532,
                      help='532 from Brown 1997')

        groupA = form.addGroup('Alveolar parameters')
        groupA.addParam('falvCO', params.FloatParam, label="Alveolar fraction of cardiac output", default=1)
        groupA.addParam('Velf_alv', params.FloatParam, label="Alveolar ELF volume [mL]", default=36,
                      help='Fronius 2012')
        groupA.addParam('Surf_alv', params.FloatParam, label="Alveolar surface [cm2]", default=13000*100,
                      help='13000*100 from Weibel 2009')

        groupB = form.addGroup('Bronchial parameters')
        groupB.addParam('fbrCO', params.FloatParam, label="Bronchial fraction of cardiac output", default=0.025,
                      help='Gaohua 2015')
        groupB.addParam('fbrVlun', params.FloatParam, label="Bronchial fraction of lung tissue volume", default=0.27,
                      help='Boger 2016 (Supplement)')
        groupB.addParam('helf_trach', params.FloatParam, label="Bronchial ELF heights for interpolation, trachea [cm]",
                      default=10*1e-4, help='Hartung 2020, Fig.4 in Suppl. 1 Appendix')
        groupB.addParam('helf_termbr', params.FloatParam,
                      label="Bronchial ELF heights for interpolation, terminal bronchioles [cm]",
                      default=1.8*1e-4, help='Hartung 2020, Fig.4 in Suppl. 1 Appendix')

        form.addSection('Lung structure')
        groupT= form.addGroup('Trachea')
        groupT.addParam('tracheaLength', params.StringParam, label='Length [cm]', default='10.0',
                        help='Yeh/Schum 1980, Bull. Math. Biol.')
        groupT.addParam('tracheaDiameter', params.StringParam, label='Diameter [cm]', default='2.01',
                        help='Yeh/Schum 1980, Bull. Math. Biol.')

        groupB1 = form.addGroup('Main bronchi')
        groupB1.addParam('bronchi1Length', params.StringParam, label='Length [cm]', default='4.36',
                         help='Yeh/Schum 1980, Bull. Math. Biol.')
        groupB1.addParam('bronchi1Diameter', params.StringParam, label='Diameter [cm]', default='1.56',
                         help='Yeh/Schum 1980, Bull. Math. Biol.')

        groupB2 = form.addGroup('Bronchi')
        groupB2.addParam('bronchi2Length', params.StringParam, label='Length [cm]', default='1.78 0.965 0.995 1.01 0.890 0.962 0.867',
                         help='Define as many as required separated by spaces. Yeh/Schum 1980, Bull. Math. Biol.')
        groupB2.addParam('bronchi2Diameter', params.StringParam, label='Diameter [cm]', default='1.13 0.827 0.651 0.574 0.435 0.373 0.322',
                         help='Define as many as required separated by spaces. Yeh/Schum 1980, Bull. Math. Biol.')

        groupB3 = form.addGroup('Bronchioles')
        groupB3.addParam('bronchi3Length', params.StringParam, label='Length [cm]',
                         default='0.667 0.556 0.446 0.359 0.275 0.212',
                         help='Define as many as required separated by spaces. Yeh/Schum 1980, Bull. Math. Biol.')
        groupB3.addParam('bronchi3Diameter', params.StringParam, label='Diameter [cm]',
                         default='0.257 0.198 0.156 0.118 0.092 0.073',
                         help='Define as many as required separated by spaces. Yeh/Schum 1980, Bull. Math. Biol.')

        groupB4 = form.addGroup('Terminal bronchioles')
        groupB4.addParam('bronchi4Length', params.StringParam, label='Length [cm]',
                         default='0.168',
                         help='Yeh/Schum 1980, Bull. Math. Biol.')
        groupB4.addParam('bronchi4Diameter', params.StringParam, label='Diameter [cm]',
                         default='0.060',
                         help='Yeh/Schum 1980, Bull. Math. Biol.')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('runParams')
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def runParams(self):
        self.lungParams = PKPhysiologyLungParameters()
        self.lungParams.Qco=self.Qco.get()
        self.lungParams.OWlun=self.OWlun.get()
        self.lungParams.falvCO=self.falvCO.get()
        self.lungParams.Velf_alv=self.Velf_alv.get()
        self.lungParams.Surf_alv=self.Surf_alv.get()
        self.lungParams.fbrCO=self.fbrCO.get()
        self.lungParams.fbrVlun=self.fbrVlun.get()
        self.lungParams.helf_trach=self.helf_trach.get()
        self.lungParams.helf_termbr=self.helf_termbr.get()

        self.lungParams.tracheaLength=self.tracheaLength.get()
        self.lungParams.tracheaDiameter=self.tracheaDiameter.get()
        self.lungParams.bronchi1Length=self.bronchi1Length.get()
        self.lungParams.bronchi1Diameter=self.bronchi1Diameter.get()
        self.lungParams.bronchi2Length=','.join(self.bronchi2Length.get().split())
        self.lungParams.bronchi2Diameter=','.join(self.bronchi2Diameter.get().split())
        self.lungParams.bronchi3Length=','.join(self.bronchi3Length.get().split())
        self.lungParams.bronchi3Diameter=','.join(self.bronchi3Diameter.get().split())
        self.lungParams.bronchi4Length=self.bronchi4Length.get()
        self.lungParams.bronchi4Diameter=self.bronchi4Diameter.get()

        self.lungParams.write(self._getPath("lung_parameters.pkpd"))

    def createOutputStep(self):
        self._defineOutputs(outputLungParameters=self.lungParams)

    def _citations(self):
        return ['Hartung2020']