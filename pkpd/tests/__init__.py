# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pkpd.objects import PKPDDataSet

try:
    from unittest.runner import _WritelnDecorator # Python 2.7+
except ImportError:
    from unittest import _WritelnDecorator # Python <2.6

PKPDDataSet(name='Gabrielsson_PK01', folder='Gabrielsson_PK01',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK02', folder='Gabrielsson_PK02',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK03', folder='Gabrielsson_PK03',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK04', folder='Gabrielsson_PK04',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK05', folder='Gabrielsson_PK05',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK06', folder='Gabrielsson_PK06',
        files={
               'experiment_ev': 'experiment_ev.pkpd',
               'experiment_iv': 'experiment_iv.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK07', folder='Gabrielsson_PK07',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK08', folder='Gabrielsson_PK08',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK09', folder='Gabrielsson_PK09',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK10', folder='Gabrielsson_PK10',
        files={
               'experimentIV': 'experimentIV.pkpd',
               'experimentPO': 'experimentPO.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK11', folder='Gabrielsson_PK11',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK12', folder='Gabrielsson_PK12',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK13', folder='Gabrielsson_PK13',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK14', folder='Gabrielsson_PK14',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK15', folder='Gabrielsson_PK15',
        files={
               'experiment': 'experiment.pkpd',
               'experimentPD': 'experimentPD.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK16', folder='Gabrielsson_PK16',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK17', folder='Gabrielsson_PK17',
        files={
               'experiment': 'experiment.pkpd',
               'experiment2': 'experiment2.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK18', folder='Gabrielsson_PK18',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK19', folder='Gabrielsson_PK19',
        files={
               'experiment': 'experiment.pkpd',
               'experimentIndividual': 'experimentIndividual.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK20', folder='Gabrielsson_PK20',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK22', folder='Gabrielsson_PK22',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK23', folder='Gabrielsson_PK23',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK25', folder='Gabrielsson_PK25',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK39', folder='Gabrielsson_PK39',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PK43', folder='Gabrielsson_PK43',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PD03', folder='Gabrielsson_PD03',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Gabrielsson_PD11', folder='Gabrielsson_PD11',
        files={
               'experiment': 'experiment.pkpd'
               })

PKPDDataSet(name='Dissolution', folder='Dissolution',
        files={
               'experiment': 'dissolution.pkpd',
               'experiment12Test': 'dissolution12Test.pkpd',
               'experiment12Ref': 'dissolution12Ref.pkpd',
               'invivo12': 'invivo12.pkpd',
               'excel': 'experiment.xlsx',
               'excelShah1998': 'experimentShah1998.xlsx',
               'excelShah1998avg': 'experimentShah1998Avg.xlsx',
               'exceldissolutionMean': 'dissolutionMeanProfile.xlsx',
               'excelinvivo': 'invivo.xlsx',
               'excelinvivoLong': 'invivoLong.xlsx',
               'excelinvivoMean': 'invivoMean.xlsx',
               'webplot': 'dataFromWebPlotDigitizer.csv',
               })

PKPDDataSet(name='model',  folder='model',
        files={
               'classesSelection': 'gold/classes_selection.sqlite',
               'modelGoldSqlite': 'gold/model_gold.sqlite',
               'modelGoldXml': 'gold/model_gold.xml',
        })
