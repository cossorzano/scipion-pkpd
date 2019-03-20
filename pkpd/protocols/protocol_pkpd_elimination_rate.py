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

from .protocol_pkpd_exponential_fit import ProtPKPDExponentialFit

# TESTED in test_workflow_gabrielsson_pk01.py
# TESTED in test_workflow_gabrielsson_pk02.py
# TESTED in test_workflow_gabrielsson_pk06.py
# TESTED in test_workflow_gabrielsson_pk07.py

class ProtPKPDEliminationRate(ProtPKPDExponentialFit):
    """ Fit a single exponential to the input data.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'elimination rate'

    def _defineParams(self, form):
        ProtPKPDExponentialFit._defineParams(self, form, fullForm=False)
