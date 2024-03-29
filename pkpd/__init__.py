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
"""
PKPD functions
"""

# Pending:
# Batch effects, Reese2013

import distutils.spawn
import os
import subprocess
import uuid
import pwem as em
from pyworkflow.utils import Environ
from pyworkflow.utils.path import cleanPath
from .constants import *


_references = ['']
_logo = ''
__version__ = '1.0.2'


class Plugin(em.Plugin):
    _homeVar = PKDP_HOME

    # @classmethod
    # def _defineVariables(cls):
    #     cls._defineEmVar(PKDP_HOME, 'pkdp-1.0.0')

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch pkdp """
        environ = Environ(os.environ)

        # environ.update({
        #     'PATH': Plugin.getHome(),
        #     'LD_LIBRARY_PATH': str.join(cls.getHome(), 'pkdplib')
        #                        + ":" + cls.getHome(),
        # }, position=Environ.BEGIN)

        return environ

    @classmethod
    def getRscript(cls):
        return distutils.spawn.find_executable("Rscript")

    @classmethod
    def runRscript(cls, scriptString):
        Rscript = cls.getRscript()
        if Rscript is None:
            raise Exception("Cannot find Rscript in the path")

        fnRoot = str(uuid.uuid4()).replace('-','_')
        fnScript = '/tmp/scipionRScript_%s.r'%fnRoot
        fh = open(fnScript,'w')
        fh.write(scriptString)
        fh.close()

        print("Running: ",Rscript,"--vanilla",fnScript)
        process = subprocess.Popen([Rscript,"--vanilla",fnScript], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = ''.join([x.decode() for x in process.stdout.readlines()])

        cleanPath(fnScript)
        return output

    @classmethod
    def checkRPackage(cls, packageName):
        response = cls.runRscript("if ('%s' %%in%% installed.packages()) {cat(1);} else {cat(0)}"%packageName)
        return response=='1'
