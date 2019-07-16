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

import os
import pyworkflow.em
from pyworkflow.utils import Environ
from .constants import *


_references = ['']
_logo = ''


class Plugin(pyworkflow.em.Plugin):
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

    # @classmethod
    # def isVersionActive(cls):
    #     return cls.getActiveVersion().startswith(V1_0_0)
    #
    @classmethod
    def defineBinaries(cls, env):
        env.addPipModule('scipy-015', pipCmd=env._pipCmd % ('scipy', '0.15'),
                         target='scipy')
        scons = tryAddPipModule(env, 'openpyxl', '2.6.2')

def tryAddPipModule(env, moduleName, *args, **kwargs):
    """ To try to add certain pipModule.
        If it fails due to it is already add by other plugin or Scipion,
          just returns its name to use it as a dependency.
        Raise the exception if unknown error is gotten.
    """
    try:
        return env.addPipModule(moduleName, *args, **kwargs)._name
    except Exception as e:
        if str(e) == "Duplicated target '%s'" % moduleName:
            return moduleName
        else:
            raise Exception(e)


pyworkflow.em.Domain.registerPlugin(__name__)
