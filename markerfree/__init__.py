# **************************************************************************
# *
# * Authors:     Mikel Iceta (miceta@cnb.csic.es)
# *
# * National Center of Biotechnology (CNB-CSIC)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import os
from shutil import which

import pwem
import pyworkflow.utils as pwutils

from markerfree.constants import *

__version__ = "3.0.0"
# _references = []
# _logo

class Plugin(pwem.Plugin):
    _homeVar = MARKERFREE_HOME
    _url = "https://github.com/scipion-em/scipion-em-markerfree"
    _validationMsg = None

    @classmethod
    def _devineVariables(cls):
        cls._defineEmVar(MARKERFREE_HOME, cls._getMarkerfreeFolder(DEFAULT_VERSION))
        cls._defineVar(MARKERFREE_ENV_ACTIVATION, MARKERFREE_DEFAULT_ACTIVATION_CMD)
        cls._defineVar(MARKERFREE_CUDA_LIB, pwem.Config.CUDA_LIB)

    @classmethod
    def getMarkerfreeEnvActivation(cls):
        return cls.getVar(MARKERFREE_ENV_ACTIVATION)
    
    @classmethod
    def _getEMFolder(cls, version, *paths):
        return os.path.join("markerfree-%s" % version, *paths)
    
    @classmethod
    def _getMarkerfreeFolder(cls, version, *paths):
        return cls._getEMFolder(version, "Markerfree", *paths)
    
    @classmethod
    def _getProgram(cls, program):
        """ Returns the same program  if config missing
        or the path to the program based on the config file."""
        # Compose path based on config
        progFromConfig = cls.getHome("bin", program)

        # Check if IMOD from config exists
        if os.path.exists(progFromConfig):
            return progFromConfig
        else:
            return program
        
    @classmethod
    def validateInstallation(cls):
        """ Check if imod is in the path """

        if not cls._validationMsg:
            mkfr = cls._getProgram(MARKERFREE_CMD)

            cls._validationMsg = [
                "MKFR's %s command not found in path, please "
                "install it." % mkfr] if not which(
                MARKERFREE_CMD) and not os.path.exists(mkfr) else []

        return cls._validationMsg
    
    @classmethod
    def installMarkerfree(cls, env):
        MARKERFREE_INSTALLED = '%s_%s_installed' % (MARKERFREE, DEFAULT_VERSION)
        installationCmd = cls.getCondaActivationCmd()
        
        # Create an environment
        installationCmd += ' conda create -y -n %s -c conda-forge python=3.8 && ' % MARKERFREE_ENV_NAME

        # Activate env
        installationCmd += ' conda activate %s && ' % MARKERFREE_ENV_NAME

        # Install MKFR
        # TODO: COMPILE

        MKFR_commands = [(installationCmd, MARKERFREE_INSTALLED)]
        envPath = os.environ.get('PATH', "") # Keep path, Conda is likely there
        installEnvVars = {'PATH': envPath} if envPath else None

        env.addPackage(MARKERFREE,
                       version=DEFAULT_VERSION,
                       tar='void.tgz',
                       commands=MKFR_commands,
                       neededProgs=cls.getDependenciesMKFR(),
                       vars=installEnvVars,
                       default=True)
        
    @classmethod
    def getDependenciesMKFR(cls):
        # Try to activate conda
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')
        return neededProgs
    
    @classmethod
    def runMarkerfree(cls, protocol, args, cwd=None, numberOfMpi=1):
        """ Run Markerfree command from a given protocol. """
        cmd = cls.getCondaActivationCmd() + " "
        cmd += cls.getMarkerfreeEnvActivation() + " "
        cmd += f"&& export PATH={cls.getHome('build/bin')}:PATH "
        cmd += f"&& {TSALIGN_PROGRAM}"
        protocol.runJob(cmd, args, env=cls.getEnviron(), cwd=cwd, numberOfMpi=1)
