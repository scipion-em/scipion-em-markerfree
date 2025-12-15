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

import logging
from os import stat
from os.path import exists
from typing import Tuple, List
import numpy as np
from pwem.objects import Transform
from pyworkflow.object import Pointer
from tomo.protocols import ProtTomoBase
from tomo.objects import TiltSeries, TiltImage

logger = logging.getLogger(__name__)

class ProtMarkerfreeBaseTsAlign(ProtTomoBase):
    pass
