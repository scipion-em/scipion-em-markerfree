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

# MF
VERSION_0_1 = '0.1'
VERSIONS = [VERSION_0_1]
DEFAULT_VERSION = VERSION_0_1

MARKERFREE = 'markerfree'
MARKERFREE_HOME = 'MARKERFREE_HOME'
MARKERFREE_ENV_NAME = '%s-%s' % (MARKERFREE, DEFAULT_VERSION)
MARKERFREE_DEFAULT_ACTIVATION_CMD = 'conda activate %s' % MARKERFREE_ENV_NAME
MARKERFREE_CUDA_LIB = 'MARKERFREE_CUDA_LIB'
MARKERFREE_ENV_ACTIVATION = 'MARKERFREE_ENV_ACTIVATION'
MARKERFREE_CMD = 'Markerfree'

# Programs
TSALIGN_PROGRAM = 'Markerfree'

OUTPUT_TILTSERIES_NAME = "TiltSeries"
OUTPUT_ALI_TILTSERIES_NAME = "AlignedTiltSeries"
OUTPUT_TS_FAILED_NAME = "FailedTiltSeries"

MRCS_EXT = 'mrcs'
MRC_EXT = 'mrc'
XF_EXT = '.xf'
TLT_EXT = 'tlt'
TXT_EXT = 'txt'

NONE_PROCESSED_MSG = 'Unable to process any of the introduced'
NO_TS_PROCESSED_MSG = f'{NONE_PROCESSED_MSG} tilt-series.'
