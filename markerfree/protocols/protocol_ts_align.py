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
import subprocess
import traceback
import typing
from collections import Counter

from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
from pyworkflow.protocol import STEPS_PARALLEL, ProtStreamingBase
from pyworkflow.utils import Message, cyanStr, redStr

from markerfree import Plugin
from markerfree.constants import *

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfTiltSeries, TiltSeries

logger = logging.getLogger(__name__)

class ProtMarkerfreeAlignTiltSeries(ProtMarkerfreeBaseTsAlign, ProtStreamingBase):
    """Protocol to align tilt series using MarkerFree.
    """

    _label = 'align tilt series'
    _possibleOutputs = {OUTPUT_TILTSERIES_NAME: SetOfTiltSeries}
    stepsExecutionMode = STEPS_PARALLEL
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tsReadList = []

    @classmethod
    def worksInStreaming(cls):
        return True

    def _defineParams(self, form):

        form.addParallelSection(threads = 0, mpi = 1)
        form.addHidden(params.USE_GPU, params.BooleanParam,
                       default=True,
                       label="Use GPU for execution",
                       help="This protocol uses GPU."
                            "Select the one you want to use.")
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       label="Choose GPU ID",
                       help="You may have several GPUs. Set it to zero"
                            " if you do not know what we are talking about."
                            " First GPU index is 0, second 1 and so on."
                            " You can use many GPUs")
        form.addParam('inTiltSeries', params.PointerParam, label="Input TS set",
                      pointerClass='SetOfTiltSeries', important=True)

        line = form.addLine('Geometry',
                            help="These options govern TM geometry.")
        line.addParam('geom_offset', params.IntParam, label="Offset")
        line.addParam('geom_tiltAxisAngle', params.IntParam, label="Offset")
        line.addParam('geom_zAxisOffset', params.IntParam, label="Offset")
        line.addParam('geom_thickness', params.IntParam, label="Offset")
        line.addParam('geom_reconThickness', params.IntParam, label="Offset")
        line.addParam('geom_downsample', params.IntParam, label="Offset")

        form.addParam('nProjs', params.IntParam, default=10, 
                      label="Projections",
                      help="Number of projections to use in the projection"
                      "matching phase.")
        
    def stepsGeneratorStep(self) -> None:
        """
        This is the same step implemented in any streaming protocol.
        It will check its input and insert new FunctionStep's when
        the conditions are met.
        """
        closeSetStepDeps = []
        inTsSet = self.getInputTsSet()
        self.readingOutput(getattr(self, OUTPUT_TILTSERIES_NAME, None))

        while True:
            with self._lock:
                listTSInput = inTsSet.getTSIds()
            if not inTsSet.isStreamOpen() and Counter(self.tsReadList) == Counter(listTSInput):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self.closeOutputSetsStep,
                                         OUTPUT_TILTSERIES_NAME,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break
            closeSetStepDeps = []
            for ts in inTsSet.iterItems():
                tsId = ts.getTsId()
                if tsId not in self.tsReadList and ts.getSize() > 0:  # Avoid processing empty TS (before the Tis are added)
                        cInputId = self._insertFunctionStep(self.convertInputStep, tsId,
                                                            prerequisites=[],
                                                            needsGPU=False)
                        tsAlignId = self._insertFunctionStep(self.runMarkerfreeStep, tsId,
                                                             prerequisites=cInputId,
                                                             needsGPU=True)
                        cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                                          prerequisites=tsAlignId,
                                                          needsGPU=False)
                        closeSetStepDeps.append(cOutId)
                        logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                        self.tsReadList.append(tsId)

            self.refreshStreaming(inTsSet)

    def convertInputStep(self, tsId: str):
        with self._lock:
            ts = self.getCurrentTs(tsId)


    def runMarkerfreeStep(self, tsId: str):
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId}: aligning...'))
                with self._lock:
                    ts = self.getCurrentTs(tsId)
                self.genTsPaths(tsId)
                args = self._getCmd(ts)
                Plugin.runMarkerfree(self, args)
            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> {TSALIGN_PROGRAM} execution failed with the exception -> {e}'))
                logger.error(traceback.format_exc())
