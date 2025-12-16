# **************************************************************************
# *
# * Authors:     Mikel Iceta (miceta@cnb.csic.es)
# * Authors:     JL Vilas (jl.vilas@cnb.csic.es)
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
import traceback
import time
import os
from collections import Counter
from enum import Enum
from typing import List, Tuple, Union
import numpy as np

from pwem.emlib import DT_FLOAT
from pwem.protocols import EMProtocol

from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
from pyworkflow.object import Set, Pointer
from pyworkflow.protocol import STEPS_PARALLEL, LEVEL_ADVANCED, ProtStreamingBase
from pyworkflow.utils import makePath, cyanStr, redStr

from markerfree import Plugin
from markerfree.constants import *
from markerfree.convert import readXfFile

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage, Pointer

logger = logging.getLogger(__name__)

# Form variables
IN_TS_SET = 'inTsSet'

# Auxiliar variables
EVEN_SUFFIX = '_even'
ODD_SUFFIX = '_odd'
IDENTITY_MATRIX = np.eye(3)  # Store in memory instead of multiple creation

class markerfreeOutputs(Enum):
    tiltSeries = SetOfTiltSeries

class ProtMarkerfreeAlignTiltSeries(EMProtocol, ProtTomoBase, ProtStreamingBase):
    """Protocol to align tilt series using MarkerFree.
    """

    _label = 'align tilt series'
    _possibleOutputs = markerfreeOutputs
    stepsExecutionMode = STEPS_PARALLEL
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.itemTsIdReadList = []
        self.failedItems = []

    @classmethod
    def worksInStreaming(cls):
        return True

    def _defineParams(self, form):
        
        form.addSection(label="Input")
        # form.addHidden(params.USE_GPU, params.BooleanParam,
        #                default=True,
        #                label="Use GPU for execution",
        #                help="This protocol uses GPU."
        #                     "Select the one you want to use.")
        # form.addHidden(params.GPU_LIST, params.StringParam, default='0',
        #                label="Choose GPU ID",
        #                help="You may have several GPUs. Set it to zero"
        #                     " if you do not know what we are talking about."
        #                     " First GPU index is 0, second 1 and so on."
        #                     " You can use many GPUs")
        form.addParam(IN_TS_SET, params.PointerParam, label="Tilt Series",
                      pointerClass='SetOfTiltSeries', important=True)

        form.addParam('geomOffset', params.IntParam, expertLevel=LEVEL_ADVANCED, default=0, label="Offset")
        form.addParam('geomZAxisOffset', params.IntParam, expertLevel=LEVEL_ADVANCED, default=0, label="Z axis offset")
        form.addParam('geomThickness', params.IntParam, default=0, label="Thickness")
        form.addParam('geomReconThickness', params.IntParam, default=0, label="Reconstruction thickness")
        form.addParam('geomDownsample', params.IntParam, default=0, label="Downsample factor")

        form.addParam('nProjs', params.IntParam, expertLevel=LEVEL_ADVANCED, default=10, 
                      label="Projections",
                      help="Number of projections to use in the projection"
                      "matching phase.")
        '''
        form.addParam('doReconstruction', params.BooleanParam,
                      label='Reconstruct tomogram?',
                      default=True)
        form.addParam('doEvenOdd', params.BooleanParam,
                      condition='doReconstruction=True',
                      label='Reconstruct odd/even tilt-series?',
                      default=False)
        '''
        
    def stepsGeneratorStep(self) -> None:
        closeSetStepDeps = []
        inTsSet = self._getInTsSet()

        self.readingOutput()

        while True:
            with self._lock:
                listInTsIds = inTsSet.getTSIds()
            # In the if statement below, Counter is used because in the tsId comparison the order doesn’t matter
            # but duplicates do. With a direct comparison, the closing step may not be inserted because of the order:
            # ['ts_a', 'ts_b'] != ['ts_b', 'ts_a'], but they are the same with Counter.
            if not inTsSet.isStreamOpen() and Counter(self.itemTsIdReadList) == Counter(listInTsIds):
                logger.info(cyanStr('Input set closed.\n'))
                self._insertFunctionStep(self._closeOutputSet,
                                         prerequisites=closeSetStepDeps,
                                         needsGPU=False)
                break

            for ts in inTsSet.iterItems():
                tsId = ts.getTsId()
                
                if tsId not in self.itemTsIdReadList and ts.getSize() > 0:  # Avoid processing empty TS (wait for TS imgs to be added)
                    tsAlignId = self._insertFunctionStep(self.runMarkerfreeStep, tsId,
                                                            prerequisites=[],
                                                            needsGPU=True)
                    cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                                        prerequisites=tsAlignId,
                                                        needsGPU=False)
                    closeSetStepDeps.append(cOutId)
                    logger.info(cyanStr(f"Steps created for tsId = {tsId}"))
                    self.itemTsIdReadList.append(tsId)

            time.sleep(10)
            if inTsSet.isStreamOpen():
                with self._lock:
                    inTsSet.loadAllProperties()  # refresh status for the streaming

    def runMarkerfreeStep(self, tsId: str):
        if tsId not in self.failedItems:
            try:
                logger.info(cyanStr(f'tsId = {tsId}: aligning...'))
                with self._lock:
                    ts = self._getCurrentTs(tsId)
                tsIdPath = self._getExtraPath(tsId)
                makePath(tsIdPath)
                tltFn = os.path.join(tsIdPath, tsId + ".rawtlt")
                ts.generateTltFile(tltFilePath = tltFn)
                cmd = "Markerfree "
                # Input TS
                # TODO: Change for getFirstEnableItem() when TODO tomo is updated
                cmd += "-i %s " % ts.getFirstItem().getFileName()
                # Output MRC
                cmd += "-o %s " % self._getExtraOutFile(tsId, "aligned", MRC_EXT) #self._getExtraPath(tsId, tsId + "_aligned.mrc")
                # Tilt angle file
                cmd += "-a %s " % tltFn
                # Geometry -g offset, tilt axis angle, z-axis offset, 
                # thickness, projection matching reconstruction thickness, 
                # output image downsampling ratio, GPU ID
                offset = 0 #TODO
                taAngle = ts.getAcquisition().getTiltAxisAngle()
                zaOffset = 0 #TODO
                thickness = self.geomThickness.get()
                projThickness = self.geomReconThickness.get()
                dsRatio = self.geomDownsample.get()
                gpuId = 0 #TODO
                cmd += "-g %d,%d,%d,%d,%d,%d,%d " % (offset, taAngle, zaOffset, thickness,
                                                     projThickness, dsRatio, gpuId)
                # The number of images used during the projection matching
                cmd += "-p %d " % self.nProjs.get()
                # -s1 means that an xf file will be generated
                cmd += "-s 1 "

                Plugin.runMarkerfree(self, cmd)
            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> MarkerFree execution failed with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, tsId: str):
        if tsId in self.failedItems:
            self.addToOutFailedSet(tsId)
            return
        try:
            with self._lock:
                ts = self._getCurrentTs(tsId)
                self.createOutTs(ts, self._getInTsSet(returnPointer=True))
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    # --------------------------- I/O functions ------------------------------
    def getOutputFailedSet(self,
                           inputPtr: Pointer) -> SetOfTiltSeries:
        """ Create output set for failed TS. """
        inputSet = inputPtr.get()
        failedTs = getattr(self, OUTPUT_TS_FAILED_NAME, None)
        if failedTs:
                failedTs.enableAppend()
        else:
            logger.info(cyanStr('Create the set of failed TS'))
            failedTs = self._createSetOfTiltSeries(suffix='Failed')
            failedTs.copyInfo(inputSet)
            failedTs.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{OUTPUT_TS_FAILED_NAME: failedTs})
            self._defineSourceRelation(inputPtr, failedTs)

        return failedTs

    def addToOutFailedSet(self,
                          tsId: str) -> None:
        """ Just copy input item to the failed output set. """
        logger.info(cyanStr(f'Failed TS ---> {tsId}'))
        try:
            with self._lock:
                inputSet = self._getInTsSet(returnPointer=True)
                output = self.getOutputFailedSet(inputSet)
                item = self._getCurrentTs(tsId)
                newItem = item.clone()
                newItem.copyInfo(item)
                output.append(newItem)
                
                newItem.copyItems(item)
                newItem.write()

                output.update(newItem)
                output.write()
                self._store(output)
                # Explicitly close the outputs (for streaming)
                output.close()
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the failed output with '
                                f'exception {e}. Skipping... '))

    def getOutputSetOfTS(self, inPointer: Pointer):
        #TODO: implement vaga de mierda
        pass

    def createOutTs(self,
                    ts: TiltSeries,
                    inTsSetPointer: Pointer) -> None:
        tsId = ts.getTsId()
        xfFile = self._getExtraOutFile(tsId, "aligned", XF_EXT)
        if os.path.exists(xfFile) and os.stat(xfFile).st_size != 0:
            tltFile = self.getTltFilePath(tsId)
            aliMatrix = readXfFile(xfFile)
            tiltAngles = self.formatAngleList(tltFile) #TODO: como se adapta esto a esta vaina jejeje?
            # Set of tilt-series
            outTsSet = self.getOutputSetOfTS(inTsSetPointer)  #TODO #TODO #TODO
            # Tilt-series
            outTs = TiltSeries()
            outTs.copyInfo(ts)
            outTs.setAlignment2D()
            outTsSet.append(outTs)
            # Tilt-images
            stackIndex = 0
            for ti in ts.iterItems(orderBy=TiltImage.INDEX_FIELD):
                outTi = TiltImage()
                outTi.copyInfo(ti)
                if ti.isEnabled():
                    tiltAngle, newTransformArray = self._getTrDataEnabled(stackIndex,
                                                                          aliMatrix,
                                                                          tiltAngles)
                    stackIndex += 1
                else:
                    tiltAngle, newTransformArray = self._getTrDataDisabled(ti)
                self._updateTiltImage(ti, outTi, newTransformArray, tiltAngle)  #TODO
                self.setTsOddEven(tsId, outTi, binGenerated=False)  #TODO
                outTs.append(outTi)
            # Data persistence
            outTs.write()
            outTsSet.update(outTs)
            outTsSet.write()
            self._store(outTsSet)
            # Close explicitly the outputs (for streaming)
            self.closeOutputsForStreaming()
        else:
            logger.error(f'tsId = {tsId} -> Output file {xfFile} was not generated or is empty. Skipping... ')

    def closeOutputsForStreaming(self):
        # Explicitly close outputs (needed for streaming)
        for outputName in self._possibleOutputs.keys():
            output = getattr(self, outputName, None)
            if output:
                output.close()

    # --------------------------- INFO functions ------------------------------
    def _validate(self) -> List[str]:
        errorMsg = []
        return errorMsg
    
    def readingOutput(self) -> None:
        outTsSet = getattr(self, self._possibleOutputs.tiltSeries.name, None)
        if outTsSet:
            for item in outTsSet:
                self.itemTsIdReadList.append(item.getTsId())
            self.info(cyanStr(f'TsIds processed: {self.itemTsIdReadList}'))
        else:
            self.info(cyanStr('No tilt-series have been processed yet'))
            
    # --------------------------- UTILS functions -----------------------------

    @staticmethod
    def _getTrDataEnabled(stackIndex: int,
                          alignmentMatrix: np.ndarray,
                          tiltAngleList: List[float]) -> Tuple[float, np.ndarray]:
        newTransform = alignmentMatrix[:, :, stackIndex]
        newTransformArray = np.array(newTransform)
        tiltAngle = float(tiltAngleList[stackIndex])
        return tiltAngle, newTransformArray

    @staticmethod
    def _getTrDataDisabled(ti: TiltImage) -> Tuple[float, np.ndarray]:
        tiltAngle = ti.getTiltAngle()
        if ti.hasTransform():
            newTransformArray = ti.getTransform().getMatrix()
        else:
            newTransformArray = IDENTITY_MATRIX
        return tiltAngle, newTransformArray

    def getTltFilePath(self, tsId):
        return self._getExtraOutFile(tsId, suffix="", ext=RAWTLT_EXT)

    @staticmethod
    def _getOutTsFileName(tsId, suffix=None, ext=MRC_EXT):
        return f'{tsId}_{suffix}.{ext}' if suffix else f'{tsId}.{ext}'

    def _getExtraOutFile(self, tsId, suffix=None, ext=MRC_EXT):
        return self._getExtraPath(tsId,
                                  self._getOutTsFileName(tsId, suffix=suffix, ext=ext))

    def _getCurrentTs(self, tsId: str) -> TiltSeries:
        return self._getInTsSet().getItem(TiltSeries.TS_ID_FIELD, tsId)

    def _getInTsSet(self, returnPointer: bool = False) -> Union[SetOfTiltSeries, Pointer]:
        inTsPointer = getattr(self, IN_TS_SET)
        return inTsPointer if returnPointer else inTsPointer.get()
    
    def _getCurrentItem(self, tsId: str, doLock: bool = True) -> TiltSeries:
        if doLock:
            with self._lock:
                return self._getInTsSet().getItem(TiltSeries.TS_ID_FIELD, tsId)
        else:
            return self._getInTsSet().getItem(TiltSeries.TS_ID_FIELD, tsId)
