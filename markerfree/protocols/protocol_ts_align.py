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
from os import stat
from os.path import exists

from collections import Counter
from enum import Enum
from typing import List, Tuple, Union
import numpy as np

from pwem.emlib import DT_FLOAT
from pwem.protocols import EMProtocol
from pwem.objects.data import Transform
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
FAILED_TS = 'FailedTiltSeries'

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

        form.addParam('geomOffset', params.FloatParam, expertLevel=LEVEL_ADVANCED, default=0.0, label="Offset")
        form.addParam('geomZAxisOffset', params.FloatParam, expertLevel=LEVEL_ADVANCED, default=0.0, label="Z axis offset")
        form.addParam('geomThickness', params.IntParam, default=200, label="Thickness")
        form.addParam('geomReconThickness', params.IntParam, default=300, label="Reconstruction thickness")
        form.addParam('geomDownsample', params.FloatParam, default=1.0, label="Downsample factor")

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
                    #TODO: Add exclude views with a convertInputStep
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
                ts = self.getTsFromTsId(tsId,doLock=True)
                tsIdPath = self._getExtraPath(tsId)
                makePath(tsIdPath)
                tltFn = os.path.join(tsIdPath, tsId + ".tlt")
                ts.generateTltFile(tltFilePath = tltFn)
                # Input TS
                # TODO: Change for getFirstEnableItem() when TODO tomo is updated
                args = "-i %s " % ts.getFirstItem().getFileName()
                # Output MRC
                args += "-o %s " % self._getExtraOutFile(tsId, "aligned", MRC_EXT) #self._getExtraPath(tsId, tsId + "_aligned.mrc")
                # Tilt angle file
                args += "-a %s " % tltFn
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
                args += "-g %d,%d,%d,%d,%d,%d,%d " % (offset, taAngle, zaOffset, thickness,
                                                     projThickness, dsRatio, gpuId)
                # The number of images used during the projection matching
                args += "-p %d " % self.nProjs.get()
                # -s1 means that an xf file will be generated
                args += "-s 1 "

                Plugin.runMarkerfree(self, args)
            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> MarkerFree execution failed with the exception -> {e}'))
                logger.error(traceback.format_exc())

    def createOutputStep(self, tsId: str):
        if tsId in self.failedItems:
            self.createOutputFailedTs(tsId)
            print('FAILED')
        try:
            self.createOutputTs(tsId)
            print('SUCESSED')
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())
    
    def createOutputTs(self, tsId: str) -> None:
        ts = self.getTsFromTsId(tsId,doLock=True) 
        xfFile = self._getExtraPath(tsId, tsId+'_aligned'+XF_EXT)
        if exists(xfFile) and stat(xfFile).st_size != 0:
            print('entro')
            tltFn = os.path.join(self._getExtraPath(tsId), tsId + ".tlt")
            aliMatrix = readXfFile(xfFile)
            print(aliMatrix)
            tiltAngles = self.formatAngleList(tltFn)
            print(tiltAngles)
            # Set of tilt-series
            outTsSet = self.getOutputSetOfTS(self._getInTsSet())

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
                self._updateTiltImage(ti, outTi, newTransformArray, tiltAngle)
                outTs.append(outTi)
            # Data persistence
            outTs.write()
            outTsSet.update(outTs)
            outTsSet.write()
            self._store(outTsSet)
            # Close explicitly the outputs (for streaming)
            self.closeOutputsForStreaming()

    def closeOutputsForStreaming(self):
        # Close explicitly the outputs (for streaming)
        for outputName in self._possibleOutputs.keys():
            output = getattr(self, outputName, None)
            if output:
                output.close()

    @staticmethod
    def _updateTiltImage(ti: TiltImage,
                         outTi: TiltImage,
                         newTransformArray: np.ndarray,
                         tiltAngle: float) -> None:
        transform = Transform()
        if ti.hasTransform():
            previousTransform = ti.getTransform().getMatrix()
            previousTransformArray = np.array(previousTransform)
            outputTransformMatrix = np.matmul(newTransformArray, previousTransformArray)
            transform.setMatrix(outputTransformMatrix)
        else:
            transform.setMatrix(newTransformArray)

        outTi.setTransform(transform)
        outTi.setTiltAngle(tiltAngle)


    def getOutputSetOfTS(self,
                         inputPtr,
                         attrName=OUTPUT_TILTSERIES_NAME,
                         tiltAxisAngle=None,
                         suffix="") -> SetOfTiltSeries:
        """ Method to generate output of set of tilt-series.
        :param inputPtr: input set pointer (TS or tomograms)
        :param binning: binning factor
        :param attrName: output attr name
        :param tiltAxisAngle: Only applies to TS. If not None, the corresponding value of the
        set acquisition will be updated (xCorr prot)
        :param suffix: output set suffix
        """
        print('entro en getOutputSetOfTS')
        inputSet = inputPtr.get()
        outputSet = getattr(self, attrName, None)
        if outputSet:
            outputSet.enableAppend()
        else:
            outputSet = self._createSetOfTiltSeries(suffix=suffix)

            outputSet.copyInfo(inputSet)
            if tiltAxisAngle:
                outputSet.getAcquisition().setTiltAxisAngle(tiltAxisAngle)

            outputSet.setStreamState(Set.STREAM_OPEN)

            # Write set properties, otherwise it may expose the set (sqlite) without properties.
            outputSet.write()

            self._defineOutputs(**{attrName: outputSet})
            self._defineSourceRelation(inputPtr, outputSet)

        
    def createOutputFailedTs(self, tsId: str):
        logger.info(cyanStr(f'Failed TS ---> {tsId}'))
        try:
            with self._lock:
                ts = self.getTsFromTsId(tsId, doLock=False)
                inTsSet = self._getInTsSet()
                outTsSet = self.getOutputFailedSetOfTiltSeries(inTsSet)
                newTs = TiltSeries()
                newTs.copyInfo(ts)
                outTsSet.append(newTs)
                newTs.copyItems(ts)
                newTs.write()
                outTsSet.update(newTs)
                outTsSet.write()
                self._store(outTsSet)
                # Close explicitly the outputs (for streaming)
                outTsSet.close()
        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the failed output with '
                                f'exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    def getOutputFailedSetOfTiltSeries(self, inputSet):
        failedTsSet = getattr(self, FAILED_TS, None)
        if failedTsSet:
            failedTsSet.enableAppend()
        else:
            failedTsSet = SetOfTiltSeries.create(self._getPath(), template='tiltseries', suffix='Failed')
            failedTsSet.copyInfo(inputSet)
            failedTsSet.setDim(inputSet.getDim())
            failedTsSet.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{FAILED_TS: failedTsSet})
            self._defineSourceRelation(self._getInTsSet(), failedTsSet)

        return failedTsSet
    
    def getTsFromTsId(self,
                      tsId: str,
                      doLock: bool = True) -> TiltSeries:
        tsSet = self._getInTsSet()
        if doLock:
            with self._lock:
                return tsSet.getItem(TiltSeries.TS_ID_FIELD, tsId)
        else:
            return tsSet.getItem(TiltSeries.TS_ID_FIELD, tsId)
    

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

    @staticmethod
    def formatAngleList(tltFilePath):
        """ This method takes an IMOD-based angle file path and
        returns a list containing the angles for each tilt-image
        belonging to the tilt-series. """

        angleList = []

        with open(tltFilePath) as f:
            tltText = f.read().splitlines()
            for line in tltText:
                angleList.append(float(line))

        return angleList

    def getTltFilePath(self, tsId):
        return self._getExtraOutFile(tsId, suffix="", ext=TLT_EXT)

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
