//
//  CsoundObject.c
//  Audio-Mosaicing
//
//  Created by Edward Costello on 27/03/2013.
//  Copyright (c) 2013 Edward Costello. All rights reserved.
//

#import "CsoundObject.h"
#import <Accelerate/Accelerate.h>
#import "ConvenienceFunctions.h"

CsoundObject *CsoundObject_new(char *csdPath)
{
    CsoundObject *self = calloc(1, sizeof(CsoundObject));
    
    self->increment = 8;
    self->csound = csoundCreate(NULL);
    
    csoundSetHostImplementedAudioIO(self->csound, 1, 0);
    csoundSetHostData(self->csound, self);
    self->channelsCount = 2;

    char *argv[2] = {
        "csound2",
        csdPath
    };
    
    self->csoundResult = csoundCompile(self->csound, 2, argv);
    
    if (self->csoundResult == 0) {
        printf("success\n");
    }
    
    self->kalmanFilter = KalmanFilter_new(44100, self->increment);

    for (size_t i = 0; i < 4; ++i) {
        
        CsoundObject_createTable(self, 100 + i, self->kalmanFilter->filCallbackBuffer->rowCount);
    }
    CsoundObject_turnOnInstruments(self);

    return self;
}

void CsoundObject_delete(CsoundObject *self)
{
    csoundStop(self->csound);
    csoundDestroy(self->csound);
    KalmanFilter_delete(self->kalmanFilter);
    free(self);
    self = NULL;
}

void CsoundObject_createTable(CsoundObject *self, size_t tableNumber, size_t tableSize)
{
    self->csoundScore[0] = tableNumber;
    self->csoundScore[1] = 0;
    self->csoundScore[2] = -((Float64)tableSize);
    self->csoundScore[3] = 2;
    self->csoundScore[4] = 0;
    
    csoundScoreEvent(self->csound, 'f', self->csoundScore, 5);
    csoundPerformKsmps(self->csound);
}

void CsoundObject_writeTables(CsoundObject *self)
{
    size_t tableSize = self->kalmanFilter->filCallbackBuffer->rowCount;
    Matrix *filCallbackBuffer = self->kalmanFilter->filCallbackBuffer;
    Float64 *tablePointer;

    for (size_t i = 0; i < 4; ++i) {

        csoundGetTable(self->csound, &tablePointer, 100 + i);
        
        cblas_dcopy((UInt32)tableSize, &filCallbackBuffer->data[i], (UInt32)self->kalmanFilter->points, tablePointer, 1);        
    }
}

void CsoundObject_turnOnInstruments(CsoundObject *self)
{
    self->csoundScore[0] = 1;
    self->csoundScore[1] = 0;
    self->csoundScore[2] = -((Float64)1);

    csoundScoreEvent(self->csound, 'i', self->csoundScore, 3);
    
    self->csoundScore[0] = 2;
    self->csoundScore[1] = 0;
    self->csoundScore[2] = -((Float64)1);
    self->csoundScore[3] = self->increment;

    csoundScoreEvent(self->csound, 'i', self->csoundScore, 4);
}

void CsoundObject_readCallback(void *inRefCon,
                               UInt32 inNumberFrames,
                               Float64 **audioData)
{
    CsoundObject *self = (CsoundObject *) inRefCon;
    Float64 *csoundOut = (Float64 *) csoundGetSpout(self->csound);
        
    csoundPerformKsmps(self->csound);
    
    for (size_t channel = 0; channel < self->channelsCount; ++channel) {
        
        cblas_dcopy((SInt32)inNumberFrames, &csoundOut[channel], 2, audioData[channel], 1);
    }
}

void CsoundObject_processingCallback(void *inRefCon,
                                     UInt32 inNumberFrames,
                                     Float64 **audioData)
{
    CsoundObject *self = (CsoundObject *) inRefCon;
    Float64 *csoundOut = (Float64 *) csoundGetSpout(self->csound);
    
    KalmanFilter_process(self->kalmanFilter, inNumberFrames, audioData[0]);

    CsoundObject_writeTables(self);

    csoundPerformKsmps(self->csound);

    for (size_t channel = 0; channel < self->channelsCount; ++channel) {
        
        cblas_dcopy((SInt32)inNumberFrames, &csoundOut[channel], 2, audioData[channel], 1);
    }
}

