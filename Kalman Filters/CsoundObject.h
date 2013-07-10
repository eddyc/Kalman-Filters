//
//  CsoundObject.h
//  Audio-Mosaicing
//
//  Created by Edward Costello on 27/03/2013.
//  Copyright (c) 2013 Edward Costello. All rights reserved.
//


#import <MacTypes.h>
#import <stdio.h>
#import <stdlib.h>
#import <CsoundLib64/csound.h>
#import "Matrix.h"
#import "KalmanFilter.h"

#ifdef __cplusplus
extern "C"
{
#endif
    
    typedef struct CsoundObject
    {
        CSOUND *csound;
        SInt32 csoundResult;
        UInt32 channelsCount;
        Float64 csoundScore[10];
        KalmanFilter *kalmanFilter;
        size_t increment;
        
    } CsoundObject;
    
    CsoundObject *CsoundObject_new(char *csdPath);
    void CsoundObject_delete(CsoundObject *self);

    void CsoundObject_createTable(CsoundObject *self, size_t tableNumber, size_t tableSize);
    void CsoundObject_writeTables(CsoundObject *self);
    void CsoundObject_turnOnInstruments(CsoundObject *self);

    void CsoundObject_readCallback(void *inRefCon,
                                   UInt32 inNumberFrames,
                                   Float64 **audioData);
    void CsoundObject_processingCallback(void *inRefCon,
                                         UInt32 inNumberFrames,
                                         Float64 **audioData);
    #ifdef __cplusplus
}
#endif