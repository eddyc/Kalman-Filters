//
//  KalmanFilter.h
//  Kalman Filters
//
//  Created by Edward Costello on 09/07/2013.
//  Copyright (c) 2013 Edward Costello. All rights reserved.
//


#import <MacTypes.h>
#import <stdio.h>
#import <stdlib.h>
#import "Matrix.h"
#import "RingBuffer.h"

#ifdef __cplusplus
extern "C"
{
#endif
    
    typedef struct KalmanFilter
    {
        Float64 samplerate;
        size_t points;
        size_t noises;
        
        Matrix *PQ;
        Matrix *pred;
        Matrix *xQ;
        Matrix *wSigmaPts;
        Matrix *wSigmaPts_xmat;
        Matrix *wSigmaPts_zmat;
        Matrix *xSigmaPts;
        Matrix *Psqrtm;
        Matrix *xPredSigmaPts;
        Matrix *v;
        Matrix *A;
        Matrix *b;
        Matrix *c;
        Matrix *d;
        Matrix *ahat;
        Matrix *xhat;
        Matrix *zPredSigmaPts;
        Matrix *temp_8by34;
        Matrix *exSigmaPt2;
        Matrix *temp_1by34;
        Matrix *ezSigmaPt2;
        Matrix *temp_8by8;
        Matrix *temp_1by8;
        Matrix *xPred;
        Matrix *exSigmaPt;
        Matrix *PPred;
        Matrix *PxzPred;
        Matrix *K;
        Matrix *filBuffer;
        Matrix *filCallbackBuffer;

        RingBuffer *ringBuffer;
        Float64 *buffer;
        size_t frame;
        size_t increment;
        UInt32 nsp;
        Float64 filterInputVector[2];
        Float64 filterFeedbackVector[2];

    } KalmanFilter;
    
    KalmanFilter *KalmanFilter_new(size_t samplerate, size_t increment);
    void KalmanFilter_delete(KalmanFilter *self);
    
    void KalmanFilter_process(KalmanFilter *self, size_t inNumberFrames, Float64 *audioIn);

    void KalmanFilter_filter(KalmanFilter *self, size_t inNumberFrames, Float64 *input);

    void KalmanFilter_rigollRecursion(Matrix *xhat,
                                      Float64 Fs,
                                      Matrix *ahat,
                                      Matrix *c,
                                      Matrix *d,
                                      Matrix *A,
                                      Matrix *b);
    
    void KalmanFilter_scaledSymmetricSigmaPoints(Matrix *x,
                                                 Matrix *P,
                                                 UInt32 alpha,
                                                 UInt32 beta,
                                                 UInt32 kappa,
                                                 Matrix *xPts,
                                                 Matrix *wPts,
                                                 UInt32 *nPts,
                                                 Matrix *Psqrtm);
    
    
#ifdef __cplusplus
}
#endif