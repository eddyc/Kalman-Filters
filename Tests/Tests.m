//
//  Tests.m
//  Tests
//
//  Created by Edward Costello on 30/06/2013.
//  Copyright (c) 2013 Edward Costello. All rights reserved.
//

#import "Tests.h"
#import "vowels.h"
#import "Matrix.h"
#import "AudioStream.h"
#import "AudioObject.h"
#import "CsoundObject.h"
#import "RingBuffer.h"
#import <Accelerate/Accelerate.h>

@implementation Tests

- (void)setUp
{
    [super setUp];
    
    // Set-up code here.
}

- (void)tearDown
{
    // Tear-down code here.
    
    [super tearDown];
}

void filter(Matrix *input)
{
    Float64 bWeights[2] = {1, -0.98};
    Float64 aWeight = 1;
    Float64 filterInputVector[2] = {0, 0};
    Float64 filterFeedbackVector[2] = {0, 0};
    Float64 feedforward, feedback;
    
    for (size_t i = 0; i < input->elementCount; ++i) {
        
        feedback = feedforward = 0;
        filterInputVector[0] = input->data[i];
        
        for (int j = 0; j < 2; j++) {
            
            feedforward += bWeights[j] * filterInputVector[j];
            feedback += aWeight * filterFeedbackVector[j];
        }
        
        input->data[i] = feedforward - feedback;
        
        filterFeedbackVector[1] = filterFeedbackVector[0];
        filterInputVector[1] = filterInputVector[0];
    }
}


void scaledSymmetricSigmaPoints(Matrix *x,
                                Matrix *P,
                                UInt32 alpha,
                                UInt32 beta,
                                UInt32 kappa,
                                Matrix *xPts,
                                Matrix *wPts,
                                UInt32 *nPts,
                                Matrix *Psqrtm)
{
    Matrix_clear(xPts);
    
    size_t n = P->rowCount;
    *nPts = (UInt32) n * 2 + 1;
    kappa = powf(alpha, 2) * (n + kappa) - n;
    Float32 nPlusKappa = n + kappa;
    
    Matrix_scalarMultiply(P, nPlusKappa, Psqrtm);
    
    Matrix_choleskyFactorisation(Psqrtm);
    
    Matrix_subMatrixCopy(Psqrtm, 0, 0, Psqrtm->rowCount, Psqrtm->columnCount, xPts, 0, xPts->columnCount - Psqrtm->columnCount);
    Matrix_negate(Psqrtm, Psqrtm);
    Matrix_subMatrixCopy(Psqrtm, 0, 0, Psqrtm->rowCount, Psqrtm->columnCount, xPts, 0, 1);
    
    for (size_t i = 0; i < xPts->columnCount; ++i) {
        
        vDSP_vaddD(&xPts->data[i], xPts->columnCount, x->data, 1, &xPts->data[i], xPts->columnCount, xPts->rowCount);
    }
    
    
    Matrix_setElementsToScalar(wPts, 0.5);
    wPts->data[0] = kappa;
    wPts->data[wPts->elementCount - 1] = 0;
    
    Float64 divisor = n + kappa;
    
    vDSP_vsdivD(wPts->data, 1, &divisor, wPts->data, 1, wPts->elementCount);
    
    wPts->data[wPts->elementCount - 1] = wPts->data[0] + (1 - powf(alpha, 2)) + beta;
}

void ukfRigollRecursion(Matrix *xhat,
                        size_t Fs,
                        Matrix *ahat,
                        Matrix *c,
                        Matrix *d,
                        Matrix *A,
                        Matrix *b)
{
    Matrix_clear(ahat);
    size_t M = xhat->elementCount  / 2;
    
    Float64 reciprocalFs = 1./Fs;
    
    Float64 twoPi = (M_PI) * 2;
    
    Float64 scalar = twoPi * reciprocalFs;
    
    vDSP_vsmulD(xhat->data, 1, &scalar, d->data, 1, M);
    
    Matrix_cosine(d, d);
    
    
    scalar = -M_PI * reciprocalFs;
    
    vDSP_vsmulD(&xhat->data[M], 1, &scalar, c->data, 1, M);
    
    Matrix_exponent(c, c);
    
    vDSP_vmulD(c->data, 1, d->data, 1, c->data, 1, M);
    
    scalar = -2;
    
    Matrix_scalarMultiply(c, scalar, c);
    
    scalar = -2 * M_PI * reciprocalFs;
    vDSP_vsmulD(&xhat->data[M], 1, &scalar, d->data, 1, M);
    
    Matrix_exponent(d, d);
    
    //    Matrix_clear(b);
    Matrix_getRow(b, 0)[0] = c->data[0];
    Matrix_getRow(b, 1)[0] = d->data[0];
    
    cblas_dcopy((UInt32)b->columnCount, b->data, (UInt32)b->columnCount, ahat->data, 1);
    
    Float64 temp[8] = {0};
    
    for (size_t i = 1; i < M; ++i) {
        
        Matrix_identity(A);
        
        Matrix_getRow(A, 1)[0] = c->data[i];
        Matrix_getRow(A, A->rowCount - 1)[A->columnCount - 2] = c->data[i];
        Matrix_getRow(A, A->rowCount - 1)[A->columnCount - 3] = d->data[i];
        
        for (size_t j = 2; j < xhat->elementCount - 1; ++j) {
            
            Matrix_getRow(A, j)[j - 1] = c->data[i];
            Matrix_getRow(A, j)[j - 2] = d->data[i];
            
        }
        
        Matrix_getRow(b, 0)[i] = c->data[i];
        Matrix_getRow(b, 1)[i] = d->data[i];
        vDSP_mmulD(A->data, 1, ahat->data, 1, temp, 1, A->rowCount, 1, A->columnCount);
        vDSP_vaddD(temp, 1, &b->data[i], b->columnCount, ahat->data, 1, xhat->elementCount);
    }
}


- (void)DONTtestKalmanFilters
{
    size_t Fs = 8000;
    size_t samplesCount = 8000;
    size_t points = 8;
    size_t noises = points + 1;
    
    Matrix *speech = Matrix_newWithArray(1, samplesCount, vowels);
    Matrix_normalise(speech, speech);
    filter(speech);
    
    Matrix *PQ = Matrix_newIdentity(points + points + 1);
    Float64 hundred = 100.;
    Float64 fifty = 50.;
    UInt32 nsp;
    vDSP_vsmulD(PQ->data, 1, &hundred, PQ->data, 1, (points + points) * points);
    vDSP_vsmulD(Matrix_getRow(PQ, points), 1, &fifty, Matrix_getRow(PQ, points), 1, (points + points) * (points + 1));
    PQ->data[PQ->elementCount - 1] = 0.1;
    
    Matrix *fil = Matrix_newWithValues(samplesCount, points, 8, 300., 2200., 3000., 3500., 20., 20., 20., 20.);
    Matrix *pred = Matrix_new(samplesCount, points);
    Matrix *xQ = Matrix_newWithValues(1, points + noises, 8, 300., 2200., 3000., 3500., 20., 20., 20., 20.);
    Matrix *wSigmaPts = Matrix_new(1, PQ->rowCount * 2 + 2);
    Matrix *wSigmaPts_xmat = Matrix_new(points, PQ->rowCount * 2);
    Matrix *wSigmaPts_zmat = Matrix_new(1, PQ->rowCount * 2);
    Matrix *xSigmaPts = Matrix_new(PQ->rowCount, PQ->rowCount * 2 + 1);
    Matrix *Psqrtm = Matrix_new(PQ->rowCount, PQ->rowCount);
    Matrix *xPredSigmaPts = Matrix_new(points, PQ->rowCount * 2 + 1);
    Matrix *v = Matrix_new(1, points);
    Matrix *A = Matrix_new(points, points);
    Matrix *b = Matrix_new(points, points / 2);
    Matrix *c = Matrix_new(1, points / 2);
    Matrix *d = Matrix_new(1, points / 2);
    Matrix *ahat = Matrix_new(1, points);
    Matrix *xhat = Matrix_new(1, points);
    Matrix *zPredSigmaPts = Matrix_new(1, points * 4 + 3);
    Matrix *temp_8by34 = Matrix_new(8, 34);
    Matrix *exSigmaPt2 = Matrix_new(8, 34);
    Matrix *temp_1by34 = Matrix_new(1, 34);
    Matrix *ezSigmaPt2 = Matrix_new(1, 34);
    Matrix *temp_8by8 = Matrix_new(8, 8);
    Matrix *temp_1by8 = Matrix_new(1, 8);
    Matrix *xPred = Matrix_new(1, points);
    Matrix *exSigmaPt = Matrix_new(1, points);
    Matrix *PPred = Matrix_new(8, 8);
    Matrix *PxzPred = Matrix_new(1, points);
    Matrix *K = Matrix_new(1, 8);
    
    
    
    
    for (size_t i = 0; i < samplesCount - points; ++i) {
        
        scaledSymmetricSigmaPoints(xQ, PQ, 1., 0., 2., xSigmaPts, wSigmaPts, &nsp, Psqrtm);
        
        //        Matrix_saveAsHDF(xQ, "desktop", "/xQ");
        //        Matrix_saveAsHDF(PQ, "desktop", "/PQ");
        //        Matrix_saveAsHDF(xSigmaPts, "desktop", "/xSigmaPts");
        //        Matrix_saveAsHDF(wSigmaPts, "desktop", "/wSigmaPts");
        //
        //
        vDSP_vaddD(xSigmaPts->data, 1, Matrix_getRow(xSigmaPts, points), 1, xPredSigmaPts->data, 1, xPredSigmaPts->elementCount);
        
        cblas_dcopy((SInt32)points, &speech->data[i], 1, v->data, 1);
        Matrix_reverse(v);
        Matrix_negate(v, v);
        
        for (size_t j = 0; j < nsp; ++j) {
            cblas_dcopy((UInt32)xhat->elementCount, &xPredSigmaPts->data[j], (UInt32)xPredSigmaPts->columnCount, xhat->data, 1);
            
            ukfRigollRecursion(xhat, Fs, ahat, c, d, A, b);
            
            Float64 dotProduct;
            vDSP_dotprD(ahat->data, 1, v->data, 1, &dotProduct, points);
            zPredSigmaPts->data[j] = dotProduct + Matrix_getRow(xSigmaPts, points + points)[j];
        }
        
        for (size_t j = 0; j < nsp - 1; ++j) {
            
            vDSP_vsubD(xPredSigmaPts->data, xPredSigmaPts->columnCount, &xPredSigmaPts->data[j + 1], xPredSigmaPts->columnCount, &temp_8by34->data[j], temp_8by34->columnCount, temp_8by34->rowCount);
        }
        
        Matrix_elementWiseMultiply(wSigmaPts_xmat, temp_8by34, temp_8by34);
        for (size_t j = 0; j < points; ++j) {
            
            vDSP_sveD(Matrix_getRow(temp_8by34, j), 1, &Matrix_getRow(pred, i)[j], temp_8by34->columnCount);
        }
        
        vDSP_vaddD(Matrix_getRow(pred, i), 1, xPredSigmaPts->data, xPredSigmaPts->columnCount, Matrix_getRow(pred, i), 1, xPred->elementCount);
        
        Float64 scalar = -zPredSigmaPts->data[0];
        vDSP_vsaddD(&zPredSigmaPts->data[1], 1, &scalar, temp_1by34->data, 1, temp_1by34->elementCount);
        vDSP_vmulD(&wSigmaPts->data[1], 1, temp_1by34->data, 1, temp_1by34->data, 1, temp_1by34->elementCount);
        
        Float64 zPred;
        vDSP_sveD(temp_1by34->data, 1, &zPred, temp_1by34->elementCount);
        
        zPred += zPredSigmaPts->data[0];
        
        vDSP_vsubD(Matrix_getRow(pred, i), 1, xPredSigmaPts->data, xPredSigmaPts->columnCount, exSigmaPt->data, 1, exSigmaPt->elementCount);
        
        Float32 ezSigmaPt = zPredSigmaPts->data[0] - zPred;
        
        Matrix_multiply(exSigmaPt, true, exSigmaPt, false, temp_8by8);
        Matrix_scalarMultiply(temp_8by8, wSigmaPts->data[wSigmaPts->elementCount - 1], PPred);
        
        Matrix_scalarMultiply(exSigmaPt, ezSigmaPt, PxzPred);
        Matrix_scalarMultiply(PxzPred, wSigmaPts->data[wSigmaPts->elementCount - 1], PxzPred);
        
        
        Float64 S = wSigmaPts->data[wSigmaPts->elementCount - 1] * ezSigmaPt * ezSigmaPt;
        
        for (size_t j = 0; j < exSigmaPt2->columnCount; ++j) {
            
            vDSP_vsubD(Matrix_getRow(pred, i), 1, &xPredSigmaPts->data[j + 1], xPredSigmaPts->columnCount, &exSigmaPt2->data[j], exSigmaPt2->columnCount, exSigmaPt2->rowCount);
        }
        
        scalar = -zPred;
        vDSP_vsaddD(&zPredSigmaPts->data[1], 1, &scalar, ezSigmaPt2->data, 1, ezSigmaPt2->elementCount);
        
        for (size_t j = 0; j < exSigmaPt2->rowCount; ++j) {
            
            vDSP_vmulD(&wSigmaPts->data[1], 1, Matrix_getRow(exSigmaPt2, j), 1, Matrix_getRow(temp_8by34, j), 1, exSigmaPt2->columnCount);
        }
        
        Matrix_multiply(temp_8by34, false, exSigmaPt2, true, temp_8by8);
        Matrix_add(temp_8by8, PPred, PPred);
        vDSP_vmulD(&wSigmaPts->data[1], 1, ezSigmaPt2->data, 1, temp_1by34->data, 1, temp_1by34->elementCount);
        
        Float64 dotProduct;
        
        vDSP_dotprD(temp_1by34->data, 1, ezSigmaPt2->data, 1, &dotProduct, temp_1by34->elementCount);
        S += dotProduct;
        
        Matrix_multiply(temp_1by34, false, exSigmaPt2, true, temp_1by8);
        vDSP_vaddD(PxzPred->data, 1, temp_1by8->data, 1, PxzPred->data, 1, PxzPred->elementCount);
        Matrix_scalarDivide(PxzPred, S, K);
        
        Matrix_scalarMultiply(K, S, temp_1by8);
        Matrix_multiply(K, true, temp_1by8, false, temp_8by8);
        
        for (size_t j = 0; j < points; ++j) {
            
            vDSP_vsubD(Matrix_getRow(temp_8by8, j), 1, Matrix_getRow(PPred, j), 1, Matrix_getRow(PQ, j), 1, points);
        }
        
        Float32 inovation = speech->data[i + points] - zPred;
        Matrix_scalarMultiply(K, inovation, temp_1by8);
        vDSP_vaddD(Matrix_getRow(pred, i), 1, temp_1by8->data, 1, Matrix_getRow(fil, i), 1, temp_1by8->elementCount);
        
        cblas_dcopy((SInt32)points, Matrix_getRow(fil, i), 1, xQ->data, 1);
    }
    
    Matrix_saveAsHDF(fil, "desktop", "/fil");
    
    Matrix_delete(A);
    Matrix_delete(b);
    Matrix_delete(c);
    Matrix_delete(d);
    Matrix_delete(K);
    Matrix_delete(Psqrtm);
    Matrix_delete(PxzPred);
    Matrix_delete(wSigmaPts);
    Matrix_delete(wSigmaPts_xmat);
    Matrix_delete(wSigmaPts_zmat);
    Matrix_delete(xSigmaPts);
    Matrix_delete(PPred);
    Matrix_delete(xPredSigmaPts);
    Matrix_delete(zPredSigmaPts);
    Matrix_delete(temp_8by34);
    Matrix_delete(exSigmaPt);
    Matrix_delete(exSigmaPt2);
    Matrix_delete(temp_1by34);
    Matrix_delete(ezSigmaPt2);
    Matrix_delete(temp_8by8);
    Matrix_delete(temp_1by8);
    Matrix_delete(speech);
    Matrix_delete(PQ);
    Matrix_delete(fil);
    Matrix_delete(pred);
    Matrix_delete(xQ);
}

- (void)testAudioStream
{
    AudioObject *audioObject = AudioObject_new();
    AudioObject_openAudioFile(audioObject, "/Users/eddyc/Desktop/Test.wav");
    CsoundObject *csoundObject = CsoundObject_new("/Users/eddyc/Dropbox/SoftwareProjects/Kalman Filters/Kalman Filters/CsoundObject.csd");
    AudioStream *audioStream = AudioStream_new(audioObject,
                                               AudioObject_readCallbackMono,
                                               csoundObject,
                                               CsoundObject_processingCallback);
     
    AudioStream_process(audioStream, 10);
    
    AudioObject_delete(audioObject);
    AudioStream_delete(audioStream);
}

- (void)DONTtestKalmanFilterStream
{
    size_t Fs = 8000;
    size_t samplesCount = 8000;
    size_t points = 8;
    size_t noises = points + 1;
    
    Matrix *speech = Matrix_newWithArray(1, samplesCount, vowels);
    Matrix_normalise(speech, speech);
    filter(speech);
    
    Matrix *PQ = Matrix_newIdentity(points + points + 1);
    Float64 hundred = 100.;
    Float64 fifty = 50.;
    UInt32 nsp;
    vDSP_vsmulD(PQ->data, 1, &hundred, PQ->data, 1, (points + points) * points);
    vDSP_vsmulD(Matrix_getRow(PQ, points), 1, &fifty, Matrix_getRow(PQ, points), 1, (points + points) * (points + 1));
    PQ->data[PQ->elementCount - 1] = 0.1;
    
    size_t inNumberFrames = 512;

    Matrix *fil = Matrix_new(samplesCount, points);
    Matrix *filBuffer = Matrix_new(1, points);
    Matrix *callbackFil = Matrix_new(inNumberFrames, points);

    Matrix *pred = Matrix_new(1, points);
    Matrix *xQ = Matrix_newWithValues(1, points + noises, 8, 300., 2200., 3000., 3500., 20., 20., 20., 20.);
    Matrix *wSigmaPts = Matrix_new(1, PQ->rowCount * 2 + 2);
    Matrix *wSigmaPts_xmat = Matrix_new(points, PQ->rowCount * 2);
    Matrix *wSigmaPts_zmat = Matrix_new(1, PQ->rowCount * 2);
    Matrix *xSigmaPts = Matrix_new(PQ->rowCount, PQ->rowCount * 2 + 1);
    Matrix *Psqrtm = Matrix_new(PQ->rowCount, PQ->rowCount);
    Matrix *xPredSigmaPts = Matrix_new(points, PQ->rowCount * 2 + 1);
    Matrix *v = Matrix_new(1, points);
    Matrix *A = Matrix_new(points, points);
    Matrix *b = Matrix_new(points, points / 2);
    Matrix *c = Matrix_new(1, points / 2);
    Matrix *d = Matrix_new(1, points / 2);
    Matrix *ahat = Matrix_new(1, points);
    Matrix *xhat = Matrix_new(1, points);
    Matrix *zPredSigmaPts = Matrix_new(1, points * 4 + 3);
    Matrix *temp_8by34 = Matrix_new(8, 34);
    Matrix *exSigmaPt2 = Matrix_new(8, 34);
    Matrix *temp_1by34 = Matrix_new(1, 34);
    Matrix *ezSigmaPt2 = Matrix_new(1, 34);
    Matrix *temp_8by8 = Matrix_new(8, 8);
    Matrix *temp_1by8 = Matrix_new(1, 8);
    Matrix *xPred = Matrix_new(1, points);
    Matrix *exSigmaPt = Matrix_new(1, points);
    Matrix *PPred = Matrix_new(8, 8);
    Matrix *PxzPred = Matrix_new(1, points);
    Matrix *K = Matrix_new(1, 8);
    

    Float64 *ioData = calloc(inNumberFrames, sizeof(Float64));
    RingBuffer *ringBuffer = RingBuffer_new(16);
    Float64 *buffer = calloc(16, sizeof(Float64));

    size_t bufferPointer = 0;
    for (size_t frame = 0; frame < (8000 / inNumberFrames); ++frame) {
        
        cblas_dcopy((SInt32)inNumberFrames, &speech->data[frame * inNumberFrames], 1, ioData, 1);
        
        if (frame == 0) {
            
            RingBuffer_write(ringBuffer, ioData, ringBuffer->capacity);
        }
        
        for (size_t i = (frame == 0 ? ringBuffer->capacity : 0); i < inNumberFrames; ++i, ++bufferPointer) {
            
            scaledSymmetricSigmaPoints(xQ, PQ, 1., 0., 2., xSigmaPts, wSigmaPts, &nsp, Psqrtm);
  
            vDSP_vaddD(xSigmaPts->data, 1, Matrix_getRow(xSigmaPts, points), 1, xPredSigmaPts->data, 1, xPredSigmaPts->elementCount);
            
            RingBuffer_copy(ringBuffer, buffer, points + 1);
            RingBuffer_skip(ringBuffer, 1);
            RingBuffer_write(ringBuffer, &ioData[i], 1);
            
            cblas_dcopy((SInt32)points, buffer, 1, v->data, 1);
//            cblas_dcopy((SInt32)points, &speech->data[bufferPointer], 1, v->data, 1);

            Matrix_reverse(v);
            Matrix_negate(v, v);
            
            for (size_t j = 0; j < nsp; ++j) {
                cblas_dcopy((UInt32)xhat->elementCount, &xPredSigmaPts->data[j], (UInt32)xPredSigmaPts->columnCount, xhat->data, 1);
                
                ukfRigollRecursion(xhat, Fs, ahat, c, d, A, b);
                
                Float64 dotProduct;
                vDSP_dotprD(ahat->data, 1, v->data, 1, &dotProduct, points);
                zPredSigmaPts->data[j] = dotProduct + Matrix_getRow(xSigmaPts, points + points)[j];
            }
            
            for (size_t j = 0; j < nsp - 1; ++j) {
                
                vDSP_vsubD(xPredSigmaPts->data, xPredSigmaPts->columnCount, &xPredSigmaPts->data[j + 1], xPredSigmaPts->columnCount, &temp_8by34->data[j], temp_8by34->columnCount, temp_8by34->rowCount);
            }
            
            Matrix_elementWiseMultiply(wSigmaPts_xmat, temp_8by34, temp_8by34);
            for (size_t j = 0; j < points; ++j) {
                
                vDSP_sveD(Matrix_getRow(temp_8by34, j), 1, &pred->data[j], temp_8by34->columnCount);
            }
            
            vDSP_vaddD(pred->data, 1, xPredSigmaPts->data, xPredSigmaPts->columnCount, pred->data, 1, xPred->elementCount);
            
            Float64 scalar = -zPredSigmaPts->data[0];
            vDSP_vsaddD(&zPredSigmaPts->data[1], 1, &scalar, temp_1by34->data, 1, temp_1by34->elementCount);
            vDSP_vmulD(&wSigmaPts->data[1], 1, temp_1by34->data, 1, temp_1by34->data, 1, temp_1by34->elementCount);
            
            Float64 zPred;
            vDSP_sveD(temp_1by34->data, 1, &zPred, temp_1by34->elementCount);
            
            zPred += zPredSigmaPts->data[0];
            
            vDSP_vsubD(pred->data, 1, xPredSigmaPts->data, xPredSigmaPts->columnCount, exSigmaPt->data, 1, exSigmaPt->elementCount);
            
            Float32 ezSigmaPt = zPredSigmaPts->data[0] - zPred;
            
            Matrix_multiply(exSigmaPt, true, exSigmaPt, false, temp_8by8);
            Matrix_scalarMultiply(temp_8by8, wSigmaPts->data[wSigmaPts->elementCount - 1], PPred);
            
            Matrix_scalarMultiply(exSigmaPt, ezSigmaPt, PxzPred);
            Matrix_scalarMultiply(PxzPred, wSigmaPts->data[wSigmaPts->elementCount - 1], PxzPred);
            
            
            Float64 S = wSigmaPts->data[wSigmaPts->elementCount - 1] * ezSigmaPt * ezSigmaPt;
            
            for (size_t j = 0; j < exSigmaPt2->columnCount; ++j) {
                
                vDSP_vsubD(pred->data, 1, &xPredSigmaPts->data[j + 1], xPredSigmaPts->columnCount, &exSigmaPt2->data[j], exSigmaPt2->columnCount, exSigmaPt2->rowCount);
            }
            
            scalar = -zPred;
            vDSP_vsaddD(&zPredSigmaPts->data[1], 1, &scalar, ezSigmaPt2->data, 1, ezSigmaPt2->elementCount);
            
            for (size_t j = 0; j < exSigmaPt2->rowCount; ++j) {
                
                vDSP_vmulD(&wSigmaPts->data[1], 1, Matrix_getRow(exSigmaPt2, j), 1, Matrix_getRow(temp_8by34, j), 1, exSigmaPt2->columnCount);
            }
            
            Matrix_multiply(temp_8by34, false, exSigmaPt2, true, temp_8by8);
            Matrix_add(temp_8by8, PPred, PPred);
            vDSP_vmulD(&wSigmaPts->data[1], 1, ezSigmaPt2->data, 1, temp_1by34->data, 1, temp_1by34->elementCount);
            
            Float64 dotProduct;
            
            vDSP_dotprD(temp_1by34->data, 1, ezSigmaPt2->data, 1, &dotProduct, temp_1by34->elementCount);
            S += dotProduct;
            
            Matrix_multiply(temp_1by34, false, exSigmaPt2, true, temp_1by8);
            vDSP_vaddD(PxzPred->data, 1, temp_1by8->data, 1, PxzPred->data, 1, PxzPred->elementCount);
            Matrix_scalarDivide(PxzPred, S, K);
            
            Matrix_scalarMultiply(K, S, temp_1by8);
            Matrix_multiply(K, true, temp_1by8, false, temp_8by8);
            
            for (size_t j = 0; j < points; ++j) {
                
                vDSP_vsubD(Matrix_getRow(temp_8by8, j), 1, Matrix_getRow(PPred, j), 1, Matrix_getRow(PQ, j), 1, points);
            }
            
            Float32 inovation = buffer[points] - zPred;
//            Float32 inovation = speech->data[bufferPointer + points] - zPred;
            Matrix_scalarMultiply(K, inovation, temp_1by8);
            vDSP_vaddD(pred->data, 1, temp_1by8->data, 1, filBuffer->data, 1, temp_1by8->elementCount);
            
            cblas_dcopy((SInt32)points, filBuffer->data, 1, xQ->data, 1);
            cblas_dcopy((SInt32)points, filBuffer->data, 1, Matrix_getRow(callbackFil, i), 1);

//            printf("b,[%zd] %f\ns,[%zd] %f\n",bufferPointer, buffer[0], bufferPointer, speech->data[bufferPointer]);
        }
        
        cblas_dcopy((SInt32)callbackFil->elementCount, callbackFil->data, 1, Matrix_getRow(fil, inNumberFrames * frame), 1);

    }
    
    
    Matrix_saveAsHDF(fil, "desktop", "/callbackFil");
    
    Matrix_delete(A);
    Matrix_delete(b);
    Matrix_delete(c);
    Matrix_delete(d);
    Matrix_delete(K);
    Matrix_delete(Psqrtm);
    Matrix_delete(PxzPred);
    Matrix_delete(wSigmaPts);
    Matrix_delete(wSigmaPts_xmat);
    Matrix_delete(wSigmaPts_zmat);
    Matrix_delete(xSigmaPts);
    Matrix_delete(PPred);
    Matrix_delete(xPredSigmaPts);
    Matrix_delete(zPredSigmaPts);
    Matrix_delete(temp_8by34);
    Matrix_delete(exSigmaPt);
    Matrix_delete(exSigmaPt2);
    Matrix_delete(temp_1by34);
    Matrix_delete(ezSigmaPt2);
    Matrix_delete(temp_8by8);
    Matrix_delete(temp_1by8);
    Matrix_delete(speech);
    Matrix_delete(PQ);
    Matrix_delete(fil);
    Matrix_delete(pred);
    Matrix_delete(xQ);
}

- (void)DONTtestKalmanStream
{
    size_t samplesCount = 8000;
    Matrix *speech = Matrix_newWithArray(1, samplesCount, vowels);
    Matrix_normalise(speech, speech);
    filter(speech);
    
    size_t increment = 2;
    KalmanFilter *kalmanFilter = KalmanFilter_new(samplesCount, 2);
    size_t inNumberFrames = 512;
    Matrix *fil = Matrix_new(samplesCount / increment, kalmanFilter->points);
    
    
    for (size_t frame = 0; frame < (samplesCount / inNumberFrames); ++frame) {
        
        KalmanFilter_process(kalmanFilter, inNumberFrames, &speech->data[frame * inNumberFrames]);
        
//        Matrix_print(kalmanFilter->filCallbackBuffer);
        cblas_dcopy((UInt32)kalmanFilter->filCallbackBuffer->elementCount, kalmanFilter->filCallbackBuffer->data, 1, Matrix_getRow(fil, (inNumberFrames / increment) * frame), 1);
    }
    
    Matrix_saveAsHDF(fil, "desktop", "/fil");
    Matrix_delete(fil);
}

@end
