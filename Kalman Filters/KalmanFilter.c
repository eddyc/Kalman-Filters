//
//  KalmanFilter.c
//  Kalman Filters
//
//  Created by Edward Costello on 09/07/2013.
//  Copyright (c) 2013 Edward Costello. All rights reserved.
//

#import "KalmanFilter.h"

KalmanFilter *KalmanFilter_new(size_t samplerate, size_t increment)
{
    KalmanFilter *self = calloc(1, sizeof(KalmanFilter));
    
    self->samplerate = (Float64) samplerate / (Float64) increment;
    self->increment = increment;
    self->points = 8;
    self->noises = self->points + 1;
    
    self->PQ = Matrix_newIdentity(self->points + self->points + 1);
    Float64 hundred = 100.;
    Float64 fifty = 50.;
    self->nsp = 0;
    vDSP_vsmulD(self->PQ->data, 1, &hundred, self->PQ->data, 1, (self->points + self->points) * self->points);
    vDSP_vsmulD(Matrix_getRow(self->PQ, self->points), 1, &fifty, Matrix_getRow(self->PQ, self->points), 1, (self->points + self->points) * (self->points + 1));
    self->PQ->data[self->PQ->elementCount - 1] = 0.1;
    
    self->pred = Matrix_new(1, self->points);
    size_t nyquist = self->samplerate/2;
//    self->xQ = Matrix_newWithValues(1, self->points + self->noises, 8, nyquist * 1./8., nyquist * 3./8., nyquist * 5./8., nyquist * 7./8., 20., 20., 20., 20.);
    
    self->xQ = Matrix_newWithValues(1, self->points + self->noises, 8, 200., 2200., 3000., 3500., 20., 20., 20., 20.);
    self->wSigmaPts = Matrix_new(1, self->PQ->rowCount * 2 + 2);
    self->wSigmaPts_xmat = Matrix_new(self->points, self->PQ->rowCount * 2);
    self->wSigmaPts_zmat = Matrix_new(1, self->PQ->rowCount * 2);
    self->xSigmaPts = Matrix_new(self->PQ->rowCount, self->PQ->rowCount * 2 + 1);
    self->Psqrtm = Matrix_new(self->PQ->rowCount, self->PQ->rowCount);
    self->xPredSigmaPts = Matrix_new(self->points, self->PQ->rowCount * 2 + 1);
    self->v = Matrix_new(1, self->points);
    self->A = Matrix_new(self->points, self->points);
    self->b = Matrix_new(self->points, self->points / 2);
    self->c = Matrix_new(1, self->points / 2);
    self->d = Matrix_new(1, self->points / 2);
    self->ahat = Matrix_new(1, self->points);
    self->xhat = Matrix_new(1, self->points);
    self->zPredSigmaPts = Matrix_new(1, self->points * 4 + 3);
    self->temp_8by34 = Matrix_new(8, 34);
    self->exSigmaPt2 = Matrix_new(8, 34);
    self->temp_1by34 = Matrix_new(1, 34);
    self->ezSigmaPt2 = Matrix_new(1, 34);
    self->temp_8by8 = Matrix_new(8, 8);
    self->temp_1by8 = Matrix_new(1, 8);
    self->xPred = Matrix_new(1, self->points);
    self->exSigmaPt = Matrix_new(1, self->points);
    self->PPred = Matrix_new(8, 8);
    self->PxzPred = Matrix_new(1, self->points);
    self->K = Matrix_new(1, 8);
    self->filterFeedbackVector[0] = self->filterFeedbackVector[1] = 0;
    self->filterInputVector[0] = self->filterInputVector[1] = 0;
    
    size_t inNumberFrames = 512 / increment;
    
    self->filBuffer = Matrix_newWithValues(1, self->points, 8, 300., 2200., 3000., 3500., 20., 20., 20., 20.);
    self->filCallbackBuffer = Matrix_new(inNumberFrames, self->points);
    self->ringBuffer = RingBuffer_new(16);
    self->buffer = calloc(16, sizeof(Float64));
    self->frame = 0;
    return self;
}

void KalmanFilter_delete(KalmanFilter *self)
{
    free(self->buffer);
    RingBuffer_delete(self->ringBuffer);
    
    Matrix_delete(self->PQ);
    Matrix_delete(self->pred);
    Matrix_delete(self->xQ);
    Matrix_delete(self->wSigmaPts);
    Matrix_delete(self->wSigmaPts_xmat);
    Matrix_delete(self->wSigmaPts_zmat);
    Matrix_delete(self->xSigmaPts);
    Matrix_delete(self->Psqrtm);
    Matrix_delete(self->xPredSigmaPts);
    Matrix_delete(self->v);
    Matrix_delete(self->A);
    Matrix_delete(self->b);
    Matrix_delete(self->c);
    Matrix_delete(self->d);
    Matrix_delete(self->ahat);
    Matrix_delete(self->xhat);
    Matrix_delete(self->zPredSigmaPts);
    Matrix_delete(self->temp_8by34);
    Matrix_delete(self->exSigmaPt2);
    Matrix_delete(self->temp_1by34);
    Matrix_delete(self->ezSigmaPt2);
    Matrix_delete(self->temp_8by8);
    Matrix_delete(self->temp_1by8);
    Matrix_delete(self->xPred);
    Matrix_delete(self->exSigmaPt);
    Matrix_delete(self->PPred);
    Matrix_delete(self->PxzPred);
    Matrix_delete(self->K);
    Matrix_delete(self->filBuffer);
    Matrix_delete(self->filCallbackBuffer);
    
    free(self);
    self = NULL;
}

void KalmanFilter_process(KalmanFilter *self,
                          size_t inNumberFrames,
                          Float64 *audioIn)
{
//    KalmanFilter_filter(self, inNumberFrames, audioIn);
    if (self->frame == 0) {
        
        RingBuffer_write(self->ringBuffer, audioIn, self->ringBuffer->capacity);
        self->frame = 1;
    }
    
    for (size_t i = (self->frame == 1 ? self->ringBuffer->capacity : 0), filIncrement = (self->frame == 1 ? self->ringBuffer->capacity/self->increment : 0); i < inNumberFrames; i += self->increment, ++filIncrement) {
        self->frame = 2;
        KalmanFilter_scaledSymmetricSigmaPoints(self->xQ, self->PQ, 1., 0., 2., self->xSigmaPts, self->wSigmaPts, &self->nsp, self->Psqrtm);

        vDSP_vaddD(self->xSigmaPts->data, 1, Matrix_getRow(self->xSigmaPts, self->points), 1, self->xPredSigmaPts->data, 1, self->xPredSigmaPts->elementCount);
        
        RingBuffer_copy(self->ringBuffer, self->buffer, self->points + 1);
        RingBuffer_skip(self->ringBuffer, 1);
        RingBuffer_write(self->ringBuffer, &audioIn[i], 1);
        
        cblas_dcopy((SInt32)self->points, self->buffer, 1, self->v->data, 1);
        
        Matrix_reverse(self->v);
        Matrix_negate(self->v, self->v);
        
        for (size_t j = 0; j < self->nsp; ++j) {
            cblas_dcopy((UInt32)self->xhat->elementCount, &self->xPredSigmaPts->data[j], (UInt32)self->xPredSigmaPts->columnCount, self->xhat->data, 1);
            
            KalmanFilter_rigollRecursion(self->xhat, self->samplerate, self->ahat, self->c, self->d, self->A, self->b);
            
            Float64 dotProduct;
            vDSP_dotprD(self->ahat->data, 1, self->v->data, 1, &dotProduct, self->points);
            self->zPredSigmaPts->data[j] = dotProduct + Matrix_getRow(self->xSigmaPts, self->points + self->points)[j];
        }
        
        for (size_t j = 0; j < self->nsp - 1; ++j) {
            
            vDSP_vsubD(self->xPredSigmaPts->data, self->xPredSigmaPts->columnCount, &self->xPredSigmaPts->data[j + 1], self->xPredSigmaPts->columnCount, &self->temp_8by34->data[j], self->temp_8by34->columnCount, self->temp_8by34->rowCount);
        }
        
        Matrix_elementWiseMultiply(self->wSigmaPts_xmat, self->temp_8by34, self->temp_8by34);
        
        for (size_t j = 0; j < self->points; ++j) {
            
            vDSP_sveD(Matrix_getRow(self->temp_8by34, j), 1, &self->pred->data[j], self->temp_8by34->columnCount);
        }
        
        vDSP_vaddD(self->pred->data, 1, self->xPredSigmaPts->data, self->xPredSigmaPts->columnCount, self->pred->data, 1, self->xPred->elementCount);
        
        Float64 scalar = -self->zPredSigmaPts->data[0];
        vDSP_vsaddD(&self->zPredSigmaPts->data[1], 1, &scalar, self->temp_1by34->data, 1, self->temp_1by34->elementCount);
        vDSP_vmulD(&self->wSigmaPts->data[1], 1, self->temp_1by34->data, 1, self->temp_1by34->data, 1, self->temp_1by34->elementCount);
        
        Float64 zPred;
        vDSP_sveD(self->temp_1by34->data, 1, &zPred, self->temp_1by34->elementCount);
        
        zPred += self->zPredSigmaPts->data[0];
        
        vDSP_vsubD(self->pred->data, 1, self->xPredSigmaPts->data, self->xPredSigmaPts->columnCount, self->exSigmaPt->data, 1, self->exSigmaPt->elementCount);
        
        Float32 ezSigmaPt = self->zPredSigmaPts->data[0] - zPred;
        
        Matrix_multiply(self->exSigmaPt, true, self->exSigmaPt, false, self->temp_8by8);
        Matrix_scalarMultiply(self->temp_8by8, self->wSigmaPts->data[self->wSigmaPts->elementCount - 1], self->PPred);
        
        Matrix_scalarMultiply(self->exSigmaPt, ezSigmaPt, self->PxzPred);
        Matrix_scalarMultiply(self->PxzPred, self->wSigmaPts->data[self->wSigmaPts->elementCount - 1], self->PxzPred);
        
        
        Float64 S = self->wSigmaPts->data[self->wSigmaPts->elementCount - 1] * ezSigmaPt * ezSigmaPt;
        
        for (size_t j = 0; j < self->exSigmaPt2->columnCount; ++j) {
            
            vDSP_vsubD(self->pred->data, 1, &self->xPredSigmaPts->data[j + 1], self->xPredSigmaPts->columnCount, &self->exSigmaPt2->data[j], self->exSigmaPt2->columnCount, self->exSigmaPt2->rowCount);
        }
        
        scalar = -zPred;
        vDSP_vsaddD(&self->zPredSigmaPts->data[1], 1, &scalar, self->ezSigmaPt2->data, 1, self->ezSigmaPt2->elementCount);
        
        for (size_t j = 0; j < self->exSigmaPt2->rowCount; ++j) {
            
            vDSP_vmulD(&self->wSigmaPts->data[1], 1, Matrix_getRow(self->exSigmaPt2, j), 1, Matrix_getRow(self->temp_8by34, j), 1, self->exSigmaPt2->columnCount);
        }
        
        Matrix_multiply(self->temp_8by34, false, self->exSigmaPt2, true, self->temp_8by8);
        Matrix_add(self->temp_8by8, self->PPred, self->PPred);
        vDSP_vmulD(&self->wSigmaPts->data[1], 1, self->ezSigmaPt2->data, 1, self->temp_1by34->data, 1, self->temp_1by34->elementCount);
        
        Float64 dotProduct;
        
        vDSP_dotprD(self->temp_1by34->data, 1, self->ezSigmaPt2->data, 1, &dotProduct, self->temp_1by34->elementCount);
        S += dotProduct;
        
        Matrix_multiply(self->temp_1by34, false, self->exSigmaPt2, true, self->temp_1by8);
        vDSP_vaddD(self->PxzPred->data, 1, self->temp_1by8->data, 1, self->PxzPred->data, 1, self->PxzPred->elementCount);
        Matrix_scalarDivide(self->PxzPred, S, self->K);
        Matrix_scalarMultiply(self->K, S, self->temp_1by8);
        Matrix_multiply(self->K, true, self->temp_1by8, false, self->temp_8by8);
        
        for (size_t j = 0; j < self->points; ++j) {
            
            vDSP_vsubD(Matrix_getRow(self->temp_8by8, j), 1, Matrix_getRow(self->PPred, j), 1, Matrix_getRow(self->PQ, j), 1, self->points);
        }
        
        Float32 inovation = self->buffer[self->points] - zPred;
        Matrix_scalarMultiply(self->K, inovation, self->temp_1by8);
        vDSP_vaddD(self->pred->data, 1, self->temp_1by8->data, 1, self->filBuffer->data, 1, self->temp_1by8->elementCount);
        
        cblas_dcopy((SInt32)self->points, self->filBuffer->data, 1, self->xQ->data, 1);
        cblas_dcopy((SInt32)self->points, self->filBuffer->data, 1, Matrix_getRow(self->filCallbackBuffer, filIncrement), 1);
    }
}


void KalmanFilter_filter(KalmanFilter *self, size_t inNumberFrames, Float64 *input)
{
    Float64 bWeights[2] = {1, -0.98};
    Float64 aWeight = 1;

    Float64 feedforward, feedback;
    
    for (size_t i = 0; i < inNumberFrames; ++i) {
        
        feedback = feedforward = 0;
        self->filterInputVector[0] = input[i];
        
        for (int j = 0; j < 2; j++) {
            
            feedforward += bWeights[j] * self->filterInputVector[j];
            feedback += aWeight * self->filterFeedbackVector[j];
        }
        
        input[i] = feedforward - feedback;
        
        self->filterFeedbackVector[1] = self->filterFeedbackVector[0];
        self->filterInputVector[1] = self->filterInputVector[0];
    }
}


inline void KalmanFilter_scaledSymmetricSigmaPoints(Matrix *x,
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

inline void KalmanFilter_rigollRecursion(Matrix *xhat,
                                         Float64 Fs,
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

