//
//  FFT.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 15/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__FFT__
#define __FFT_Chain__FFT__

#include <Accelerate/Accelerate.h>
#include <math.h>
#include <iostream>

#include "fftw3.h"
#include "PrintfDebug.h"
#include "GlobalTypedefs.h"

#define FTWVDSPCONVTYPE BRUTE

class FFT{
public:
    FFT();
    ~FFT();
    void    Setup(const u_short_int *s_fft_length);
    double* ProcessFFTR2C(double * buffer);
    double* ProcessFFTR2CBatchMag(double * buffer);
    double* ProcessFFTR2CBatchMagFFT(double * buffer, fftw_complex* fft_write_buffer);
    u_short_int GetFFTVectorLength();
    u_short_int* GetFFTVectorLengthAdr();

    
private:
    inline double* SquaredMagnitude(DSPDoubleSplitComplex *buffer);
    inline double* SquareRoot(double* buffer);
    inline double* ConvToDB(double * buffer);

    DSPDoubleSplitComplex* FFTWConvertVDSP(fftw_complex* fftw);
    double zero_ref;
    double *write_buffer;
    double *mag_buf;
    u_short_int *fft_length;
    u_short_int fft_length_d2_p1;
    fftw_complex *r2c_out;
    fftw_plan r2c_plan;
    double* r2c_in;
    DSPDoubleSplitComplex *vdsp_cpx_buf;
};

enum CONVTYPE{
    REINTERPRET = 0,
    BRUTE
};

#endif /* defined(__FFT_Chain__FFT__) */
