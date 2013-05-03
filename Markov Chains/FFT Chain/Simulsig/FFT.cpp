//
//  FFT.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 15/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "FFT.h"

FFT::FFT()
{
    zero_ref = 1.0;
    r2c_in = NULL;
    write_buffer = NULL;
    mag_buf = NULL;
    vdsp_cpx_buf = NULL;
}

FFT::~FFT()
{
    if(r2c_plan)     fftw_destroy_plan(r2c_plan);
    if(r2c_in)       free(r2c_in);
    if(write_buffer) free(write_buffer);
    if(mag_buf)      free(mag_buf);
    if(vdsp_cpx_buf) free(vdsp_cpx_buf);
}

void FFT::Setup(const u_short_int *s_fft_length)
{
    fft_length = const_cast<u_short_int*>(s_fft_length);
    fft_length_d2_p1 = ((*fft_length)/2) + 1;
 
    write_buffer = (double*)malloc(sizeof(double) * (*fft_length));
    mag_buf  = (double*)malloc(sizeof(double)* fft_length_d2_p1);
    r2c_in   = (double*)malloc(sizeof(double*)*(*fft_length));
    r2c_out  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_length_d2_p1);
    r2c_plan = fftw_plan_dft_r2c_1d((*fft_length), r2c_in, r2c_out, FFTW_ESTIMATE);
    
    vdsp_cpx_buf = (DSPDoubleSplitComplex*)malloc(sizeof(DSPDoubleSplitComplex));
    vdsp_cpx_buf->imagp = (double*)malloc(sizeof(double)*fft_length_d2_p1);
    vdsp_cpx_buf->realp = (double*)malloc(sizeof(double)*fft_length_d2_p1);

}

double* FFT::ProcessFFTR2CBatchMag(double *buffer)
{
    memcpy(r2c_in, buffer, (*fft_length)*sizeof(double*));
    fftw_execute(r2c_plan);
    switch(FTWVDSPCONVTYPE){
        case REINTERPRET: vdsp_cpx_buf = reinterpret_cast<DSPDoubleSplitComplex*>(r2c_out); break;
        case BRUTE      : FFTWConvertVDSP(r2c_out); break;
        default: break;
    }
    mag_buf = SquaredMagnitude(vdsp_cpx_buf);
    mag_buf = SquareRoot(mag_buf);
    return mag_buf;
}

double* FFT::ProcessFFTR2CBatchMagFFT(double *buffer, fftw_complex *fft_write_buffer)
{
    memcpy(r2c_in, buffer, (*fft_length)*sizeof(double));
    fftw_execute(r2c_plan);
    for(int i = 0 ; i < fft_length_d2_p1; i++){
        fft_write_buffer[i][0] = r2c_out[i][0];
        fft_write_buffer[i][1] = r2c_out[i][1];
    }
    // Two conversion types as BRUTE previously had problems.
    switch(FTWVDSPCONVTYPE){
        case REINTERPRET: vdsp_cpx_buf = reinterpret_cast<DSPDoubleSplitComplex*>(r2c_out); break;
        case BRUTE      : FFTWConvertVDSP(r2c_out); break;
        default: break;
    }
    mag_buf = SquaredMagnitude(vdsp_cpx_buf);
    mag_buf = SquareRoot(mag_buf);
    return mag_buf;
}

DSPDoubleSplitComplex* FFT::FFTWConvertVDSP(fftw_complex *fftw)
{
    for(int i = 0; i < fft_length_d2_p1; i++){
        vdsp_cpx_buf->realp[i] = fftw[i][0];
        vdsp_cpx_buf->imagp[i] = fftw[i][1];
    }
    return vdsp_cpx_buf;
}

double* FFT::ProcessFFTR2C(double *buffer)
{
    memcpy(r2c_in, buffer, (*fft_length)*sizeof(double));
    fftw_execute(r2c_plan);
    for(int i = 0 ; i < (*fft_length); i++){
        write_buffer[i] = r2c_out[i][0];
    }
    vDSP_vabsD(write_buffer, 1, write_buffer, 1, fft_length_d2_p1);
    return write_buffer;
}

inline double* FFT::SquareRoot(double *buffer)
{
    for(int i = 0 ; i < fft_length_d2_p1; i++){
        buffer[i] = sqrt(buffer[i]);
    }
    return buffer;
}

inline double* FFT::SquaredMagnitude(DSPDoubleSplitComplex * s_buffer)
{
    vDSP_zvmagsD(s_buffer, 1, mag_buf, 1, fft_length_d2_p1);
    return mag_buf;
}

inline double* FFT::ConvToDB(double *buffer)
{
    for(int i = 0 ; i < fft_length_d2_p1; i++){
       buffer[i] = log10(buffer[i]);
    }
    vDSP_vdbconD(buffer, 1, &zero_ref, buffer, 1, fft_length_d2_p1, 1);
    return buffer;
}

u_short_int FFT::GetFFTVectorLength()
{
    return fft_length_d2_p1;
}

u_short_int* FFT::GetFFTVectorLengthAdr()
{
    return &fft_length_d2_p1;
}