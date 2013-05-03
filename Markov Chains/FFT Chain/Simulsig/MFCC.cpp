//
//  MFCC.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 28/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "MFCC.h"

MFCC::MFCC()
{
    mel_bank = NULL;
    mel_bank_c_freqs = NULL;
    mfcc_temp_vec = NULL;
    final_output = NULL;
    mel_bank_c_freqs = NULL;
}

MFCC::~MFCC()
{
    if(mel_bank){
        for(int i = 0 ; i < (*n_mel_bands); i++){
            free(mel_bank[i]);
        }
        free(mel_bank);
    }    
    if(mfcc_temp_vec)    free(mfcc_temp_vec);
    if(final_output)     free(final_output);
    if(mel_bank_c_freqs) free(mel_bank_c_freqs);
}

void MFCC::Setup(const u_short_int *s_n_bands,const u_short_int s_mn_freq,const u_short_int s_mx_freq,const u_short_int *s_fft_len,const int s_sampr)
{
    n_mel_bands = const_cast<u_short_int*>(s_n_bands);
    n_mel_b_d2_p1 = (*n_mel_bands)/2 + 1;
    fft_len = const_cast<u_short_int*>(s_fft_len);
    mn_freq = s_mn_freq;
    mx_freq = s_mx_freq;
    fft_len_d2 = (*fft_len) / 2;
    fft_len_d2_p1 = fft_len_d2 + 1;
    sampr = s_sampr;
    mel_applied_filter = (double**)malloc(sizeof(double*)*(*n_mel_bands));
    if(!mel_applied_filter){
        printf("Applied Filter Memory allocation error.\n");
    }
    for(int i = 0 ; i < (*n_mel_bands); i++){
        mel_applied_filter[i] = (double*)malloc(sizeof(double)*fft_len_d2_p1);
        memset(mel_applied_filter[i], 0, sizeof(double)*fft_len_d2_p1);
        if(!mel_applied_filter[i]){
            printf("Applied Filter allocation error.\n");
        }
    }
    mfcc_temp_vec = (double*)malloc(sizeof(double)*(*n_mel_bands));
    memset(mfcc_temp_vec, 0, sizeof(double)*(*n_mel_bands));
    final_output = (double*)malloc(sizeof(double)*n_mel_b_d2_p1);
    memset(final_output, 0, sizeof(double)*n_mel_b_d2_p1);
    dct_p = fftw_plan_r2r_1d(n_mel_b_d2_p1, mfcc_temp_vec, final_output, FFTW_REDFT10, FFTW_ESTIMATE);
    BuildFilterBank();
}

void MFCC::BuildFilterBank()
{
    printf("*-------------------Building Filter Bank----------------------*\n");

    double mel_bands_f = (*n_mel_bands) + 2;
    mel_bank_c_freqs = (double*)malloc(sizeof(double)*mel_bands_f);
    memset(mel_bank_c_freqs, 0, sizeof(double)*mel_bands_f);
    double mx_mel = FtoMel(mx_freq);
    double mn_mel = FtoMel(mn_freq);
    double mel_freq_dif = mx_mel - mn_mel;
    double band_w = mel_freq_dif / (mel_bands_f-1);
    debug_print(("MinFreq   : %d\n",mn_freq));
    debug_print(("MaxFreq   : %d\n",mx_freq));
    debug_print(("MaxMel    : %f\n",mx_mel));
    debug_print(("MinMel    : %f\n",mn_mel));
    debug_print(("Freq Diff : %f\n",mel_freq_dif));
    debug_print(("Band W    : %f\n",band_w));
    mel_bank_c_freqs[0] = mn_mel;
    
    debug_print(("\n*-------Mel filters fc-------------*\n"));
    for(int i = 1 ; i < mel_bands_f; i++){
        mel_bank_c_freqs[i] = mel_bank_c_freqs[i-1] + band_w;
        debug_print(("%f ",mel_bank_c_freqs[i]));
        if((i%5==4)) debug_print(("\n"));
    }
    debug_print(("\n*-------Mel fc lin scale-----------*\n"));
    for(int i = 0 ; i < mel_bands_f; i++){
        mel_bank_c_freqs[i] = MtoFreq(mel_bank_c_freqs[i]);
        debug_print(("%f ",mel_bank_c_freqs[i]));
        if((i%5==4)) debug_print(("\n"));
    }
    debug_print(("\n*-------Mel fc scaled fft length---*\n"));
    // now need to scale this appropriately so that the length
    // of each filterbank matrix is the fft size
    for(int i = 0 ; i < mel_bands_f; i++){
        mel_bank_c_freqs[i] *= (*fft_len);
        mel_bank_c_freqs[i] /= sampr;
        debug_print(("%f ",mel_bank_c_freqs[i]));
        if((i%5==4)) debug_print(("\n"));
    }
    printf("\n*----------------------------------*\n");
    printf("Building Filterbank Matrix...\n");
    mel_bank = (double**)malloc(sizeof(double*)*(*n_mel_bands));
    if(!mel_bank){
        printf("Filterbank Memory allocation error.\n");
        return;
    }
    for(int i = 0 ; i < (*n_mel_bands); i++){
        mel_bank[i] = (double*)malloc(sizeof(double)*fft_len_d2);
        if(!mel_bank[i]){
            printf("Filterbank Memory allocation error.\n");
            return;
        }
    }
    // apply triangle filter function
    for(int i = 0 ; i < (*n_mel_bands); i++){
        BuildFilter(mel_bank[i], mel_bank_c_freqs, i+1, fft_len_d2);
    }
    debug_print(("Filterbank build complete\n"));
    printf("*---------------------------------------------------------------*\n\n");

}


void MFCC::BuildFilter(double* result,double* filter_bank,const int filter_index,u_short_int size)
{
    // Eq ref'd in dissertation.
    double k0 = filter_bank[filter_index-1];
    double k1 = filter_bank[filter_index];
    double k2 = filter_bank[filter_index+1];

    for(int i = 0 ; i < size; i++){
        if(i < k0)            { result[i] = 0;                            continue; }
        if(i >= k0 && i <= k1){ result[i] = (((double)i - k0)/(k1 - k0)); continue; }
        if(i >= k1 && i <= k2){ result[i] = ((k2-i)/(k2-k1));             continue; }
        if(i >= k2)           { result[i] = 0;                            continue; }
    }
}

double* MFCC::CalcMFCC(double *buffer)
{
    // amplitude spectrum log
    // for each mel freq bank. multiply power signal by it.
    for(int i = 0 ; i < (*n_mel_bands); i++){
        vDSP_vmulD(buffer, 1, mel_bank[i], 1, mel_applied_filter[i], 1, fft_len_d2_p1);       
    }
    memset(mfcc_temp_vec, 0, sizeof(double)*(*n_mel_bands));
    // add up all of these. could do this within the loop
    for(int i = 0 ; i < (*n_mel_bands); i++){
        for(int j = 0 ; j < fft_len_d2_p1; j++){
            mfcc_temp_vec[i] += mel_applied_filter[i][j];
            if(isnan(mel_applied_filter[i][j])) printf("NANO %d : %d \n",i,j);
        }
    }
    double t = 0.0;
    for(int i = 0; i < (*n_mel_bands); i++){
        if(mfcc_temp_vec[i] != 0.0 ){
            t = mfcc_temp_vec[i];
           // printf("%f \n",mfcc_temp_vec[i]);
            mfcc_temp_vec[i] = log(mfcc_temp_vec[i]);
        }
    }
    fftw_execute(dct_p);
    
    final_output[0] = 0.0;
    return final_output;
}

double MFCC::FtoMel(double freq)
{
    return (2595 * log10(1.0 + ((double)freq/700.0)));
}

double* MFCC::GetMF()
{
    return mfcc_temp_vec;
}

double MFCC::MtoFreq(double mel)
{
    return 700.0 * (pow(10,mel/2595.0)-1.0);
}

std::vector<double>* MFCC::LinFreqGraph()
{
    for(int i = 0 ; i < 22050; i++){
        lin_freq_graph.push_back((double)i);
    }
    return &lin_freq_graph;
}

std::vector<double>* MFCC::MelFreqGraph()
{
    for(int i = 0 ; i < 22050; i++){
        double mel = 2595 * log10(1.0 + ((double)i/700.0));
        mel_freq_graph.push_back(mel);
    }
    return &mel_freq_graph;
}

double** MFCC::GetFilterBank()
{
    return mel_bank;
}

int MFCC::GetNumBands()
{
    return (*n_mel_bands);
}

u_short_int MFCC::GetVecLength()
{
    return fft_len_d2;
}

u_short_int* MFCC::GetVelLengthAdr()
{
    return &fft_len_d2;
}


