//
//  Concat.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 12/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "Player.h"

Player::Player()
{
    
}

Player::~Player()
{
    
}

void Player::Setup(int fft_size,int hop_size)
{
    database.ReadFile(0);
    database2.ReadFile(0);
    hop_overlap =  hop_size;
    vector_index = 0;
    vector_index_offset = 0;
    single_frame = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * database.GetFFTLength()/2);
    all_frames = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* database.GetNumFrames());
    fft_interleaved = (double*)malloc(sizeof(double)*database.GetFFTLength());
    
    real_signal = (double*)malloc(sizeof(double)*(database.GetFFTLength()-2));
    ifft = fftw_plan_dft_c2r_1d(fft_size, single_frame, real_signal, FFTW_BACKWARD);
    
    //ham.MakeWindow(fft_size);

}

void Player::Play()
{
    double* mfcc_one = (double*)malloc(sizeof(double)*database.GetNumMFCCCoefs());
    int s = 0;
    for(int i = 0; i < database.GetNumFrames(); i++){
        mfcc_one = database.GetMFCCMatrix()[i];
        double nearest = EuclidDistance(mfcc_one, database2.GetMFCCMatrix()[0]);
        int nearest_index = 0;
        for(int j = 0; j < database2.GetNumFrames();j++){
            double di = EuclidDistance(mfcc_one, database2.GetMFCCMatrix()[j]);
            if(di < nearest){
                nearest = di;
                nearest_index = j;
            }
        }
        fft_interleaved = database2.GetFFTMatrix()[nearest_index];
        IntlvdToFFTW(fft_interleaved);
        InverseFFTW();
        //
        
        RealHopOvrAddPushback();
        //RealPushback();
      //  FFTFramePushback(fft_interleaved);
        s = i;
    }
    debug_print(("I : %d\n",s));
    WriteSignal();
  //  VecToFFTW();
}

void Player::FFTFramePushback(double *fft_values)
{
    for(int i = 0 ; i < database.GetFFTLength();i+=2){
        frames.push_back(std::make_pair(fft_values[i], fft_values[i+1]));
    }
}

void Player::WriteSignal()
{
    ///Users/TomRushmore/Final Project/Signals
    ///Users/thomasrushmore/Documents/Goldsmiths/Third Year/Final Project/Signals/sig.txt"
    std::string dir = MACHINEPATH;
    dir.append("/Signals/sig.txt");
    std::ofstream signal(dir.c_str());
    signal << std::fixed << std::setprecision(10);

    printf("Writing signal summary...\n");
    debug_print(("Real Frames Size = %ld \n",real_frames.size()));
    for(int i = 0 ; i < real_frames.size(); i++){
        signal << real_frames[i] << " ";
    }
    signal.close();
    printf("Signal write complete\n");
}

void Player::RealHopOvrAddPushback()
{
    int ho =4096 - hop_overlap;
    if(vector_index_offset == 0){
        for(int i = 0 ; i < database.GetFFTLength()-2;i++){
            real_frames.push_back(real_signal[i]);
            vector_index++;
        }
        vector_index_offset = 1;
    } else {
        for(int i = 0 ; i < database.GetFFTLength()-2;i++){
            if(i < ho){
                real_frames[vector_index - ho + i] += real_signal[i];
            } else {
                real_frames.push_back(real_signal[i]);
                vector_index++;
            }
        }
    }
    
}

void Player::RealPushback()
{
    for(int i = 0 ; i < database.GetFFTLength()-2;i++){
        real_frames.push_back(real_signal[i]);
    }

}

double* Player::InverseFFTW()
{
    fftw_execute(ifft);
    double max = 0.0;
    for(int i = 0 ; i < database.GetFFTLength()-2;i++){
        if(real_signal[i] > max){
            max = real_signal[i];
        }
    }
    return real_signal;
}
fftw_complex* Player::IntlvdToFFTW(double *interleaved_fft)
{
    for(int i = 0 ; i < database.GetFFTLength()/2; i++){
        single_frame[i][0] = interleaved_fft[i*2];
        single_frame[i][1] = interleaved_fft[i*2+1];

    }
    return single_frame;
}

fftw_complex* Player::VecToFFTW()
{
    for(int i = 0 ; i < frames.size();i++){
        all_frames[i][0] = frames[i].first;
        all_frames[i][1] = frames[i].second;
    }
    return all_frames;
}

double Player::EuclidDistance(double *a, double *b)
{
    double running_sum = 0.0;
    for(int i = 0 ; i < database.GetNumMFCCCoefs(); i++){
        running_sum += pow((a[i] - b[i]),2);
    }
    running_sum = sqrt(running_sum);
    return running_sum;
}