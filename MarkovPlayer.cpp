//
//  MarkovPlayer.cpp
//  Simulsig
//
//  Created by Thomas Rushmore on 02/04/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "MarkovPlayer.h"

MarkovPlayer::MarkovPlayer()
{
    vector_index = 0;
    all_time_doms.resize(0);
}

MarkovPlayer::~MarkovPlayer()
{
    
}

void MarkovPlayer::Setup(MarkovChain* s_mkv_ch,double** s_fft_matrix,int s_n_sequences,int s_max_track_length,int s_n_coef,int s_n_active_clusters,int s_fft_length,int s_grouping,int s_hop_overlap,int s_ramp)
{
    mkv_ch = s_mkv_ch;
    n_sequences = s_n_sequences;
    max_track_length = s_max_track_length;
    fft_matrix = s_fft_matrix;
    n_coef = s_n_coef;
    n_active_clusters = s_n_active_clusters;
    fft_length = s_fft_length ;
    fft_length_m2 = fft_length-2;
    fft_length_2_m1 = (fft_length_m2/2)+1;
    grouping = s_grouping;
    hop_overlap = s_hop_overlap;
    rem_frames = int((22050)) % fft_length_m2;
    ramp = int((s_ramp-2)/4);
    forward_ramp = (double*)malloc(sizeof(double)*ramp);
}

void MarkovPlayer::GenerateSampleBuffers()
{
    time_dom_signals = (double**)malloc(sizeof(double)*n_active_clusters);
    for(int i = 0 ; i < n_active_clusters; i++){        
        time_dom_signals[i] = (double*)malloc(sizeof(double)*fft_length_m2*grouping);
    }
    real_signal = (double*)malloc(sizeof(double)*fft_length_m2);
    memset(real_signal, 0, sizeof(double)*fft_length_m2);
    single_frame = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_length_2_m1);
    ifft = fftw_plan_dft_c2r_1d(fft_length_m2, single_frame, real_signal, FFTW_ESTIMATE);
    for(int i = 0 ; i < n_active_clusters; i++){
        int time_offset = 0;
        int other_offset = 0;
        for(int j = 0 ; j < grouping; j++){
            single_frame_interleaved = &fft_matrix[i][time_offset];
            IntlvdToFFTW(single_frame_interleaved);
            InverseFFTW();
            for(int p = 0; p < fft_length_m2; p++){
                time_dom_signals[i][p+other_offset] = real_signal[p]/fft_length_m2;
            }
            time_offset = (j+1)*fft_length;
            other_offset =(j+1)*fft_length_m2;
        }
        RealHopOvrAddPushback(i);
    }
    WriteSignal();
}

void MarkovPlayer::RealHopOvrAddPushback(int sample_num)
{
    int idx = 0;
    cur_track.clear();
    vector_index = 0;
    int temp_hop = 0;
    
    for(int i = 0 ; i < fft_length_m2;i++){
        cur_track.push_back(time_dom_signals[sample_num][idx++]);
        vector_index++;
    }
    for(int j = 0; j <  grouping - 2; j++){
        for(int i = 0 ; i < fft_length_m2;i++){
            if(i < temp_hop){
                cur_track[vector_index - temp_hop + i - 1] = time_dom_signals[sample_num][idx++];
            } else {
                cur_track.push_back(time_dom_signals[sample_num][idx++]);
                vector_index++;
            }
            
        }
      }
    for(int i = 0 ; i < rem_frames;i++){
        cur_track.push_back(time_dom_signals[sample_num][idx++]);
    }
    final_signals.push_back(cur_track);
}



fftw_complex* MarkovPlayer::IntlvdToFFTW(double *interleaved_fft)
{
    for(int i = 0 ; i < fft_length_2_m1; i++){
        single_frame[i][0] = interleaved_fft[i*2];
        single_frame[i][1] = interleaved_fft[i*2+1];
    }
    return single_frame;
}

double* MarkovPlayer::InverseFFTW()
{
    fftw_execute(ifft);
    double max = 0.0;
    for(int i = 0 ; i < fft_length;i++){
        if(real_signal[i] > max){
            max = real_signal[i];
        }
    }
    return real_signal;
}

void MarkovPlayer::GenerateChainSequences(const MKVTYPE type)
{
    single_chain_play.resize(0);
    int r_idx = 0;
    for(int j = 0 ; j < n_sequences; j++){
        for(int i = 0 ; i < max_track_length; i++){
            type == HOMOG ? r_idx = mkv_ch->letsplay() : r_idx = mkv_ch->stdletsplay();
            single_chain_play.push_back(r_idx);
        }
        mkv_ch->letsplayreset();
        chain_plays.push_back(single_chain_play);
        single_chain_play.clear();
    }
    for(int i = 0 ; i < chain_plays.size(); i++){
        printf("Sequence %d : \n",i);
        for(int j = 0 ; j < chain_plays[i].size(); j++){
            printf("%d ",chain_plays[i][j]);
        }
        printf("\n\n");
    }
}

void MarkovPlayer::GenerateSignals()
{
    // get the chain sequence.
    // create vector for the new sample to be created
    // for each sequence value get the corresponding signal
    // append the signal to the end
    // write to file
    abs_signal.resize(0);
    BuildRamp();
    long signal_length = final_signals[0].size();
    int write_offset = 0;
    int ramp_offset = 0;
    for(int i = 0; i < n_sequences; i++){
        std::vector<double> newsig;
        write_offset = 0;
        ramp_offset = 0;
        for(int j = 0 ; j < max_track_length; j++){
            int seq_idx = chain_plays[i][j];
            if(!j){
                for(int q = 0; q < signal_length; q++){
                    if(q > signal_length - ramp){
                        newsig.push_back(final_signals[seq_idx][q] * forward_ramp[ramp-1-ramp_offset]);
                        ramp_offset++;
                    } else {
                        newsig.push_back(final_signals[seq_idx][q]);
                        
                    }
                    write_offset++;
                 }
                ramp_offset = 0;
            } else {
                for(int q = 0; q < signal_length; q++){
                    if(q < ramp){
                        newsig[write_offset-ramp+q] += final_signals[seq_idx][q] * forward_ramp[q];
                    } else {
                        if(q > signal_length - ramp){
                            newsig.push_back(final_signals[seq_idx][q] * forward_ramp[ramp-1-ramp_offset]);
                            ramp_offset++;
                        } else {
                            newsig.push_back(final_signals[seq_idx][q]);
                        }
                        write_offset++;
                    }
                }
                ramp_offset = 0;
            }
        }
        abs_signal.push_back(newsig);
        newsig.clear();
    }
}

void MarkovPlayer::WriteFinalSignals()
{
    std::string dir = MACHINEPATH;
    const int format=SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    const int channels=1;
    const int sampleRate=44100;
    dir.append("Final Tracks/");
    char buf[10];
    for(int i = 0 ; i < n_sequences; i++){
        printf("Writing signals %d...\n",i);
        sprintf(buf, "%d", i);
        std::string audio_path = dir;
        audio_path.append(buf);
        audio_path.append(".wav");
        SndfileHandle audio_file(audio_path, SFM_WRITE, format, channels, sampleRate);
        audio_file.write(&abs_signal[i][0], abs_signal[0].size());
        printf("Signal write complete\n");
    }
       
}

void MarkovPlayer::BuildRamp()
{
    double ramp_offval = 1.0 / ramp;
    for(int i = 0 ; i < ramp; i++){
        forward_ramp[i] = i * ramp_offval;
    }
}

void MarkovPlayer::WriteSignal()
{
    std::string dir = MACHINEPATH;
    dir.append("/Signals/signal_");
    char buf[10];
    for(int i = 0 ; i < n_active_clusters; i++){
        std::string path = dir;
        sprintf(buf, "%d", i);
        path.append(buf);
        path.append(".txt");
        std::ofstream signal(path.c_str());
        signal << std::fixed << std::setprecision(10);
        debug_print(("Writing signal summary %d...\n",i));
        for(int j = 0 ; j < cur_track.size(); j++){
            signal << final_signals[i][j] << " ";
        }
        signal.close();
        debug_print(("Signal write complete\n"));
    }
}

/*
 void MarkovPlayer::WriteFinalSignals()
 {
 std::string dir = MACHINEPATH;
 const int format=SF_FORMAT_WAV | SF_FORMAT_PCM_16;
 const int channels=1;
 const int sampleRate=44100;
 
 dir.append("/Signals/Finals/signal_");
 char buf[10];
 double* write_buf = (double*)malloc(sizeof(double)*abs_signal[0].size());
 for(int i = 0 ; i < n_sequences; i++){
 std::string path = dir;
 sprintf(buf, "%d", i);
 path.append(buf);
 path.append(".txt");
 std::ofstream signal(path.c_str());
 signal << std::fixed << std::setprecision(10);
 printf("Writing signal summary %d...\n",i);
 for(int j = 0 ; j < abs_signal[0].size(); j++){
 signal << abs_signal[i][j] << " ";
 write_buf[j] = abs_signal[i][j];
 }
 std::string audio_path = dir;
 audio_path.append(buf);
 audio_path.append(".wav");
 SndfileHandle outfile(audio_path, SFM_WRITE, format, channels, sampleRate);
 outfile.write(&write_buf[i], abs_signal[0].size());
 signal.close();
 printf("Signal write complete\n");
 }
 free(write_buf);
 
 }
*/