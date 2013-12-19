//
//  FeatureExtraction.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 15/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "FeatureExtraction.h"

FeatureExtraction::FeatureExtraction()
{
    total_num_frames = 0;
    fft_len          = 4096;
    hop_ovrlp        = 0;
    n_mfcc_coefs     = 0;
    win_chunk = NULL;
    n_win_chunk = NULL;
}

FeatureExtraction::~FeatureExtraction()
{
    
}

bool FeatureExtraction::Setup(const int s_fft_size,const  int s_hop_ovrlp,const int n_mel_bands,const int sample_rate,const double min_freq,const double max_freq,const int samp_chunk_size,const int mfcc_grouping)
{
    if(!IsPowerTwo(s_fft_size) && s_fft_size < type_max_lim(u_short_int)){
        printf("SETUP ERROR. Analysis size must be a power of 2, currently %d\n",s_fft_size);
        return false;
    } else {
        fft_len = s_fft_size;
    }
    s_hop_ovrlp <= 0 ? hop_ovrlp = 1 : hop_ovrlp = s_hop_ovrlp;
 
    if(n_mel_bands > type_max_lim(u_short_int) || n_mel_bands < 1){
        printf("SETUP ERROR. Number of mel bands must be in range 0 < n < %d . currently : %d",type_max_lim(u_short_int),n_mel_bands);
        return false;
    } else {
        n_mfcc_coefs = n_mel_bands;
        fin_n_mfcc_coefs = n_mfcc_coefs/2;
    }
    sample_chunk_size = samp_chunk_size;
    sample_chunk_offset = (fft_len * mfcc_grouping) - sample_chunk_size;
    file_reader.SetChunkSize(&fft_len);
    file_reader.SetHopSize(&hop_ovrlp);
    double * ham = hamming_window.MakeWindow(&fft_len);
    hamming_directory = file_writer.WriteVector(ham, &fft_len, "Hamming");
    ham = NULL;
    fft.Setup(&fft_len);
    database.SetFixedFramesLength(&fft_len);
    mfcc.Setup(&n_mfcc_coefs, min_freq , max_freq,&fft_len,sample_rate);
    //file_writer.WriteMelBanks(const_cast<const double**>(mfcc.GetFilterBank()), 20000, NMELBANDS/2);
    file_reader.SetChunkReadOffset(sample_chunk_offset);
    ConsoleDebug();
    return true;
}

void FeatureExtraction::ProcessSingleFile(const char *file_path,const PROCESSTYPE type)
{
    // Calculate the number of frames to process
    // Setup of Hop Size for audio reader
    int num_frames = file_reader.Load(file_path);
    int orig = num_frames;
    printf("*------------Processing file-------------*\n");
    debug_print(("File size : %d samples\n",num_frames));
    num_frames = (int)floor((num_frames/(double)fft_len));
    debug_print(("Num frames of %d fft length without hop    : %d frames\n",fft_len,num_frames));
    num_frames = orig / hop_ovrlp;
    debug_print(("Num frames of %d fft length with %d ms hop : %d frames\n",fft_len,hop_ovrlp,num_frames));
    debug_print(("Setting hop ovrlp in FE : %d\n",hop_ovrlp));
    file_reader.SetHopSize(&hop_ovrlp);
    
    // WRITE Enum leftover from previous Feature Extraction Class
    // where only the Magnitudes etc could be written
    if(WRITE == PAIR){
        // Create the fftw types and set the database and filereader.
        fft_write = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft.GetFFTVectorLength());
        fft_write_no_w = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft.GetFFTVectorLength());
        database.SetFixFFTMFCCLength(fft.GetFFTVectorLengthAdr(), &fin_n_mfcc_coefs);
        total_num_frames += num_frames;
        database.AddNumFramesPerTrack(num_frames);
        database.SetBufferSizes(file_reader.GetBufferSize());
        int frame_length = floor(num_frames/MFCCGROUPING);
        frame_length *= (MFCCGROUPING*fft_len);
        frame_length /= fft_len;
        // Windowed segment FFT used for MFCC calc
        // Unwindowed segment for later IFT in reconstruction phase
        for(int i = 0 ; i < frame_length  ; i++){
            win_chunk = file_reader.GetChunk();
            n_win_chunk = win_chunk;
            fft.ProcessFFTR2CBatchMagFFT(n_win_chunk, fft_write_no_w);
            win_chunk = hamming_window.Process(win_chunk);
            win_chunk = fft.ProcessFFTR2CBatchMagFFT(win_chunk, fft_write);
            win_chunk = mfcc.CalcMFCC(win_chunk);
            database.AddFFTMFCCSplit(fft_write_no_w, win_chunk);
        }
        fftw_free(fft_write);
        fftw_free(fft_write_no_w);
    }
    printf("*----------------------------------------*\n\n");
}

void FeatureExtraction::WriteDatabase(long num_tracks,WRITER type)
{
     file_writer.DBWriteFFTMFCCSplit(database.GetFFTVector(), database.GetMFCCVector(),database.GetTrackNumFramesVector(), GetFinMBLen(),fft_len, total_num_frames,0,num_tracks,type,database.GetBufferSizes());
}

void FeatureExtraction::ConsoleDebug()
{
    printf("*---------Feature Extraction Setup-----------*\n");
    debug_print(("Feature Extraction Setup\n"));
    debug_print(("Frame Size            : %d samples\n",fft_len));
    debug_print(("Chunk Size            : %d samples\n",file_reader.GetChunkSize()));
    debug_print(("FFT Vector Length     : %d samples\n",fft.GetFFTVectorLength()));
    debug_print(("DB Write Frame Length : %d samples\n",database.GetFramesLength()));
    debug_print(("Hamming win directory : %s \n",hamming_directory.c_str()));
    printf("*--------------------------------------------*\n\n");
}

void FeatureExtraction::ProcessFile(const char *file_path)
{
    int num_frames = file_reader.Load(file_path);
    num_frames = (int)ceil((num_frames/(double)fft_len));
    debug_print(("Num frames without hop : %d",num_frames));
    num_frames *= (int)(fft_len / hop_ovrlp);
    debug_print(("Num frames with hop : %d",num_frames));
    
    for(int i = 0; i < num_frames; i++){
        double * ptr = file_reader.GetChunk();
        file_writer.WriteVector(ptr, fft_len);
        ptr = hamming_window.Process(ptr);
        file_writer.WriteVector(ptr, fft_len);
        double * vec = fft.ProcessFFTR2C(ptr);
        file_writer.WriteVector(vec, fft_len);
    }
}

void FeatureExtraction::SetWriteFilePath(const char *write_file_path)
{
    directory = const_cast<char *>(write_file_path);
    file_writer.SetFilePath(directory);
}

//http://stackoverflow.com/questions/600293/how-to-check-if-a-number-is-a-power-of-2

bool FeatureExtraction::IsPowerTwo(int number)
{
    return (number & (number - 1)) == 0;
}

// used for dissertation value print outs
void FeatureExtraction::ProcessDatabase(const char *file_path,int file_num)
{
    printf("*------------Processing Database-------------*\n");
    int num_frames = file_reader.Load(file_path);
    int orig = num_frames;
    debug_print(("Buffer size : %d samples\n",num_frames));
    num_frames = (int)ceil((num_frames/(double)fft_len));
    debug_print(("Num %d size frames without hop : %d\n",fft_len,num_frames));
    num_frames = orig / hop_ovrlp;
    debug_print(("Num %d size frames with %d ms hop : %d\n",fft_len,hop_ovrlp,num_frames));
    int ovr = fft_len -1;
    file_reader.SetHopSize(&hop_ovrlp);
    
    double *b_ptr = NULL;
    if(WRITE == MAGN){
        for(int i = 0; i < num_frames; i++){
            b_ptr = file_reader.GetChunk();
            b_ptr = hamming_window.Process(b_ptr);
            b_ptr = fft.ProcessFFTR2CBatchMag(b_ptr);
            database.AddVectorMag(b_ptr);
        }
    } else
        if(WRITE == MFCC){
            //  database.SetFixedFramesLength(13);
            for(int i = 0; i < num_frames - 1; i++){
                b_ptr = file_reader.GetChunk();
                if(i >= num_frames-ovr){
                    //     ptr = ZeroPad(ptr);
                }
                hamming_window.Process(b_ptr);
                b_ptr = fft.ProcessFFTR2CBatchMag(b_ptr);
                b_ptr = mfcc.CalcMFCC(b_ptr);
                database.AddVectorMag(b_ptr);
            }
        } else
            if(WRITE == PAIR){
                fft_write = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft.GetFFTVectorLength());
                // needs to be done before hand ( all files will have the same)
                database.SetFixFFTMFCCLength(fft.GetFFTVectorLengthAdr(),&n_mfcc_coefs);
                database.AddNumFramesPerTrack(num_frames - 1);
                for(int i = 0 ; i < num_frames - 1; i++){
                    b_ptr = file_reader.GetChunk();
                    b_ptr = hamming_window.Process(b_ptr);
                    b_ptr = fft.ProcessFFTR2CBatchMagFFT(b_ptr, fft_write);
                    b_ptr = mfcc.CalcMFCC(b_ptr);
                    //database.AddFFTMFCCPair(fft_write, b_ptr);
                    database.AddFFTMFCCSplit(fft_write, b_ptr);
                }
                file_writer.DBWriteFFTMFCCSplit(database.GetFFTVector(), database.GetMFCCVector(),database.GetTrackNumFramesVector(),n_mfcc_coefs,fft_len, num_frames -1 ,file_num,0,WRITEFILE,database.GetBufferSizes());
                free(fft_write);
            }
    printf("*--------------------------------------------*\n\n");
    
}

u_short_int FeatureExtraction::GetFinMBLen()
{
    return fin_n_mfcc_coefs;
}