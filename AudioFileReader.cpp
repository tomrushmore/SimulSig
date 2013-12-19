//
//  AudioFileReader.cpp
//  test
//
//  Created by Thomas Rushmore on 10/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "AudioFileReader.h"

AudioFileReader::AudioFileReader()
{
    chunk_size = NULL;
    chunk_iterator = 0;
    buffer = NULL;
    buffer_size = 0;
    chunk_buffer = NULL;
}

AudioFileReader::~AudioFileReader()
{
    if(chunk_buffer){
       // free(chunk_buffer);
    }
    if(buffer){
        free(buffer);
    }
}

int AudioFileReader::Load(const char *s_file_path)
{

    printf("*--------------Loading Audio File-----------------*\n");
    if(buffer){
        printf("REALLOCATING BUFFER\n");
        free(buffer);
    }
    file_path = const_cast<char*>(s_file_path);
    SNDFILE *sf;
    SF_INFO info;
    sf = sf_open(file_path, SFM_READ, &info);
    printf("Audio File directory : %s \n",file_path);
    if (sf == NULL){
        printf("Failed to open the file.\n");
        exit(-1);
    }
    printf("File Loaded\n");
    // following code taken from libsnd file documentation
    int f,sr,c,num_items,num;
    f = static_cast<int>(info.frames);
    sr = info.samplerate;
    c = info.channels;
    num_items = f*c;
    buffer_size = num_items;
    chunk_ptr = NULL;
    
    buffer = (double *) malloc(buffer_size*sizeof(double));
    num = static_cast<int>(sf_read_double(sf,buffer,num_items));
    sf_close(sf);
    chunk_iterator = 0;
    del = 0;
    debug_print(("Sample Rate  :%d\n",sr));
    debug_print(("Num Channels :%d\n",c));
    debug_print(("Num Samples  :%ld\n",buffer_size));
    printf("*-----------Audio File Load Complete--------------*\n\n");
    
    return num_items;
}

int AudioFileReader::WriteFile(const char *s_file_path)
{
    const int format=SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    const int channels=1;
    const int sampleRate=48000;
    const char* outfilename="fooat.wav";
    SndfileHandle outfile(outfilename, SFM_WRITE, format, channels, sampleRate);
    const int size = sampleRate*3;
    outfile.write(buffer, size);
    
    return 0;
}

double* AudioFileReader::GetBuffer()
{
    return buffer;
}

double AudioFileReader::GetSample(int index)
{
    return buffer[index];
}

long AudioFileReader::GetBufferSize()
{
    return buffer_size;
}

void AudioFileReader::SetChunkSize(const u_short_int *set_chunk_size)
{
    chunk_size = const_cast<u_short_int*>(set_chunk_size);
    chunk_buffer = (double*)malloc(sizeof(double)*(*chunk_size));
}

double* AudioFileReader::GetChunk()
{
    // offset the read head. Used when desired chunk size
    // is less than the MFCC/FFT length grouping amount.
    // Ex. 6 * 4096 samp = 24576. Subtract to get 22050 chunk
    if((del%MFCCGROUPING)==0 && del != 0){
        chunk_iterator-=  chnk_r_offset;
    }
    chunk_ptr = &buffer[chunk_iterator];    
    chunk_iterator <= buffer_size - (*hop_size)  ? chunk_iterator += (*hop_size)  : chunk_iterator = 0;
    del++;
    return chunk_ptr;
}

double* AudioFileReader::ChunkCopy()
{
    for(int i = 0 ; i < (*chunk_size) ; i++){
        chunkcopy[i] = chunk_buffer[i];
    }
    return chunkcopy;
}

double* AudioFileReader::GetCopy()
{
    return chunkcopy;
}

u_short_int AudioFileReader::GetChunkSize()
{
    return *chunk_size;
}

void AudioFileReader::SetHopSize(const u_short_int *s_hop_size)
{
        hop_size = const_cast<u_short_int*>(s_hop_size);
}

u_short_int AudioFileReader::GetHopSize()
{
    return *hop_size;
}

void AudioFileReader::SetChunkReadOffset(int s_chnk_r_offset)
{
    chnk_r_offset = s_chnk_r_offset;
}

