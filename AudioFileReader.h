//
//  AudioFileReader.h
//  test
//
//  Created by Thomas Rushmore on 10/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __test__AudioFileReader__
#define __test__AudioFileReader__

#include <iostream>

#include "sndfile.hh"
#include "PrintfDebug.h"
#include "GlobalTypedefs.h"
#include "SetupMacros.h"

class AudioFileReader{
public:
    AudioFileReader();
    ~AudioFileReader();
    int Load(const char * s_file_path);
    int WriteFile(const char * s_file_path);
    double * GetBuffer();
    double GetSample(int index);
    void SetChunkSize(const u_short_int *s_chunk_size);
    void SetChunkReadOffset(int s_chnk_r_offset);
    u_short_int GetChunkSize();
    double * GetChunk();
    double* ChunkCopy();
    double* GetCopy();
    void SetHopSize(const u_short_int *s_hop_size);
    u_short_int GetHopSize();
    long GetBufferSize();
    
private:
    char * file_path;
    double * buffer;
    long buffer_size;
    u_short_int *chunk_size;
    int chunk_iterator;
    double * chunk_buffer;
    u_short_int *hop_size;
    double *chunkcopy;
    double * chunk_ptr;
    int chnk_r_offset;

    int del;
};
#endif /* defined(__test__AudioFileReader__) */
