//
//  MarkovMaster.cpp
//  Simulsig
//
//  Created by Thomas Rushmore on 27/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "MarkovMaster.h"

MarkovMaster::MarkovMaster()
{
    SetStates();
    Train();
    PlayChain();
}

MarkovMaster::~MarkovMaster()
{
    
}

void MarkovMaster::SetStates()
{
    state_mfccs.ReadMFCCFFTStatePairs();
    mkv_ch.SetStates(state_mfccs.GetNumActiveClusters(), state_mfccs.GetMFCCMatrix(), state_mfccs.GetNumMFCCCoef());
}

void MarkovMaster::Train()
{
    database.ReadFileGroup(0,state_mfccs.GetGrouping(),1);
    mfcc_train_matrix = database.GetMFCCMatrix();
    mkv_ch.StdBuildChainMatrix();
    for(int i = 0 ; i < database.GetNumTracks();i++){
        printf("Track Frame Length %d : %d \n",i,database.GetTrackFrameLengths(i));
    }
    CalcHoppedTrackLen(database.GetNumTracks());

    int offset = 0;
    int group_length = 0;
    printf("Longest Track Length : %d\n",database.GetLongestTrackLength());
    for(int i = 0 ; i < database.GetNumTracks(); i++){
        group_length = trk_n_frame_hop[i];
        for(int j = 0; j < group_length - 1 ; j++){
            mkv_ch.StdstepTrain(mfcc_train_matrix[j+offset], mfcc_train_matrix[j+offset+1]);
        }
        offset += group_length;
        mkv_ch.reset();
    }
    //mkv_ch.MatrixPrinter();
    //mkv_ch.PrintMatrix();
}

void MarkovMaster::PlayChain()
{  
    mkv_pl.Setup(&mkv_ch,state_mfccs.GetFFTMatrix(),NUMMKVGEN,max_h_len,state_mfccs.GetNumMFCCCoef(),state_mfccs.GetNumActiveClusters(),state_mfccs.GetFFTLength(),state_mfccs.GetGrouping(),state_mfccs.GetHopOverlap(),state_mfccs.GetFFTLength());
    mkv_pl.GenerateSampleBuffers();
    mkv_pl.GenerateChainSequences(NHOMOG);
    mkv_pl.GenerateSignals();
    mkv_pl.WriteFinalSignals();
}

void MarkovMaster::CalcHoppedTrackLen(int num_tracks)
{
    int frame_size = (database.GetFFTLength() -2) + (state_mfccs.GetGrouping()-1) * ((database.GetFFTLength()-2) - state_mfccs.GetHopOverlap() + 2) ;
    trk_n_frame_hop.push_back(floor(database.GetBufSizes(0)/frame_size));
    int max = trk_n_frame_hop[0];
    for(int i = 1 ; i < num_tracks; i++){
        trk_n_frame_hop.push_back(floor(database.GetBufSizes(i)/frame_size));
        if(trk_n_frame_hop[i] > max){
            max = trk_n_frame_hop[i];
        }
    }
    max_h_len = max;
}


