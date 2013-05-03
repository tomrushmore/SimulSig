//
//  main.cpp
//  test
//
//  Created by Thomas Rushmore on 08/11/2012.
//  Copyright (c) 2012 Thomas Rushmore. All rights reserved.
//

#include "MainIncludes.h"

#define APPLYMDS false

enum PhaseSelect{
    Extract = 0,
    Cluster,
    Play,
    ExCl,
    MDSTest,
    Mkv,
    AllPhase
};

int main(int argc, const char * argv[])
{
    PhaseSelect Phase = ExCl;

    // Extract and write the features from the corpus to file.
    if(Phase == Extract || Phase == ExCl || Phase == AllPhase){
        DatabaseProcessor dbp(SAMPMACHINEPATH,"Analysis Corpus/",MACHINEPATH);
        bool s_cmplt = dbp.FeExSetup(FFTSIZE,HOPOVRLAPMS,NMELBANDS,SR,MINFREQ,MAXFREQ,NSAMP,MFCCGROUPING);
        if(!s_cmplt) return 0;
        dbp.Process(WRITEFILE);
    }
    if(Phase == Cluster || Phase == ExCl || Phase == AllPhase){
        clock_t start;
        double diff;
        start = clock();
        FeatureStateAnalysis cluster;
        cluster.KdTreeCluster(NUMCLUSTERS, MAXCLUSTITERATIONS);
        diff = (std::clock()-start)/(double)CLOCKS_PER_SEC;
        printf("Clustering Time  : %f seconds\n",diff);
        // MDS struggles with large matrices.
        if(APPLYMDS){
            MDS ms;
            ms.SetInputMatrix(cluster.GetGroupedMFCCMat(), cluster.GetGroupedMFCCDimen(),cluster.GetGroupedMFCCSize());
            ms.ApplyMDS(3);
        }
    }
    
    if(Phase == Mkv || Phase == AllPhase){
        // Extract the features from the training corpus.
        DatabaseProcessor dbpt(SAMPMACHINEPATH,"Training Corpus/",MACHINEPATH);
        bool s_cmplt = dbpt.FeExSetup(FFTSIZE,HOPOVRLAPMS,NMELBANDS,SR,MINFREQ,MAXFREQ,NSAMP,MFCCGROUPING);
        if(!s_cmplt) return 0;
        dbpt.Process(NWRITEFILE);
        MarkovMaster mkv_mst;
    }
    
    // Unused classic audio-mosaicing class.
    if(Phase == Play){
        Player concat_player;
        concat_player.Setup(FFTSIZE,HOPOVRLAPMS);
        concat_player.Play();
    }
    // Used to test MDS function against known results
    if(Phase == MDSTest){
        int s = 3;
        double matrix[3*3] = {1,5,10,
                              6,4,2,
                              1,2,10};
        double** ma = (double**)malloc(sizeof(double*)*s);
        for(int i = 0 ; i < s ; i++){
            ma[i] = (double*)malloc(sizeof(double)*s);
        }
        for(int i = 0 ; i < s; i++){
            for(int j = 0; j < s; j++){
                ma[i][j] = matrix[i*s+j];
            }
        }
        MDS ms;
        ms.SetInputMatrix(ma, s,s);
        ms.ApplyMDS(2);
    }
    return 0;
}
