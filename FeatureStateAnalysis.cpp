//
//  FeatureStateAnalysis.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 10/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "FeatureStateAnalysis.h"

FeatureStateAnalysis::FeatureStateAnalysis()
{

}

FeatureStateAnalysis::~FeatureStateAnalysis()
{
    printf("\n*-------Clustering Phase Destructor------*\n");
}

void FeatureStateAnalysis::KdTreeCluster(const u_short_int num_clusters, const u_short_int num_iterations)
{
    printf("\n*-----------------Clustering Phase-----------------*\n");
    RetrieveDatabase();
    kMClustering.SetClusterParameters(num_clusters,num_iterations,database.GetNumFrames(), database.GetNumMFCCCoefs(),MFCCGROUPING);
    kMClustering.KdTreeCluster(database.GetMFCCMatrix());
    state_mfccs.SetupWrite(kMClustering.GetNActiveClusters(), database.GetNumMFCCCoefs()*MFCCGROUPING, database.GetFFTLength(),database.GetNumTracks());
    state_mfccs.SetMFCCFFTPairs(database.GetMFCCMatrix(), database.GetFFTMatrix(), kMClustering.GetBestFitIndices());
    state_mfccs.WriteMFCCFFTPairs();
    printf("*------------Clustering Phase Complete-------------*\n\n");
}

void FeatureStateAnalysis::Cluster(const u_short_int num_clusters,const u_short_int num_iterations)
{
    RetrieveDatabase();
    kMClustering.SetClusterParameters(num_clusters,num_iterations,database.GetNumFrames(), database.GetNumMFCCCoefs(),MFCCGROUPING);
    kMClustering.Cluster(database.GetMFCCMatrix());
    state_mfccs.SetupWrite(kMClustering.GetNActiveClusters(), database.GetNumMFCCCoefs(), database.GetFFTLength(),database.GetNumTracks());
    state_mfccs.SetMFCCFFTPairs(database.GetMFCCMatrix(), database.GetFFTMatrix(), kMClustering.GetBestFitIndices());
    state_mfccs.WriteMFCCFFTPairs();
}

void FeatureStateAnalysis::RetrieveDatabase()
{
    database.ReadFileGroup(0,MFCCGROUPING,0);
}

double** FeatureStateAnalysis::GetGroupedMFCCMat()
{
    return database.GetMFCCMatrix();
}

int FeatureStateAnalysis::GetGroupedMFCCDimen()
{
    return kMClustering.GetNumCoefGrouped();
}

int FeatureStateAnalysis::GetGroupedMFCCSize()
{
    return database.GetGroupedS();
}
