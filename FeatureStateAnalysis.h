//
//  HMMFeatureStateAnalysis.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 10/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__FeatureStateAnalysis__
#define __FFT_Chain__FeatureStateAnalysis__

#include <iostream>
#include <string>

#include "Database.h"
#include "kMeansClustering.h"
#include "StateMFCCs.h"
#include "PrintfDebug.h"
#include "GlobalTypedefs.h"

class FeatureStateAnalysis{
public:
    FeatureStateAnalysis();
    ~FeatureStateAnalysis();
    void RetrieveDatabase();
    void Cluster(const u_short_int num_clusters,const u_short_int num_iterations);
    void KdTreeCluster(const u_short_int num_clusters,const u_short_int num_iterations);
    double** GetGroupedMFCCMat();
    int GetGroupedMFCCDimen();
    int GetGroupedMFCCSize();
private:
    Database database;
    kMeansClustering kMClustering;
    StateMFCCs state_mfccs;
};
#endif /* defined(__FFT_Chain__FeatureStateAnalysis__) */
