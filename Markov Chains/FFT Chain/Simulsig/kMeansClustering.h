//
//  Clustering.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 09/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__kMeansClustering__
#define __FFT_Chain__kMeansClustering__

#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "kdTree.h"
#include <set>
#include "PrintfDebug.h"

class kMeansClustering{
public:
    kMeansClustering();
    ~kMeansClustering();
    void SetClusterParameters(int s_num_clusters,int s_num_iterations,long s_num_mfccs,int s_mfcc_num_coef,int s_grouping);
    void Cluster(double** s_mfcc_matrix);
    void KdTreeCluster(double** s_mfcc_matrix);
    void CalcBestFitMFCC();
    int* GetBestFitIndices();
    int GetNActiveClusters();
    int GetNumCoefGrouped();

private:
    int grouping;
    int num_coef_grouped;
    int num_fram_grouped;
    
    double EuclidDistance(double* a,double* b);
    void CalcMaxDimensionVal();
    void InitializeCentroids();
    int num_clusters;
    int max_num_iterations;
    long num_mfccs;
    int mfcc_num_coef;
    
    double** centroids;
    double** mfcc_matrix;
    int* paired_centroid;
    std::vector<int> paired_vector;
    double* dimensions_max_val;
    double* dimensions_offset_val;
    int* closest_index;
    int n_active_clusters;
    
    std::vector<int>* pairings;
    
    kdTree kdtree;
    
    
};

#endif /* defined(__FFT_Chain__kMeansClustering__) */
