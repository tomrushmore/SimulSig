//
//  kdTree.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 17/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__kdTree__
#define __FFT_Chain__kdTree__

#include <iostream>
#include <nanoflann.hpp>
#include "kdTreeStruct.h"
#include "PrintfDebug.h"
#include "SetupMacros.h"
// needs to changed. DIMEN not always going to be 40
#define DIMENT FNMELBANDS * MFCCGROUPING

class kdTree{
public:
    kdTree();
    ~kdTree();
    void Setup();
    void SetDimensionSize(int dim_size);
    void PopulateMFCCTree(double** mfccs,int num_coef,long num_frames);
    void PopulateCentroidTree(double** centroids,int num_coef,int num_centroids);
    int FindClosestMFCCs(double* centroid);
    int FindClosestCentroid(double* mfcc);

private:
    PointCloud pc;
    PointCloud centroid_cloud;
    
    typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloud > ,PointCloud,
    DIMENT
    > my_kd_tree_t;
    
    my_kd_tree_t *mfcc_tree;
    my_kd_tree_t *centroid_tree;
    
    std::vector<unsigned long> ret_index;
    std::vector<double> out_dist_sqr;
};

#endif /* defined(__FFT_Chain__kdTree__) */
