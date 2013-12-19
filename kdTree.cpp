//
//  kdTree.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 17/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "kdTree.h"

kdTree::kdTree()
{
    centroid_tree = NULL;
    ret_index.resize(5);
    out_dist_sqr.resize(5);
}

kdTree::~kdTree()
{
  //  delete mfcc_tree;
  //  delete centroid_tree;
}

void kdTree::Setup()
{
    
}

void kdTree::PopulateMFCCTree(double **mfccs, int num_coef, long num_frames)
{
    debug_print(("\nPopulating MFCC Tree..."));
    pc.add_database(mfccs, num_coef, num_frames);
    debug_print(("\tComplete\n"));
    debug_print(("Indexing Tree.\n"));
 
    mfcc_tree = new my_kd_tree_t(num_coef/*dim*/, pc, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
	mfcc_tree->buildIndex();
    debug_print(("Index complete.\n"));
}

void kdTree::PopulateCentroidTree(double **centroids, int num_coef,int num_centroids)
{
    if(centroid_tree != NULL){
        delete (centroid_tree);
        centroid_cloud.pts.empty();
    }
    debug_print(("Populating Centroid Tree..."));
    centroid_cloud.add_database(centroids, num_coef, num_centroids);
    debug_print(("\tComplete\n"));
    debug_print(("Indexing Tree.\n"));
    
    centroid_tree = new my_kd_tree_t(num_coef/*dim*/, centroid_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
	centroid_tree->buildIndex();
    debug_print(("Index complete.\n"));
}

int kdTree::FindClosestMFCCs(double *centroid)
{
    mfcc_tree->knnSearch(centroid, 5, &ret_index[0], &out_dist_sqr[0]);
    for(int i = 0 ; i < 5; i++){
       // printf("%f ",pc.kdtree_get_orig_index(ret_index[i]));
    }
    return pc.kdtree_get_orig_index(ret_index[0]);
}

int kdTree::FindClosestCentroid(double *mfcc)
{
    centroid_tree->knnSearch(&mfcc[0], 1, &ret_index[0], &out_dist_sqr[0]);
    return centroid_cloud.kdtree_get_orig_index(ret_index[0]);
}