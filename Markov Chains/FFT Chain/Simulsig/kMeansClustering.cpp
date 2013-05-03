//
//  Clustering.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 09/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "kMeansClustering.h"

kMeansClustering::kMeansClustering()
{
    centroids = NULL;
    mfcc_matrix = NULL;
    dimensions_max_val = NULL;
    paired_centroid = NULL;
    grouping = 1;
    num_coef_grouped = 1;
    kdtree.Setup();
}

kMeansClustering::~kMeansClustering()
{
    if(centroids){
        printf("Freeing Cluster Matrix...\n");
        for(int i = 0 ; i < num_clusters; i++){
            free(centroids[i]);
        }
        free(centroids);
        printf("Freed.\n");
    }
    if(paired_centroid){
        printf("Freeing paired centroid vector...\n");
        free(paired_centroid);
        printf("Freed.\n");
    }
    if(dimensions_max_val){
        printf("Freeing dimensions_max_val...\n");
        free(dimensions_max_val);
        printf("Freed.\n");
    }
}

void kMeansClustering::SetClusterParameters(int s_num_clusters,int s_num_iterations,long s_num_mfccs,int s_mfcc_num_coef,int s_grouping)
{ 
    num_clusters = s_num_clusters;
    max_num_iterations = s_num_iterations;
    num_mfccs = s_num_mfccs;
    mfcc_num_coef = s_mfcc_num_coef;
    
    grouping = s_grouping;
    num_coef_grouped = mfcc_num_coef * grouping;
    num_fram_grouped = floor(num_mfccs / grouping);
    printf("COEF GROP : %d\n",num_coef_grouped);
    printf("FRAM GROP : %d\n",num_fram_grouped);
    printf("CLUS GROP : %d\n",num_clusters);

    centroids = (double**)malloc(sizeof(double*)*num_clusters);
    for(int i = 0; i < num_clusters; i++){
        centroids[i] = (double*)malloc(sizeof(double)*num_coef_grouped);
        memset(centroids[i],0,sizeof(double)*num_coef_grouped);
    }
    paired_centroid = (int*)malloc(sizeof(int*)*num_fram_grouped);
    memset(paired_centroid, 0, sizeof(int)*num_fram_grouped);
    paired_vector.resize(num_fram_grouped);
    
}

void kMeansClustering::KdTreeCluster(double **s_mfcc_matrix)
{
    mfcc_matrix = s_mfcc_matrix;
    clock_t start;
    double diff;
    start = clock();
    
    kdtree.PopulateMFCCTree(mfcc_matrix, num_coef_grouped, num_fram_grouped);
    diff = (std::clock()-start)/(double)CLOCKS_PER_SEC;
    printf("Kd-tree Index Time : %f\n",diff);
    CalcMaxDimensionVal();
    InitializeCentroids();
    
    int* mfcc_p_cluster = (int*)malloc(sizeof(int)*num_clusters);
    double** sum = (double**)malloc(sizeof(double)*num_clusters);
    for(int p = 0; p < num_clusters; p++){
        sum[p] = (double*)malloc(sizeof(double)*num_coef_grouped);
        memset(sum[p], 0.0, sizeof(double)*num_coef_grouped);
    }
    for(int i = 0 ; i < max_num_iterations; i++){
        kdtree.PopulateCentroidTree(centroids, num_coef_grouped, num_clusters);
        for(int p = 0; p < num_clusters; p++){
            memset(sum[p], 0.0, sizeof(double)*num_coef_grouped);
        }
        memset(mfcc_p_cluster, 0, sizeof(int)*num_clusters);

        int closest_centroid = 0;
        //debug_print(("Finding closest centroids %d ...\n",i));

        for(int p = 0; p < num_fram_grouped ; p++){
            closest_centroid = kdtree.FindClosestCentroid(mfcc_matrix[p]);
            //printf("Closet Centroid : %d  \n",closest_centroid);
            paired_centroid[p] = closest_centroid;
            for(int q = 0 ; q < num_coef_grouped; q++){
                sum[closest_centroid][q] += mfcc_matrix[p][q];
            }
            mfcc_p_cluster[closest_centroid]++;
        }

        double change_check = 0.0;
        for(int j = 0 ; j < num_clusters; j++){
            for(int p = 0; p < num_coef_grouped; p++){
                if(mfcc_p_cluster[j]){
                    sum[j][p] /= mfcc_p_cluster[j];
                    change_check += (centroids[j][p] - sum[j][p]);
                    centroids[j][p] = sum[j][p];
                }
//                change_check += (centroids[j][p] - sum[j][p]);
//                centroids[j][p] = sum[j][p];
            }
        }
        if(change_check == 0.0 || i == max_num_iterations-1){
            printf("Optimal centroids found after %d iterations.\n",i);
            closest_index = (int*)malloc(sizeof(int)*num_clusters);
            memset(closest_index, 0,sizeof(int)*num_clusters);
            std::set<int> active_clusters;
            
            for(int p = 0; p < num_clusters; p++){
               int sassy =  kdtree.FindClosestMFCCs(centroids[p]);
                closest_index[p] = sassy;
                active_clusters.insert(sassy);
            }
            std::set<int>::iterator it;
            int idx = 0;
            for (it=active_clusters.begin(); it!=active_clusters.end(); ++it){
                closest_index[idx] = *it;
                //printf("ACTIVE CLUST : %d \n",*it);
                idx++;
            }
            //printf("Found best fit MFCC's for each centroid.\n");
            n_active_clusters = (int)active_clusters.size();
            printf("Number active clusters : %d \n",n_active_clusters);
            for(int a = 0; a < num_clusters; a++){
                free(sum[a]);
            }
            free(sum);
            free(mfcc_p_cluster);
            return;
        }
    }
    for(int a = 0; a < num_clusters; a++){
        free(sum[a]);
    }
    free(sum);
    free(mfcc_p_cluster);
}

// i want to cluster all of the mfccs to a number of centroids.
// the mean ( but actually closest mfcc to the mean ) will be used
// as a markov chain state.
// i must keep the index of the closest to mean mfcc so I can retrieve
// the corresponding FFT for later inversion.

void kMeansClustering::Cluster(double **s_mfcc_matrix)
{
    mfcc_matrix = s_mfcc_matrix;
    CalcMaxDimensionVal();
    InitializeCentroids();
    
    int* mfcc_p_cluster = (int*)malloc(sizeof(int)*num_clusters);
    double** sum = (double**)malloc(sizeof(double)*num_clusters);
    for(int p = 0; p < num_clusters; p++){
        sum[p] = (double*)malloc(sizeof(double)*mfcc_num_coef);
        memset(sum[p], 0.0, sizeof(double)*mfcc_num_coef);
    }
    for(int i = 0 ; i < max_num_iterations; i++){
        // for each n vector, calculate the distance to every initialized
        // centroid. hold the minimum and assign it to it via
        // the paired_centroid matrix
        double dist,tmp_dist;
        int closest_centroid = 0;
        
        for(int p = 0; p < num_mfccs/MFCCGROUPING ; p++){
            dist = EuclidDistance(mfcc_matrix[p], centroids[0]);
            closest_centroid = 0;
            for(int q = 1 ; q < num_clusters; q++){
                tmp_dist = EuclidDistance(mfcc_matrix[p], centroids[q]);
                if(tmp_dist < dist){
                    dist = tmp_dist;
                    closest_centroid = q;
                }
            }
            paired_centroid[p] = closest_centroid;
            paired_vector[p] = closest_centroid;
        }
        for(int p = 0; p < num_clusters; p++){
            memset(sum[p], 0.0, sizeof(double)*mfcc_num_coef);
        }
        memset(mfcc_p_cluster, 0, sizeof(int)*num_clusters);
       
        for(int j = 0; j < num_mfccs/MFCCGROUPING; j++){
           // int sum_index = paired_centroid[j];
            int sum_index = paired_vector[j];
            for(int p = 0 ; p < mfcc_num_coef; p++){
                sum[sum_index][p] += mfcc_matrix[j][p];
            }
            mfcc_p_cluster[sum_index]++;
        }
        double change_check = 0.0;
        for(int j = 0 ; j < num_clusters; j++){
            for(int p = 0; p < mfcc_num_coef; p++){
                if(mfcc_p_cluster[j]){
                    sum[j][p] /= mfcc_p_cluster[j];
                }
                change_check += (centroids[j][p] - sum[j][p]);
                centroids[j][p] = sum[j][p];
            }
        }
        
        if(change_check == 0.0 || i == max_num_iterations -1){
            printf("Optimal centroids found after %d iterations.\n",i);
            CalcBestFitMFCC();
            for(int a = 0; a < num_clusters; a++){
                free(sum[a]);
            }
            free(sum);
            free(mfcc_p_cluster);
            return;
        }
    }
    for(int a = 0; a < num_clusters; a++){
        free(sum[a]);
    }
    free(sum);
    free(mfcc_p_cluster);
}

void kMeansClustering::CalcMaxDimensionVal()
{
   debug_print(("Calculating max dimensional val...\n"));
    dimensions_max_val = (double*)malloc(sizeof(double)*num_coef_grouped);
    memset(dimensions_max_val, 0, sizeof(double)*num_coef_grouped);
    dimensions_offset_val = (double*)malloc(sizeof(double)*num_coef_grouped);
    memset(dimensions_offset_val, 0, sizeof(double)*num_coef_grouped);
    double max;
    double min;
    for(int i = 0 ; i < num_coef_grouped; i++){
        max = mfcc_matrix[0][i];
        min = mfcc_matrix[0][i];
        for(int j = 1 ; j < num_fram_grouped; j++){
            if(j ==1)
            if(mfcc_matrix[j][i] > max){
                max = mfcc_matrix[j][i];
            }
            if(mfcc_matrix[j][i] < min){
                min = mfcc_matrix[j][i];
            }
        }
        double tot = max + abs(min);
        dimensions_max_val[i] = tot;
        dimensions_offset_val[i] = abs(min);
        //debug_print(("%f \n",dimensions_max_val[i]));
    }
    debug_print(("\tComplete\n"));
}

void kMeansClustering::InitializeCentroids()
{
    debug_print(("Initializing centroids..."));
    // for now, randomly select centroids.
    // future look at making sure centroids
    // are distributed evenly across the hyperspace.
    if(dimensions_max_val && centroids){
        double val = 0.0;
        for(int i = 0 ; i < num_clusters; i++){
            for(int j = 0 ; j < num_coef_grouped ; j++){
                val = (double)rand()/RAND_MAX;
                val *= dimensions_max_val[j];
                val -= dimensions_offset_val[j];
                centroids[i][j] = val;
            }
            val = 0.0;
        }
    }
    debug_print(("\tComplete\n"));

    //printf("Centroid random initialization \n");
    for(int i = 0 ; i < num_clusters; i++){
        for(int j = 0 ; j < mfcc_num_coef; j++){
     //       printf("%f ",centroids[i][j]);
        }
   //     printf("\n");
    }
  // printf("\n");
}

void kMeansClustering::CalcBestFitMFCC()
{
    // already have the centroids at the best point possible.
    // for each centroid, find the closest point, using pair?
    // pairedcentroid?
    for(int i = 0 ; i < num_clusters;i++){
        debug_print(("Paired Centroid : %d \n",paired_centroid[i]));
    }
    double* closest = (double*)malloc(sizeof(double)*num_clusters);
    for(int p = 0; p < num_clusters; p++){
        closest[p] = 1000;
    }
    closest_index = (int*)malloc(sizeof(int)*num_clusters);
    memset(closest_index, 0,sizeof(int)*num_clusters);
    
    double d = 0.0;
    int ind = 0;
    for(int i = 0 ; i < num_fram_grouped; i++){
        ind = paired_centroid[i];
        d = EuclidDistance(mfcc_matrix[i],centroids[ind]);
        if(d < closest[ind]){
            closest[ind] = d;
            closest_index[ind] = i;
        }
    }
    n_active_clusters = 0;
    for(int i = 0; i < num_clusters; i++){
        debug_print(("Closest MFCC to Centroid %d : %d \n",i,closest_index[i]));
        if(closest_index[i] > 0){
            n_active_clusters++;
        }
    }
   // printf("Found best fit MFCC's for each centroid.\n");
    printf("Num active out of %d potential clusters : %d\n",num_clusters,n_active_clusters);
}

int* kMeansClustering::GetBestFitIndices()
{
    return closest_index;
}

int kMeansClustering::GetNActiveClusters()
{
    return n_active_clusters;
}

double kMeansClustering::EuclidDistance(double *a, double *b)
{
    double running_sum = 0.0;
    for(int i = 0 ; i < mfcc_num_coef; i++){
        running_sum += pow((a[i] - b[i]),2);
    }
    running_sum = sqrt(running_sum);
    return running_sum;
}

int kMeansClustering::GetNumCoefGrouped()
{
    return num_coef_grouped;
}
