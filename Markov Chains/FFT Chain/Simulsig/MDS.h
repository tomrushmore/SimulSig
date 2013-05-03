//
//  MDS.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 19/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __FFT_Chain__MDS__
#define __FFT_Chain__MDS__

#include <iostream>
#include <iostream>
#include <math.h>
#include <vector>
#include <Accelerate/Accelerate.h>
#include "PrintfDebug.h"

#define W 3

typedef std::pair<double, int> eigval_index_pair;

class MDS{
public:
    MDS();
    ~MDS();
    void SetInputMatrix(double** s_input_matrix,int s_matrix_dimen,int s_matrix_row_len);
    void MakeProximitySquaredMatrix();
    void CenterMatrix();
    void DoubleCentering();
    void ApplyMDS(int n_reduced_dim);
    
    
private:
    
    void RotateMatrix(double*mat);
    double EuclidDistance(double *a,double *b);
    void MatrixMultiplty(double **a, double **b);
    void MatrixScalarMultiply(double **a,double scalar);
    double** input_matrix;
    double** proxim_matrix;
    double** temp_centering_matrix;
    double** centering_matrix_one;
    double** centering_matrix_two;
    double** double_center_matrix;
    double** identity_matrix;
    double** matrix_multiply_mat;
    int matrix_dimen;
    int matrix_row_len;
    
 
   
};

#endif /* defined(__FFT_Chain__MDS__) */
