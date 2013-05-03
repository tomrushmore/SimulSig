//
//  MDS.cpp
//  FFT Chain
//
//  Created by Thomas Rushmore on 19/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "MDS.h"

extern bool sort_pred(const eigval_index_pair& left, const eigval_index_pair& right);

MDS::MDS()
{
    
}

MDS::~MDS()
{
    
}

void MDS::SetInputMatrix(double **s_input_matrix,int s_matrix_dimen,int s_matrix_row_len)
{
    input_matrix = s_input_matrix;
    matrix_dimen = s_matrix_dimen;
    matrix_row_len = s_matrix_row_len;
    // this needs to be put in.
    proxim_matrix = (double**)malloc(sizeof(double*)*matrix_dimen);
    temp_centering_matrix = (double**)malloc(sizeof(double*)*matrix_dimen);
    centering_matrix_one = (double**)malloc(sizeof(double*)*matrix_dimen);
    centering_matrix_two = (double**)malloc(sizeof(double*)*matrix_dimen);

    double_center_matrix =(double**)malloc(sizeof(double*)*matrix_dimen);
    identity_matrix = (double**)malloc(sizeof(double*)*matrix_dimen);
    matrix_multiply_mat = (double**)malloc(sizeof(double*)*matrix_dimen);
    for(int i = 0 ; i < matrix_dimen; i++){
        proxim_matrix[i] = (double*)malloc(sizeof(double)*matrix_dimen);
        temp_centering_matrix[i] = (double*)malloc(sizeof(double)*matrix_dimen);
        centering_matrix_one[i] = (double*)malloc(sizeof(double)*matrix_dimen);
        centering_matrix_two[i] = (double*)malloc(sizeof(double)*matrix_dimen);
        double_center_matrix[i] = (double*)malloc(sizeof(double)*matrix_dimen);
        identity_matrix[i] = (double*)malloc(sizeof(double)*matrix_dimen);
        matrix_multiply_mat[i] = (double*)malloc(sizeof(double)*matrix_dimen);
    }
    
}

void MDS::MakeProximitySquaredMatrix()
{
    // first calculate proximity matr
    for(int i = 0 ; i < matrix_dimen; i++){
        for(int j = 0 ; j < matrix_dimen; j++){
            proxim_matrix[i][j] = pow(EuclidDistance(input_matrix[i], input_matrix[j]),2);
            
        }
    }
    debug_print(("Proximity Matrix\n"));
    for(int i = 0 ; i < matrix_dimen; i++){
        for(int j = 0 ; j < matrix_dimen; j++){
          //  proxim_matrix[i][j] =  pow(proxim_matrix[i][j],2);
            debug_print(("%f ",proxim_matrix[i][j]));
        }
        debug_print(("\n"));
    }
}

void MDS::CenterMatrix()
{
    for(int i = 0 ; i < matrix_dimen; i++){
        for(int j = 0; j < matrix_dimen; j++){
//            if(i == j){
//                identity_matrix[i][j] = 1.0;
//            } else {
//                identity_matrix[i][j] = 0.0;
//            }
            i == j ? identity_matrix[i][j] = 1.0 : identity_matrix[i][j] = 0.0;
        }
    }
    debug_print(("Identity Matrix\n"));
    for(int i = 0 ; i < matrix_dimen; i++){
        for(int j = 0; j < matrix_dimen; j++){
            debug_print(("%f ",identity_matrix[i][j]));
        }
        debug_print(("\n"));
    }
    debug_print(("\n\n"));
    
    double scalar = 1.0 / matrix_dimen;
    debug_print(("Centering matrix before Id multiplication\n"));
    for(int i = 0 ; i < matrix_dimen; i++){
        for(int j = 0 ; j < matrix_dimen; j++){
            temp_centering_matrix[i][j] = scalar;
            debug_print(("%f ",temp_centering_matrix[i][j]));
        }
        debug_print(("\n"));
    }
    debug_print(("\n"));
    debug_print(("Final centering Matrix\n"));
    double id_m_val = 0;
    for(int i = 0 ; i < matrix_dimen; i++){
        for(int j = 0 ; j < matrix_dimen; j++){
            i == j ? id_m_val = 1.0 : id_m_val = 0.0;
            centering_matrix_one[i][j] = id_m_val - scalar;
            centering_matrix_two[i][j] = centering_matrix_one[i][j] * -0.5;

            debug_print(("%f ",centering_matrix_one[i][j]));
        }
        debug_print(("\n"));
    }
    debug_print(("\n"));
    debug_print(("Final centering Matrix Two \n"));
    for(int i = 0 ; i < matrix_dimen; i++){
        for(int j = 0 ; j < matrix_dimen; j++){
            debug_print(("%f ",centering_matrix_two[i][j]));
        }
        debug_print(("\n"));
    }
    debug_print(("\n"));
}

void MDS::DoubleCentering()
{
    for(int i = 0; i < matrix_dimen; i++){
        for(int j = 0; j < matrix_dimen; j++){
            double_center_matrix[i][j] = proxim_matrix[i][j];
        }
    }
    MatrixMultiplty(centering_matrix_two, double_center_matrix);
    debug_print(("\nDouble Centered Matrix \n"));
    MatrixMultiplty(centering_matrix_two, centering_matrix_one);
    for(int i = 0; i < matrix_dimen; i++){
        for(int j = 0; j < matrix_dimen; j++){
            double_center_matrix[i][j] = centering_matrix_two[i][j];
        }
    }
}

void MDS::ApplyMDS(int n_reduced_dim)
{
    MakeProximitySquaredMatrix();
    CenterMatrix();
    DoubleCentering();
    
    double* work_matrix = (double*)malloc(sizeof(double*)*matrix_dimen*matrix_dimen);
    for(int i = 0 ; i < matrix_dimen;i++){
        for(int j = 0 ; j <matrix_dimen; j++){
            work_matrix[i+j*matrix_dimen] = double_center_matrix[i][j];
        }
    }
    RotateMatrix(work_matrix);

    // Ref'd in dissertation
    __CLPK_integer  lda,ldvl,ldvr,lwork,info;
    char job = 'V';
    char jobl = 'V';
    __CLPK_integer no = matrix_dimen;
    lda = no;
    ldvl = no;
    ldvr = no;
    info = 0;
    lwork = pow(matrix_dimen,3);
    
    double* wr = (double*)malloc(sizeof(double)*matrix_dimen);
    double* wi = (double*)malloc(sizeof(double)*matrix_dimen);
    double* v1 = (double*)malloc(sizeof(double)*pow(matrix_dimen,2)*2);
    double* vr = (double*)malloc(sizeof(double)*pow(matrix_dimen,2));
    double* dummy = (double*)malloc(sizeof(double)*matrix_dimen);
    memset(wr, 0, sizeof(double)*matrix_dimen);
    memset(wi, 0, sizeof(double)*matrix_dimen);
    memset(v1, 0, sizeof(double)*pow(matrix_dimen,2)*2);
    memset(vr, 0, sizeof(double)*pow(matrix_dimen,2));
    memset(dummy, 0, sizeof(double)*matrix_dimen);

    dgeev_(&jobl, &job, &no, work_matrix, &lda, wr, wi, v1, &ldvl, vr, &ldvr, dummy, &lwork, &info);
    
    if (info>=0)
    {
        debug_print(("\e-vals = "));
        for(int i = 0 ; i < matrix_dimen; i++){
            if(wi[i] == (double)0.0){
                debug_print(("%f , ",wr[i]));
            } else {
                debug_print((" ( %f , %f )",wr[i],wi[i]));
            }
        }
        debug_print(("\n\n"));
        
        struct Eigenvector{
            std::vector<std::pair<double, double>> vec;
            std::vector<double> realvec;
        };
        std::vector<Eigenvector> evectors;
        evectors.resize(matrix_dimen);
        
        int j;
        for(int i = 0; i < matrix_dimen; i++ ) {
            j = 0;
            while( j < matrix_dimen ) {
                if( wi[j] == (double)0.0 ) {
                    debug_print(("%f ", v1[i+j*matrix_dimen]));
                    evectors[j].vec.push_back(std::make_pair(v1[i+j*matrix_dimen], 0.0));
                    evectors[j].realvec.push_back(v1[i+j*matrix_dimen]);
                    j++;
                } else {
                    debug_print((" %f,%f ", v1[i+j*matrix_dimen], v1[i+(j+1)*matrix_dimen]));
                    evectors[j].vec.push_back(std::make_pair(v1[i+j*matrix_dimen], v1[i+(j+1)*matrix_dimen]));
                    debug_print((" %f,%f ", v1[i+j*matrix_dimen], -v1[i+(j+1)*matrix_dimen]));
                    evectors[j+1].vec.push_back(std::make_pair(v1[i+j*matrix_dimen], -v1[i+(j+1)*matrix_dimen]));
                    j += 2;
                }
            }
            debug_print(("\n"));
        }
        for(int i = 0 ; i < evectors[0].vec.size(); i++){
            debug_print(("e-vecs %d : [ ",i));
            for(int j = 0; j < evectors[0].vec.size(); j++){
                debug_print(("%f ",evectors[i].realvec[j]));
            }
            debug_print(("]\n"));
        }
        std::vector<eigval_index_pair> eigval_ind_pairs;
        for(int i = 0 ; i < matrix_dimen; i++){
            eigval_ind_pairs.push_back(eigval_index_pair(wr[i],i));
            debug_print(("%f ",eigval_ind_pairs[i].first));
        }
        debug_print(("\n\n"));
        std::sort(eigval_ind_pairs.begin(),eigval_ind_pairs.end(),sort_pred);
        eigval_ind_pairs.resize(n_reduced_dim);

        for(int i = 0 ; i < n_reduced_dim; i++){
            debug_print(("%f ",eigval_ind_pairs[i].first));
            debug_print(("%d ",eigval_ind_pairs[i].second));
            debug_print(("\n"));
        }
        
        double * eig_mat = (double*)malloc(sizeof(double)*n_reduced_dim*matrix_dimen);
        memset(eig_mat, 0, sizeof(double)*n_reduced_dim*matrix_dimen);
        double * eig_val = (double*)malloc(sizeof(double)*n_reduced_dim*matrix_dimen);
        memset(eig_val, 0, sizeof(double)*n_reduced_dim*matrix_dimen);
        
        for(int i = 0; i < n_reduced_dim; i++){
            for(int j = 0; j < matrix_dimen; j++){
                int sel = eigval_ind_pairs[i].second;
                debug_print(("%d ",sel));
                eig_mat[i+j*n_reduced_dim] = evectors[sel].realvec[j];
            }
        }
        for(int i = 0; i < n_reduced_dim; i++){
            for(int j = 0; j < n_reduced_dim; j++){
                if(i == j){
                    eig_val[i+j*n_reduced_dim] =wr[eigval_ind_pairs[i].second];
                } else {
                    eig_val[i+j*n_reduced_dim] = 0.0;
                }
            }
        }
        debug_print(("\n\nMMAT\n\n"));
        for(int i = 1; i < n_reduced_dim*j+1; i++){
            debug_print(("%f ",eig_mat[i-1]));
            if((i % n_reduced_dim)==0 ){
                debug_print(("\n"));
            }
        }
        debug_print(("\n\n"));
        for(int i = 1; i < n_reduced_dim*n_reduced_dim+1; i++){
            debug_print(("%4.16f ",eig_val[i-1]));
            if((i % n_reduced_dim)==0 ){
                debug_print(("\n"));
            }
        }
        // result matrix
        debug_print(("\n"));
        double * reduced_matrices = (double*)malloc(sizeof(double)*n_reduced_dim*matrix_dimen);
        memset(reduced_matrices, 0, sizeof(double)*n_reduced_dim*matrix_dimen);
        for(int i = 0; i < matrix_dimen; i++){
            for(int j = 0; j < n_reduced_dim ; j++){
                debug_print(("%f ",eig_mat[i*n_reduced_dim+j]));
                reduced_matrices[i*n_reduced_dim+j] = eig_mat[i*n_reduced_dim+j] * sqrt(eig_val[j*n_reduced_dim+j]);
            }
            debug_print(("\n"));
        }
        debug_print(("\n"));
        
        printf("Final reduction\n");
        for(int i = 1; i < n_reduced_dim*matrix_dimen+1; i++){
            printf("%f ",reduced_matrices[i-1]);
            if((i % n_reduced_dim)==0 ){
                printf("\n");
            }
        }
        free (eig_mat);
        free (eig_val);
        free (reduced_matrices);
        free (wr);
        free (wi);
        free (v1);
        free (vr);
        free (dummy);
        free(work_matrix);
    }
}

void MDS::MatrixScalarMultiply(double **a, double scalar)
{
    debug_print(("Scalar multiply\n"));
    for(int i = 0 ;i < matrix_dimen; i++){
        for(int j = 0 ; j < matrix_dimen; j++){
            a[i][j] *= scalar;
            debug_print(("%f ",a[i][j]));
        }
        debug_print(("\n"));
    }
}

void MDS::MatrixMultiplty(double **a, double **b)
{
    debug_print(("Matrix multiply\n"));
    for(int i = 0; i < matrix_dimen; i++){
        for(int j = 0; j < matrix_dimen; j++){
            double sum = 0.0;
            for(int p = 0; p < matrix_dimen; p++){
                sum += a[p][j] * b[i][p];
            }
            matrix_multiply_mat[i][j] = sum;
        }
    }
    for(int i = 0 ; i < matrix_dimen; i++){
        for(int j = 0 ; j < matrix_dimen; j++){
            a[i][j] = matrix_multiply_mat[i][j];
            debug_print(("%f ",a[i][j]));
        }
        debug_print(("\n"));
    }
    debug_print(("\n"));
}

void MDS::RotateMatrix(double *mat)
{
    for(int i = 0 ; i < matrix_dimen ; i++){
            for(int j = i + 1; j < matrix_dimen;j++) {
                double a = mat[i*matrix_dimen+j];
                double b = mat[i+j*matrix_dimen];
                mat[i*matrix_dimen+j] = b;
                mat[i+j*matrix_dimen] = a;
            }
        }
}

double MDS::EuclidDistance(double *a, double *b)
{
    double running_sum = 0.0;
    for(int i = 0 ; i < matrix_dimen; i++){
        running_sum += pow((a[i] - b[i]),2);
    }
    running_sum = sqrt(running_sum);
    return running_sum;
}

bool sort_pred(const eigval_index_pair& left, const eigval_index_pair& right)
{
    return left.first > right.first;
}
