//
//  Markov Chain.h
//  Simulsig
//
//  Created by Thomas Rushmore on 27/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __Simulsig__Markov_Chain__
#define __Simulsig__Markov_Chain__

#include <math.h>
#include "tnt.h"
#include "PrintfDebug.h"
#include "DBWrite.h"
#include "FileWriteR.h"


#if 1
    #define FCOEF 1
    #define FCOEFP1 2
#else
    #define FCOEF 0
    #define FCOEFP1 1

#endif

class MarkovChain{
public:
    MarkovChain();
    ~MarkovChain();
    void BuildChain();
    void SetStates(int set_number_states,double** mfcc_states,int s_num_mfcc_coef);
    void BuildChainMatrix(int set_depth);
    void BuildChainZeroed(int set_number_states,int set_depth);
    void BuildChainNormalized(int set_number_states,int set_depth);
    void CheckTotalProbability();
    void PrintMatrix();
    void CDF(int val,int cur_state,int cur_depth);
    void SetFirstState(int state);
    void step();
    
    void stepTrain(double* mfc_step1,double* mfc_step2);    
    void reset();
    int letsplay();
    int stdletsplay();
    void letsplayreset();
    void normalizerow(int d,int r,int c);
    int iter;
    double * current_step;
    int pi,pj,pd;
    int lastnorm;
    int lastpi;
    void SetGrouping(int s_grouping);
    void MatrixPrinter();

    void StdSetStates(int set_number_states,double** mfcc_states,int s_num_mfcc_coef);
    void StdBuildChainMatrix();
    void StdstepTrain(double* mfc_step1,double* mfc_step2);

private:
    double EuclidDistance(double *a, double *b);
    DBWrite writer;
    FileWriteR f_writer;

    TNT::Array3D<double> *matrix;
    TNT::Array2D<double> *mfcc_vectors;
    double** mfcc_matrix;

    TNT::Array1D<double> *first_states;

    int first_train;
    TNT::Array2D<int> *train_matrix;
    TNT::Array2D<double> *one_train_matrix;
    TNT::Array1D<double> *cumal_distribution;
    
    bool trig;
    int number_states;
    int depth;
    int num_mfccs;
    int num_mfcc_coef;
    int cur_learn_state;
    bool first_learnt;
    int learn_pd;
    int grouping;
    
    double closest_dist_one,closest_dist_two;
    double tmp_dist_one,tmp_dist_two;
    int closest_state_one,closest_state_two;
    double normal;
    
    
    bool CheckMatrixSetupParams(int set_number_states,int set_depth);
    bool CheckTotalProbability(int state,int column);
    void PrintAnyMatrix(TNT::Array3D<double> matrix,int D,int N,int K);
    void PrintAnyMatrix(TNT::Array2D<int> s_matrix,int N,int K);

    // http://stackoverflow.com/questions/279854/how-do-i-sort-a-vector-of-pairs-based-on-the-second-element-of-the-pair
    // thanks to Evan Teran.
    // used as a custom sort type
    struct sort_pred {
        bool operator()(const std::pair<int,float> &left, const std::pair<int,float> &right) {
            return left.second < right.second;
            }
        };
            
};

#endif /* defined(__Simulsig__Markov_Chain__) */
