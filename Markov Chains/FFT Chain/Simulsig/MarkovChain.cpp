//
//  Markov Chain.cpp
//  Simulsig
//
//  Created by Thomas Rushmore on 27/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "MarkovChain.h"
#include <string.h>
#include <algorithm>
#include <vector>

MarkovChain::MarkovChain()
{
    iter = 0;
    pi = 0;
    pj = 0;
    pd = 0;
    lastpi = 0;
    first_train = 0;
    cur_learn_state = 0;
    num_mfcc_coef = 40;
    num_mfccs = 40;
    learn_pd = 0;
    first_learnt = false;
    matrix = NULL;
    train_matrix = NULL;
    cumal_distribution = NULL;
    mfcc_vectors = NULL;
    trig = false;
}

MarkovChain::~MarkovChain()
{
    if(matrix){
        delete(matrix);
    }
    if(train_matrix){
        delete(train_matrix);
    }
    if(cumal_distribution){
        delete(cumal_distribution);
    }
    if(mfcc_vectors){
        delete(mfcc_vectors);
    }
}

// Initialise Markov chain states with the clustered MFCC groups
void MarkovChain::SetStates(int set_number_states, double **mfcc_states, int s_num_mfcc_coef)
{
    debug_print(("\nSetting markov chain states...\n"));
    number_states = set_number_states;
    num_mfcc_coef = s_num_mfcc_coef;
    debug_print(("Num Chain States : %d\n",number_states));
    debug_print(("Grouped MFCC coef length : %d\n",num_mfcc_coef));
    mfcc_matrix = mfcc_states;
    cumal_distribution = new TNT::Array1D<double>(number_states,0.0);
    first_states = new TNT::Array1D<double>(number_states,0.0);
    mfcc_vectors = new TNT::Array2D<double>(set_number_states,num_mfcc_coef,0.0);
    for(int i = 0 ; i < number_states; i++){
        for(int j = 0 ; j < num_mfcc_coef; j++){
            (*mfcc_vectors)[i][j] = mfcc_matrix[i][j];
        }
    }
    debug_print(("Complete\n"));
}

void MarkovChain::StdBuildChainMatrix()
{
    depth = 1;
    pd = 0;
    // Originally a 3D array as was previously non-homog markov chain
    // due to misunderstanding
    matrix = new TNT::Array3D<double>(1,number_states,number_states,0.0);
    train_matrix = new TNT::Array2D<int>(number_states,number_states,1);
}

void MarkovChain::BuildChainMatrix(int set_depth)
{
    debug_print(("Depth set to %d : \n",set_depth));
    depth = set_depth;
    matrix = new TNT::Array3D<double>(depth,number_states,number_states,0.0);
    train_matrix = new TNT::Array2D<int>(depth,number_states,1);
}

void MarkovChain::stepTrain(double *mfc_step1, double *mfc_step2)
{
    // find distance between training input feature and state vec
    closest_dist_one = EuclidDistance(mfc_step1, (*mfcc_vectors)[0]);
    closest_state_one = 0;
    tmp_dist_one = 0;
    
    closest_dist_two = EuclidDistance(mfc_step2, (*mfcc_vectors)[0]);
    closest_state_two = 0;
    tmp_dist_two = 0;
    // find the closest
    for(int i = 1 ; i < number_states; i++){
        tmp_dist_one = EuclidDistance(mfc_step1, (*mfcc_vectors)[i]);
        if(tmp_dist_one < closest_dist_one){
            closest_dist_one = tmp_dist_one;
            closest_state_one = i;
        }
        tmp_dist_two = EuclidDistance(mfc_step2, (*mfcc_vectors)[i]);
        if(tmp_dist_two < closest_dist_two){
            closest_dist_two = tmp_dist_two;
            closest_state_two = i;
        }
    }    
    pi = closest_state_one;
    int ran = (rand()/(double)RAND_MAX) * 53;
    pi = ran;
    pj = closest_state_two;
    if(pd ==0){
        SetFirstState(pi);
    }
    double normal = 1.0 / (*train_matrix)[pd][pi];
    lastnorm = (*train_matrix)[pd][lastpi];
    for(int i = 0 ; i < number_states; i++){
        (*matrix)[pd][pi][i] *= (*train_matrix)[pd][pi]-1;
        (*matrix)[pd][pi][i] *= normal;
    }
    (*train_matrix)[pd][pi]++;
    (*matrix)[pd][pi][pj] += normal;
    pd++;
    lastpi = pi;

}

double MarkovChain::EuclidDistance(double *a, double *b)
{
    double running_sum = 0.0;
    for(int i = FCOEF ; i < num_mfcc_coef; i++){
        running_sum += pow((a[i] - b[i]),2);
    }
    running_sum = sqrt(running_sum);
    return running_sum;
}

void MarkovChain::SetGrouping(int s_grouping)
{
    grouping = s_grouping;
}

void MarkovChain::SetFirstState(int state)
{
    printf("First Train state : %d\n",state);
    first_train++;
    double norms = 1.0 / (double)first_train;
    for(int i = 0 ; i < number_states; i++){
        (*first_states)[i] *= first_train - 1;
        (*first_states)[i] *= norms;
    }
    (*first_states)[state] += norms;

}

void MarkovChain::StdstepTrain(double *mfc_step1, double *mfc_step2)
{
    // find distance between training input feature and state vec
    closest_dist_one = EuclidDistance(mfc_step1, (*mfcc_vectors)[0]);
    closest_state_one = 0;
    tmp_dist_one = 0;
    closest_dist_two = EuclidDistance(mfc_step2, (*mfcc_vectors)[0]);
    closest_state_two = 0;
    tmp_dist_two = 0;
    // find closest state to input feature
    for(int i = 1 ; i < number_states; i++){
        tmp_dist_one = EuclidDistance(mfc_step1, (*mfcc_vectors)[i]);
        if(tmp_dist_one < closest_dist_one){
            closest_dist_one = tmp_dist_one;
            closest_state_one = i;
        }
        tmp_dist_two = EuclidDistance(mfc_step2, (*mfcc_vectors)[i]);
        if(tmp_dist_two < closest_dist_two){
            closest_dist_two = tmp_dist_two;
            closest_state_two = i;
        }
    }
    pi = closest_state_one;
    pj = closest_state_two;
    if(!trig){
        SetFirstState(pi);
        trig = true;
    }
    normal = 1.0 / (*train_matrix)[pd][pi];
    lastnorm = (*train_matrix)[pd][lastpi];
    for(int i = 0 ; i < number_states; i++){
        (*matrix)[pd][pi][i] *= (*train_matrix)[pd][pi]-1;
        (*matrix)[pd][pi][i] *= normal;
    }
    (*train_matrix)[pd][pi]++;
    (*matrix)[pd][pi][pj] += normal;
    lastpi = pi;
}

int MarkovChain::stdletsplay()
{
    // if first play, determine the starting state
    if(!first_learnt){
        std::pair<int, double> *pair_array = new std::pair<int, double>[number_states];       
        for(int i = 0 ; i < number_states; i++){
            pair_array[i].first = i;
            pair_array[i].second = (*first_states)[i];
            
        }
        std::sort(pair_array,pair_array+number_states,sort_pred());
        std::vector<double> cdf;
        double running = 0;
        for(int i = 0 ; i < number_states; i++){
            cdf.push_back(pair_array[i].second + running);
            running += pair_array[i].second;
        }
        double ran = (rand()/(double)RAND_MAX);
        std::vector<double>::iterator low = std::upper_bound(cdf.begin(), cdf.end(), ran);
        int lb = (int)(low - cdf.begin());
        cur_learn_state = pair_array[lb].first;
        first_learnt = true;
        return pair_array[lb].first;
    }
    pd = 0;
    std::pair<int, double> *pair_array = new std::pair<int, double>[number_states];
    for(int i = 0 ; i < number_states; i++){
        pair_array[i].first = i;
        pair_array[i].second = (*matrix)[learn_pd][cur_learn_state][i];
    }
    std::sort(pair_array,pair_array+number_states,sort_pred());
    // Cumalative Distribution Function
    std::vector<double> cdf;
    double running = 0;
    for(int i = 0 ; i < number_states; i++){
        cdf.push_back(pair_array[i].second + running);
        running += pair_array[i].second;
    }
    std::vector<double>::iterator low = std::upper_bound(cdf.begin(), cdf.end(), (rand()/(double)RAND_MAX));
    int lb = (int)(low - cdf.begin());
    cur_learn_state = pair_array[lb].first ;
    delete(pair_array);
    return cur_learn_state;
}

int MarkovChain::letsplay()
{
    if(!first_learnt){
        std::pair<int, double> *pair_array = new std::pair<int, double>[number_states];
        for(int i = 0 ; i < number_states; i++){
            pair_array[i].first = i;
            pair_array[i].second = (*first_states)[i];
        }
        std::sort(pair_array,pair_array+number_states,sort_pred());
        std::vector<double> cdf;
        double running = 0;
        for(int i = 0 ; i < number_states; i++){
            cdf.push_back(pair_array[i].second + running);
            running += pair_array[i].second;
        }
        double ran = (rand()/(double)RAND_MAX);
        std::vector<double>::iterator low = std::lower_bound(cdf.begin(), cdf.end(), ran);
        int lb = (int)(low - cdf.begin());
        cur_learn_state = pair_array[lb].first;
        first_learnt = true;
        return pair_array[lb].first;
    }
    pd = 0;
    std::pair<int, double> *pair_array = new std::pair<int, double>[number_states];
    for(int i = 0 ; i < number_states; i++){
        pair_array[i].first = i;
        pair_array[i].second = (*matrix)[learn_pd][cur_learn_state][i];
    }
    std::sort(pair_array,pair_array+number_states,sort_pred());
    std::vector<double> cdf;
    double running = 0;
    for(int i = 0 ; i < number_states; i++){
        cdf.push_back(pair_array[i].second + running);
        running += pair_array[i].second;
    }
    double ran = (rand()/(double)RAND_MAX) ;
    std::vector<double>::iterator low = std::lower_bound(cdf.begin(), cdf.end(), ran);
    int lb = (int)(low - cdf.begin());
    learn_pd++;
    cur_learn_state = pair_array[lb].first ;// will equal last take picked
    delete(pair_array);
    return cur_learn_state;
}

void MarkovChain::letsplayreset()
{
    first_learnt = false;
    learn_pd = 0;
}

void MarkovChain::reset()
{
    pi = 0;
    pj = 0;
    pd = 0;
    lastpi = 0;
    trig = false;
}

void MarkovChain::BuildChainZeroed(int set_number_states, int set_depth)
{
    if(!CheckMatrixSetupParams(set_number_states, set_depth)){
        return;
    }
    matrix = new TNT::Array3D<double>(depth,number_states,number_states,0.0);
    train_matrix = new TNT::Array2D<int>(depth,number_states,1);
    cumal_distribution = new TNT::Array1D<double>(number_states,0.0);
    first_states = new TNT::Array1D<double>(number_states,0.0);
    
}

void MarkovChain::BuildChainNormalized(int set_number_states, int set_depth)
{
    if(!CheckMatrixSetupParams(set_number_states, set_depth)){
        return;
    }
    double norm_probability = 1.0/(double)number_states;
    printf("Normalized probability value = %3.10f\n",norm_probability);
    matrix = new TNT::Array3D<double>(depth,number_states,number_states,norm_probability);
    
    
}

bool MarkovChain::CheckMatrixSetupParams(int set_number_states, int set_depth)
{
    if(set_number_states <= 0){
        printf("0 states selected. Must have more than 1 state\n");
        return false;
    }
    if(set_depth < 1){
        printf("Depth cannot be set to 0\n");
        return false;
    }
    number_states = set_number_states;
    depth = set_depth;
    return true;
}

void MarkovChain::CheckTotalProbability()
{
    for(int i = 0 ; i < depth ; i++){
        for(int j = 0 ; j < number_states; j++){
            CheckTotalProbability(j, i);
        }
    }
}

bool MarkovChain::CheckTotalProbability(int state, int column)
{
    // use this to check that transitions from i to all j ( next states ) sum to 1
    double sum = 0;
    for(int i = 0 ; i < number_states; i++){
        sum += (*matrix)[column][state][i];
    }
    if(sum < 0.999999999999995781152506424405){
        printf("Total transition probability = %3.30f, must be equal to 1.0\n",sum);
        return false;
        
    }
    return true;
}

void MarkovChain::PrintMatrix()
{
    PrintAnyMatrix(*matrix, depth, number_states, number_states);
    std::cout << "Train Matrix" << std::endl;
    PrintAnyMatrix(*train_matrix,  1, number_states);
}

void MarkovChain::MatrixPrinter()
{
    for(int i = 0 ; i < depth ; i++)
    {
        for(int j = 0 ; j < number_states ; j++)
        {
            for(int l = 0 ; l < number_states ; l++)
            {
             //   printf("%2.3f ",(*matrix)[i][j][l]);

                writer.AddMarkov((*train_matrix)[j][l]);
            }
        }
    }
    f_writer.DBWriteMarkov(writer.GetMarkov());

}

void MarkovChain::PrintAnyMatrix(TNT::Array3D<double> matrix,int D,int N,int K)
{
    printf("Print Any\n");
    printf("D : %d N : %d K : %d\n",D,N,K);
    printf("\n");
    for(int i = 0 ; i < D ; i++)
    {
        for(int j = 0 ; j < N ; j++)
        {
            for(int l = 0 ; l < K ; l++)
            {
                printf("%2.3f ",matrix[i][j][l]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

void MarkovChain::PrintAnyMatrix(TNT::Array2D<int> s_matrix,int N,int K)
{
    printf("Print Any\n");
    printf("N : %d K : %d\n",N,K);
    printf("\n");

        for(int j = 0 ; j < N ; j++)
        {
            for(int l = 0 ; l < K ; l++)
            {
                printf("%d ",s_matrix[j][l]);
            }
            printf("\n");
        }
        printf("\n");
    
}

void MarkovChain::StdSetStates(int set_number_states, double **mfcc_states, int s_num_mfcc_coef)
{
    debug_print(("\nSetting standard markov chain states..."));
    // do some error checks
    number_states = set_number_states;
    debug_print(("Num states: %d \n",number_states));
    num_mfcc_coef = s_num_mfcc_coef;
    debug_print(("Num Coef : %d \n",num_mfcc_coef));
    mfcc_matrix = mfcc_states;
    cumal_distribution = new TNT::Array1D<double>(number_states,0.0);
    first_states = new TNT::Array1D<double>(number_states,0.0);
    mfcc_vectors = new TNT::Array2D<double>(set_number_states,num_mfcc_coef,0.0);
    for(int i = 0 ; i < number_states; i++){
        for(int j = 0 ; j < num_mfcc_coef; j++){
            (*mfcc_vectors)[i][j] = mfcc_matrix[i][j];
        }
    }
}

