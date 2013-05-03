//
//  HammingWindow.cpp
//  test
//
//  Created by Thomas Rushmore on 10/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#include "HammingWindow.h"

HammingWindow::HammingWindow()
{
    ham_win = NULL;
    result_vector = NULL;
}

HammingWindow::~HammingWindow()
{
    if(ham_win){
        free(ham_win);
    }
    if(result_vector){
        free(result_vector);
    }
}

double* HammingWindow::MakeWindow(const u_short_int *s_win_len)
{
    if((*s_win_len) <=0){
        printf("Hamming window must be > 0\n");
    } else {
        win_len = const_cast<u_short_int*>(s_win_len);
        ham_win = (double*)malloc(sizeof(double)*(*win_len));
        result_vector  = (double*)malloc(sizeof(double)*(*win_len));
        vDSP_hamm_windowD(ham_win, (*win_len), 0);
    }
    return ham_win;
}

double* HammingWindow::Process(double *input_vector)
{
    vDSP_vmulD(input_vector,1,ham_win,1,result_vector,1,(*win_len));
    return result_vector;
}

double* HammingWindow::GetWindow()
{
    return ham_win;
}
