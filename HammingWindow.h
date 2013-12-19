//
//  HammingWindow.h
//  test
//
//  Created by Thomas Rushmore on 10/02/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef __test__HammingWindow__
#define __test__HammingWindow__

#include <Accelerate/Accelerate.h>
#include <iostream>

#include "PrintfDebug.h"
#include "GlobalTypedefs.h"

class HammingWindow{
public:
    HammingWindow();
    ~HammingWindow();
    double* Process(double* input_vector);
    double* MakeWindow(const u_short_int *s_win_len);
    double* GetWindow();
private:
    double* ham_win;
    double* result_vector;
    u_short_int *win_len;
};
#endif /* defined(__test__HammingWindow__) */
