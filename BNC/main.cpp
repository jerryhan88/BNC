//
//  main.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 20/7/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include <iostream>
#include "gurobi_c++.h"
using namespace std;

int main(int argc, const char * argv[]) {
    int n = 10;
    double *x = new double[n];
    double *y = new double[n];
    
    int i;
    for (i = 0; i < n; i++) {
        x[i] = ((double) rand()) / RAND_MAX;
        y[i] = ((double) rand()) / RAND_MAX;
    }
    
    GRBEnv *env = NULL;
    GRBVar **vars = NULL;
    
    vars = new GRBVar*[n];
    

    
    return 0;
}
