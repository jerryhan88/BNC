//
//  Problem.hpp
//  BNC
//
//  Created by Chung-Kyun HAN on 1/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Problem_hpp
#define Problem_hpp

#include <stdio.h>
#include <vector>
#include <set>
#include <cmath>
#include <float.h>

#include "Etc.hpp"


class Problem {
public:
    double bv, bw, bu;
    //
    std::vector<int> K, PD, S, N;
    std::set<int> P, D;
    int o, d, *h_k, *n_k, **c_ij;
    double *r_k, *v_k, *w_k;
    double *al_i, *be_i;
    double **t_ij;
    double M;
    //
    Problem(int,
            double *, double *, double *,
            int,
            double **, double **,
            double, double, double);
    ~Problem();
};

Problem gen_problemInstance(int numTasks, int num_rrPoints, int maxReward,
                            double bv, double bw, double bu);


#endif /* Problem_hpp */

