//
//  main.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 20/7/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

//
//#include "PDPTW_ILP.hpp"
#include <iostream>
#include "Problem.hpp"
#include "Etc.hpp"
#include "BnB/Node.hpp"

#define BUF_SIZE 1024


Problem ex1() {
    srand(1);
    int numTasks = 6, num_rrPoints = 2, maxReward = 3;
    double bv = 2.5, bw = 3.0, bu = 1.0 * 3;
    
    Problem pi = gen_problemInstance(numTasks, num_rrPoints, maxReward,
                                     bv, bw, bu);
    return pi;
}



int main(int argc, const char * argv[]) {

    FilePathOrganizer fpo(argv);
    Problem prob = ex1();
    
    Node n0("n0", &prob, &fpo);
    n0.calc_bound();
    
    return 0;
}
