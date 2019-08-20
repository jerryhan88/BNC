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
#include "BnB/Tree.hpp"
#include "MathematicalModel.hpp"

#define BUF_SIZE 1024


Problem ex1() {
    srand(1);
//    int numTasks = 6, num_rrPoints = 2, maxReward = 3;
    int numTasks = 3, num_rrPoints = 2, maxReward = 3;
    double bv = 2.5, bw = 3.0, bu = 1.0 * 3;
    
    Problem pi = gen_problemInstance(numTasks, num_rrPoints, maxReward,
                                     bv, bw, bu);
    return pi;
}


int main(int argc, const char * argv[]) {
 
    
    FilePathOrganizer fpo(argv);
    Problem prob = ex1();
    
    MathematicalModel mm('I', &prob);
    (*mm.grbModel).optimize();
    
    Tree tree(&prob, &fpo);
    tree.run_BnB();

    return 0;
}



/*
 Priority queue example
 
 template<typename T> void print_queue(T& q) {
 while(!q.empty()) {
 std::cout << q.top() << " ";
 q.pop();
 }
 std::cout << '\n';
 }
 
 std::priority_queue<int> q;
 
 for(int n : {1,8,5,6,3,4,0,9,7,2})
 q.push(n);
 
 print_queue(q);
 
 std::priority_queue<int, std::vector<int>, std::greater<int> > q2;
 
 for(int n : {1,8,5,6,3,4,0,9,7,2})
 q2.push(n);
 
 print_queue(q2);
 
 // Using lambda to compare elements.
 auto cmp = [](int left, int right) { return left > right;};
 std::priority_queue<int, std::vector<int>, decltype(cmp)> q3(cmp);
 
 for(int n : {1,8,5,6,3,4,0,9,7,2})
 q3.push(n);
 
 print_queue(q3);
 */
