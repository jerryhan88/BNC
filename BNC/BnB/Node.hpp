//
//  Node.hpp
//  BNC
//
//  Created by Chung-Kyun HAN on 2/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Node_hpp
#define Node_hpp

#include <stdio.h>
#include <string>

#include "../MathematicalModel.hpp"
#include "../Problem.hpp"
#include "../Etc.hpp"

class Node {
public:
    std::string nid;
    int numIter;
    double upperBound;
    double **x_ij;
    MathematicalModel *mm;
    FilePathOrganizer *fpo;
    bool isIntegral;
    //
    Node(std::string, Problem *, FilePathOrganizer *);
    ~Node();
    //
    void calc_bound();
private:
    bool valid_IntegerSolution();
    bool gen_cuts_update_model();
};


#endif /* Node_hpp */
