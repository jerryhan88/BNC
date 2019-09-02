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
#include "../Cut.hpp"
#include "../Etc.hpp"
#include "../NetworkFlow/gomoryHuAlgo.hpp"

#define INFEASIBLE_MODEL GRB_INFEASIBLE

class Node {
public:
    std::string nid;
    int numIter;
    double upperBound = -1.0;
    double **x_ij;
    MathematicalModel *mm;
    FilePathOrganizer *fpo;
    bool isIntegral;
    int mf_i = -1, mf_j = -1;  // Index of the most frational x_ij
    //
    cut_composer *cc;
    //
    Node(std::string, Problem *, FilePathOrganizer *);
    Node(){};
    ~Node();
    //
    int calc_bound();
    void search_mostFractionalVarIndex();
    Node* clone(std::string);
private:
    bool valid_IntegerSolution();
    bool gen_cuts_update_model();
};


#endif /* Node_hpp */
