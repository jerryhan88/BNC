//
//  MathematicalModel.hpp
//  BNC
//
//  Created by Chung-Kyun HAN on 1/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef MathematicalModel_hpp
#define MathematicalModel_hpp

#include <stdio.h>

#include "gurobi_c++.h"
#include "Problem.hpp"
#include "Cut.hpp"

class MathematicalModel{
public:
    const Problem *prob;
    GRBEnv *env;
    GRBVar **x_ij, *u_i;
    GRBModel *grbModel;
    //
    MathematicalModel(char, Problem *);
    MathematicalModel(){};
    ~MathematicalModel();
    //
    void add_SEC(std::set<std::set<int>> SEC, int SEC_no,
                 std::string nid, int numIter);
    void get_x_ij(double **);
    MathematicalModel* clone();
    void add_intConstr(int i, int j, int rhs);
private:
    void def_dvs(char);
    void def_FC_cnsts();
    void def_AT_cnsts();
    void def_CP_cnsts();
    void def_objF();
};


#endif /* MathematicalModel_hpp */
