//
//  Node.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 2/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Node.hpp"

#define LOGGING_MODE true
#define DEFAULT_BUFFER_SIZE 2048

Node::Node(std::string nid, Problem *prob, FilePathOrganizer *fpo) {
    this->nid = nid;
    this->fpo = fpo;
    numIter = 0;
    mm = new MathematicalModel('L', prob);
    x_ij = new double *[(*prob).N.size()];
    for (int i: (*prob).N) {
        x_ij[i] = new double[(*prob).N.size()];
    }
}

void Node::calc_bound() {
    bool isAnyCutGenerated;
    while (true) {
        (*(*mm).grbModel).optimize();
        int modelState = (*(*mm).grbModel).get(GRB_IntAttr_Status);
        if (modelState == GRB_INFEASIBLE) {
            break;
        } else if (modelState != GRB_OPTIMAL) {
            (*(*mm).grbModel).computeIIS();
            (*(*mm).grbModel).write("/Users/ckhan/workspace/BNC/BNC/test.ilp");
//    (*mm.grbModel).write("/Users/ckhan/workspace/BNC/BNC/test_LP_cut.lp");
        }
        (*mm).get_x_ij(x_ij);
        isIntegral = valid_IntegerSolution();
        if (isIntegral)
            break;
        upperBound = (*(*mm).grbModel).get(GRB_DoubleAttr_ObjVal);
        isAnyCutGenerated = gen_cuts_update_model();
        if (!isAnyCutGenerated)
            break;
        numIter += 1;
    }
}

bool Node::valid_IntegerSolution() {
    bool isAllIntegers = true;
    for (int i: (*(*mm).prob).N) {
        for (int j: (*(*mm).prob).N) {
            if (x_ij[i][j] == 0.0 || x_ij[i][j] == 1.0) {
                continue;
            } else {
                isAllIntegers = false;
                break;
            }
        }
        if (!isAllIntegers)
            break;
    }
    return isAllIntegers;
}

bool Node::gen_cuts_update_model() {
    std::vector<std::set<int>> SEC1 = get_SEC(1, x_ij, (*(*mm).prob));
    std::vector<std::set<int>> SEC2 = get_SEC(2, x_ij, (*(*mm).prob));
    if (LOGGING_MODE) {
        char buf[DEFAULT_BUFFER_SIZE];
        // nid,numIter,upperBound,SEC1,SEC2"
        sprintf(buf, "%s,%d,%f,%d,%d",
                nid.c_str(), numIter, upperBound,
                (int) SEC1.size(), (int) SEC2.size());
        (*fpo).appendLog(buf);
    }
    int num_cuts = 0;
    (*mm).add_SEC(SEC1, 1);
    num_cuts += (int) SEC1.size();
    (*mm).add_SEC(SEC2, 2);
    num_cuts += (int) SEC2.size();
    //
    return num_cuts > 0 ? true : false;
}

Node::~Node() {
    for (int i = 0; i < (*(*mm).prob).N.size(); i++)
        delete [] x_ij[i];
    delete [] x_ij;
    delete mm;
}
