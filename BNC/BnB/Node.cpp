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
    se1_cut = new SE1_cut("SE1", prob);
    se2_cut = new SE2_cut("SE2", prob);
}


Node* Node::clone(std::string new_nid) {
    Node *n = new Node();
    n->nid = new_nid;
    n->fpo = fpo;
    n->numIter = 0;
    n->mm = mm->clone();
    n->x_ij = new double *[(*(*mm).prob).N.size()];
    for (int i: (*(*mm).prob).N) {
        n->x_ij[i] = new double[(*(*mm).prob).N.size()];
    }
    n->ES_SEC1.insert(ES_SEC1.begin(), ES_SEC1.end());
    n->ES_SEC2.insert(ES_SEC2.begin(), ES_SEC2.end());
    return n;
}

int Node::calc_bound() {
    bool isAnyCutGenerated;
    int modelState;
    while (true) {
        (*(*mm).grbModel).optimize();
        modelState = (*(*mm).grbModel).get(GRB_IntAttr_Status);
        if (modelState == GRB_INFEASIBLE) {
            (*(*mm).grbModel).computeIIS();
            (*(*mm).grbModel).write("/Users/ckhan/workspace/BNC/BNC/test.ilp");
            break;
        } else if (modelState != GRB_OPTIMAL) {
            (*(*mm).grbModel).computeIIS();
            (*(*mm).grbModel).write("/Users/ckhan/workspace/BNC/BNC/test.ilp");
        }
        (*mm).get_x_ij(x_ij);
        upperBound = (*(*mm).grbModel).get(GRB_DoubleAttr_ObjVal);
        isIntegral = valid_IntegerSolution();
        if (isIntegral)
            break;
        isAnyCutGenerated = gen_cuts_update_model();
        if (!isAnyCutGenerated)
            break;
        numIter += 1;
    }
    return modelState;
}

void Node::search_mostFractionalVarIndex() {
    float mostFractionalVal = 0.5, fractionalVal;
    for (int i: (*(*mm).prob).N) {
        for (int j: (*(*mm).prob).N) {
            fractionalVal = abs(x_ij[i][j] - 0.5);
            if (fractionalVal < mostFractionalVal) {
                mostFractionalVal = fractionalVal;
                mf_i = i;
                mf_j = j;
            }
        }
    }
}

bool Node::valid_IntegerSolution() {
    bool isAllIntegers = true;
    float fractionalVal;
    for (int i: (*(*mm).prob).N) {
        for (int j: (*(*mm).prob).N) {
            fractionalVal = abs(x_ij[i][j] - 0.5);;
            if (fractionalVal == 0.5) {
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
    std::set<std::set<int>> SEC1 = get_SEC(1,
                                              ES_SEC1,
                                              x_ij, (*(*mm).prob));
    std::set<std::set<int>> SEC2 = get_SEC(2,
                                              ES_SEC2,
                                              x_ij, (*(*mm).prob));
    if (LOGGING_MODE) {
        char buf[DEFAULT_BUFFER_SIZE];
        // nid,numIter,upperBound,SEC1,SEC2"
        sprintf(buf, "%s,%d,%f,%d,%d",
                nid.c_str(), numIter, upperBound,
                (int) SEC1.size(), (int) SEC2.size());
        (*fpo).appendLog(buf);
    }
    int num_cuts = 0;
    (*mm).add_SEC(SEC1, 1, nid, numIter);
    ES_SEC1.insert(SEC1.begin(), SEC1.end());
    num_cuts += (int) SEC1.size();
    (*mm).add_SEC(SEC2, 2, nid, numIter);
    ES_SEC2.insert(SEC2.begin(), SEC2.end());
    num_cuts += (int) SEC2.size();
    //
    return num_cuts > 0 ? true : false;
}

Node::~Node() {
    for (int i = 0; i < (*(*mm).prob).N.size(); i++)
        delete [] x_ij[i];
    delete [] x_ij;
    if (nid == "*") {
        delete mm->env;
    }
    delete mm;
}
