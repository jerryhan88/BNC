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
    std::vector< std::string > ch_names;
//    ch_names.push_back("SE0");
//    ch_names.push_back("SE1");
//    ch_names.push_back("SE2");
    //
    std::string _header("nid,numIter,upperBound");
    for (std::string cn: ch_names) {
        _header += "," + cn;
    }
    _header += ",note";
    char header[_header.size() + 1];
    std::strcpy(header, _header.c_str());
    createCSV(this->fpo->logPath, header);
    //
    cc = new cut_composer(prob, ch_names);
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
    n->cc = cc->clone();
    return n;
}

int Node::calc_bound() {
    bool isAnyCutGenerated;
    int modelState;
    while (true) {
        (*(*mm).grbModel).optimize();
        modelState = (*(*mm).grbModel).get(GRB_IntAttr_Status);
        if (modelState == GRB_INFEASIBLE) {
            break;
        } else if (modelState == GRB_UNBOUNDED) {
            // Check the reason
            break;
        } else if (modelState != GRB_OPTIMAL) {
            (*(*mm).grbModel).computeIIS();
            (*(*mm).grbModel).write(fpo->ilpPath);
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
//    printf("numIter: %d\n", numIter);
//    for (int i: mm->prob->N) {
//        for (int j: mm->prob->N) {
//            printf("x[%d][%d]: %f\n", i, j, x_ij[i][j]);
//        }
//    }
    
    cc->solve_seperation_problem(x_ij);
    std::vector<int> numViolatedCnsts = cc->get_numViolatedCnsts();
    if (LOGGING_MODE) {
        char buf[DEFAULT_BUFFER_SIZE];
        // nid,numIter,upperBound,..."
        sprintf(buf, "%s,%d,%f",
                nid.c_str(), numIter, upperBound);
        std::string _log(buf);
        for (int nc: numViolatedCnsts) {
            _log += "," + std::to_string(nc);
        }
        char log[_log.size() + 1];
        std::strcpy(log, _log.c_str());
        (*fpo).appendLog(log);
    }
    unsigned int num_cuts = 0;
    for (int nc: numViolatedCnsts) {
        num_cuts += nc;
    }
    cc->add_cuts(mm->grbModel, mm->x_ij);

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
    delete cc;
}
