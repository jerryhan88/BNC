//
//  Tree.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 2/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Tree.hpp"
#define DEFAULT_BUFFER_SIZE 1024


Tree::Tree(Problem *prob, FilePathOrganizer *fpo) {
    this->prob = prob;
    this->fpo = fpo;
}

Tree::~Tree() {
}

MathematicalModel* Tree::get_incumbentMM() {
    assert(incumbent != NULL);
    return incumbent->mm;
}

void Tree::run_BnB() {
    Node *root = new Node("*", prob, fpo);
    (*root).calc_bound();
    pq.push(root);


    bool isSearchFinished = false;
    while (!isSearchFinished) {
        isSearchFinished = branching();
    }
}


bool Tree::branching() {
    Node *n0 = NULL;
    if (!pq.empty()) {
        n0 = pq.top();
        pq.pop();
    } else {
        return true;
    }
    //
    int modelState;
    if ((*n0).isIntegral) {
        if (incumbent == NULL) {
            incumbent = n0;
            log_simpleNote(incumbent, "firstIncumbent");
        } else {
            if ( (*incumbent).upperBound < (*n0).upperBound ) {
                delete incumbent;
                incumbent = n0;
                log_simpleNote(incumbent, "updatedIncumbent");
            }
        }
    } else {
//        if ((*n0).nid == "*10")
//            printf("TEST");
        (*n0).search_mostFractionalVarIndex();
        log_simpleNote(n0, "branch(" + std::to_string(n0->mf_i) + "-" + std::to_string(n0->mf_j) + ")");

        for (int rhs = 0; rhs < 2; rhs++) {
            Node *n1 = (*n0).clone(n0->nid + std::to_string(rhs));
//            if ((*n1).nid == "*11")
//                printf("TEST");
            n1->mm->add_intConstr(n0->mf_i, n0->mf_j, rhs);
            modelState = (*n1).calc_bound();
//            if (modelState == INFEASIBLE_MODEL) {
//                n1->mm->grbModel->computeIIS();
//                n1->mm->grbModel->write("/Users/ckhan/workspace/BNC/BNC/temp.ilp");
//            }
            
            if (modelState != INFEASIBLE_MODEL) {
                if (incumbent != NULL && (*n1).upperBound <= (*incumbent).upperBound) {
                    log_simpleNote(n1, "pruned");
                    return false;
                } else {
                    pq.push(n1);
                    log_simpleNote(n1, "pushedPQ");
                }
            } else {
                log_simpleNote(n1, "infeasible");
            }
        }
        if ((*n0).nid != "*")
            delete n0;
    }
    return false;
}

void Tree::log_simpleNote(Node *n, std::string note) {
    std::string _row;
    _row += (*n).nid;
    _row += ",-";
    if ((*n).upperBound != -1) {
        _row += "," + std::to_string((*n).upperBound);
    } else {
        _row += ",-";
    }
    for (int i = 0; i < n->cc->chs.size(); i++)
        _row += ",-";
    _row += "," + note;
    char buf[_row. size() + 1];
    std::strcpy(buf, _row.c_str());
    appendRow(fpo->logPath, buf);
}
