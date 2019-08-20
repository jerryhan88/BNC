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
    char buf[DEFAULT_BUFFER_SIZE];
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
            sprintf(buf, "%s,-,%f,-,-,firstIncumbent",
                    (*incumbent).nid.c_str(), (*incumbent).upperBound);
            (*fpo).appendLog(buf);
        } else {
            if ( (*incumbent).upperBound < (*n0).upperBound ) {
                delete incumbent;
                incumbent = n0;
                sprintf(buf, "%s,-,%f,-,-,updatedIncumbent",
                        (*incumbent).nid.c_str(), (*incumbent).upperBound);
                (*fpo).appendLog(buf);
            }
        }
    } else {
//        if ((*n0).nid == "*1101001111000001")
//            printf("TEST");
        (*n0).search_mostFractionalVarIndex();
        sprintf(buf, "%s,-,-,-,-,branch(%d-%d)",
                (*n0).nid.c_str(), n0->mf_i, n0->mf_j);
        (*fpo).appendLog(buf);

        for (int rhs = 0; rhs < 2; rhs++) {
            Node *n1 = (*n0).clone(n0->nid + std::to_string(rhs));
            n1->mm->add_intConstr(n0->mf_i, n0->mf_j, rhs);
            modelState = (*n1).calc_bound();
            if (modelState != INFEASIBLE_MODEL) {
                if (incumbent != NULL && (*n1).upperBound <= (*incumbent).upperBound) {
                    sprintf(buf, "%s,-,-,-,-,pruned",
                             (*n1).nid.c_str());
                    (*fpo).appendLog(buf);
                    return false;
                } else {
                    pq.push(n1);
                    sprintf(buf, "%s,-,-,-,-,pushedPQ",
                            (*n1).nid.c_str());
                    (*fpo).appendLog(buf);
                }
            } else {
                sprintf(buf, "%s,-,-,-,-,infeasible",
                        (*n1).nid.c_str());
                (*fpo).appendLog(buf);
            }
        }
        if ((*n0).nid != "*")
            delete n0;
        printf("test");
    }
    return false;
}

Tree::~Tree() {
    
}
