//
//  Tree.hpp
//  BNC
//
//  Created by Chung-Kyun HAN on 2/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Tree_hpp
#define Tree_hpp

#include <stdio.h>
#include <queue>
#include "../Problem.hpp"
#include "../Etc.hpp"
#include "Node.hpp"


class Tree {
public:
    FilePathOrganizer *fpo;
    Problem *prob;
    //
    std::priority_queue<Node *, std::vector<Node *>,
    auto(*) (Node *, Node *) -> bool > pq{
                                            [](Node *left, Node *right) -> bool {
                                                if ((*left).upperBound > (*right).upperBound) {
                                                    return true;
                                                } else {
                                                    return false;
                                                }
                                            }
                                        };
    std::vector<Node *> handledNodes;
    Node *incumbent = NULL;
    //
    Tree(Problem *, FilePathOrganizer *);
    ~Tree();
    //
    void run_BnB();
private:
    bool branching();
};



#endif /* Tree_hpp */


