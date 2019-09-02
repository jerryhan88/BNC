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
#include <string>
#include "../Problem.hpp"
#include "../Etc.hpp"
#include "../MathematicalModel.hpp"
#include "Node.hpp"


class Tree {
public:
    Problem *prob;
    FilePathOrganizer *fpo;
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
    MathematicalModel* get_incumbentMM();
    void log_simpleNote(Node *n, std::string note);
    void run_BnB();
private:
    bool branching();
};


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


#endif /* Tree_hpp */
