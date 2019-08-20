//
//  gomoryHuAlgo.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 20/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "gomoryHuAlgo.hpp"

struct NodeT;
struct EdgeT;

struct NodeT {
public:
    std::set<int> nodes;
    std::vector<EdgeT*> edges;
    int nid = -1;
    //
    NodeT() {};
    NodeT(std::vector<int> vec) {
        this->nodes.insert(vec.begin(), vec.end());
    };
    ~NodeT() {
        nodes.clear();
        edges.clear();
    };
    //
    void addEdge(EdgeT *e) {
        edges.push_back(e);
    }
    NodeT* gen_newNodeT(std::set<int> givenSet) {
        std::vector<int> intersect;
        std::set_intersection(givenSet.begin(), givenSet.end(),
                              nodes.begin(), nodes.end(),
                              std::inserter(intersect, intersect.begin()));
        return new NodeT(intersect);
    }
    void update_last() {
        assert(nodes.size() == 1);
        this->nid = *(this->nodes.begin());
    }
};

struct EdgeT {
public:
    NodeT *n0, *n1;
    double c;
    int nid0, nid1;
    //
    EdgeT(NodeT *n0, NodeT *n1, double c) {
        this->n0 = n0;
        this->n1 = n1;
        this->c = c;
        //
        this->n0->addEdge(this);
        this->n1->addEdge(this);
    }
    void update_last() {
        assert(n0->nid != -1);
        assert(n1->nid != -1);
        //
        nid0 = n0->nid;
        nid1 = n1->nid;
    }
};

double** get_GomoryHuTree(double **graph, int numNodes) {
    // Step 1
    std::vector<NodeT*> V_T;
    NodeT *vt0 = new NodeT();
    for (int i = 0; i < numNodes; i ++)
        vt0->nodes.insert(i);
    V_T.push_back(vt0);
    std::vector<EdgeT*> E_T;
    while (true) {
        // Step 2
        int X_i = -1;
        for (int i = 0; i < V_T.size(); i++) {
            if (V_T[i]->nodes.size() >= 2) {
                X_i = i;
                break;
            }
        }
        
        if (X_i == -1) {
            // Step 6
            for (NodeT *n: V_T) {
                n->update_last();
            }
            for (EdgeT *e: E_T) {
                e->update_last();
            }
            break;
        }
        NodeT *X = V_T[X_i];
        
        // Step 3
        std::vector< std::vector<int> > S;
        double **T;
        T = new double*[V_T.size()];
        std::map<NodeT *, int> mapper_NodeT;
        for (int i = 0; i < V_T.size(); i++) {
            T[i] = new double[V_T.size()];
            mapper_NodeT.insert(std::pair<NodeT *, int>(V_T[i], i));
            for (int j = 0; j < V_T.size(); j++) {
                T[i][j] = 0.0;
            }
        }
        for (EdgeT *e: E_T) {
            int n0_id = mapper_NodeT.find(e->n0)->second;
            int n1_id = mapper_NodeT.find(e->n1)->second;
            T[n0_id][n1_id] = e->c;
            T[n1_id][n0_id] = e->c;
        }
        for (EdgeT *e: X->edges) {
            std::vector<int> S_C;
            //
            int n0_id = mapper_NodeT.find(e->n0)->second;
            int n1_id = mapper_NodeT.find(e->n1)->second;
            T[n0_id][n1_id] = 0.0;
            T[n1_id][n0_id] = 0.0;
            int s = e->n0 != X ? n0_id : n1_id;
            bool visited[V_T.size()];
            std::memset(visited, false, sizeof(visited));
            reachability_dfs(T, (int) V_T.size(), s, visited);
            for (int i = 0; i < V_T.size(); i++) {
                if (visited[i]) {
                    for (int nid: V_T[i]->nodes) {
                        S_C.push_back(nid);
                    }
                }
            }
            S.push_back(S_C);
        }
        for (int i = 0; i < V_T.size(); i++)
            delete[] T[i];
        delete T;
        //
        std::vector< std::vector<int> > V_G_;
        std::vector<int> V_G_membership;
        std::vector< std::array<int, 2> > E_G_;
        std::vector<double> c_;
        for (int nid : X->nodes) {
            V_G_.push_back(std::vector<int>(1, nid));
            V_G_membership.push_back((int) S.size());
        }
        for (int i = 0; i < S.size(); i++) {
            V_G_.push_back(S[i]);
            V_G_membership.push_back(i);
        }
        
        for (int i = 0; i < V_G_.size(); i++) {
            for (int j = 0; j < V_G_.size(); j++) {
                if (i == j) {
                    continue;
                }
                double c = 0.0;
                if (V_G_membership[i] == S.size() && V_G_membership[j] == S.size()) {
                    assert(V_G_[i].size() == 1);
                    assert(V_G_[j].size() == 1);
                    //
                    int u = V_G_[i][0];
                    int v = V_G_[j][0];
                    c += graph[u][v];
                } else if (V_G_membership[i] == S.size() && V_G_membership[j] != S.size()) {
                    assert(V_G_[i].size() == 1);
                    int u = V_G_[i][0];
                    for (int v: V_G_[j]) {
                        c += graph[u][v];
                    }
                } else if (V_G_membership[i] != S.size() && V_G_membership[j] == S.size()) {
                    assert(V_G_[j].size() == 1);
                    int v = V_G_[j][0];
                    for (int u: V_G_[i]) {
                        c += graph[u][v];
                    }
                } else if (V_G_membership[i] != S.size() && V_G_membership[j] != S.size()) {
                    for (int u: V_G_[i]) {
                        for (int v: V_G_[j]) {
                            c += graph[u][v];
                        }
                    }
                } else {
                    assert(false);
                }
                if (c == 0.0)
                    continue;
                E_G_.push_back(std::array<int, 2> {i, j});
                c_.push_back(c);
            }
        }
        // Step 4
        
        double **G_ = gen_graphInstance((int) V_G_.size(), E_G_, c_);
        int s = -1, t = -1;
        int _s, _t;
        std::set<int>::iterator it = X->nodes.begin();
        _s = *it;
        it++;
        _t = *it;
        for (int i = 0; i < V_G_.size(); i++) {
            if (V_G_[i].size() != 1)
                continue;
            if(V_G_[i][0] == _s)
                s = i;
            if(V_G_[i][0] == _t)
                t = i;
            if (s != -1 && t != -1)
                break;
        }
        std::vector<std::array<int, 2>> cut = get_minCut(G_, (int) V_G_.size(), s, t);
        double minCapa = 0.0;
        for (auto e : cut) {
            minCapa += G_[e[0]][e[1]];
            G_[e[0]][e[1]] = 0.0;
            G_[e[1]][e[0]] = 0.0;
        }
        bool visited[V_G_.size()];
        std::memset(visited, false, sizeof(visited));
        reachability_dfs(G_, (int) V_G_.size(), s, visited);
        std::set<int> A, B;
        
        for (int i = 0; i < V_G_.size(); i++) {
            if (visited[i]) {
                A.insert(V_G_[i].begin(), V_G_[i].end());
            } else {
                B.insert(V_G_[i].begin(), V_G_[i].end());
            }
        }
        
        // Step 5
        NodeT *A_X = X->gen_newNodeT(A);
        NodeT *B_X = X->gen_newNodeT(B);
        for (NodeT *vt: {A_X, B_X}) {
            V_T.push_back(vt);
        }
        for (EdgeT *e: E_T) {
            if (e->n0 == X) {
                std::vector<int> difference;
                std::set_difference(e->n1->nodes.begin(), e->n1->nodes.end(),
                                    A.begin(), A.end(),
                                    std::inserter(difference, difference.begin()));
                if (difference.size() == 0) {
                    e->n0 = A_X;
                    A_X->addEdge(e);
                } else {
                    e->n0 = B_X;
                    B_X->addEdge(e);
                }
            } else if (e->n1 == X) {
                std::vector<int> difference;
                std::set_difference(e->n0->nodes.begin(), e->n0->nodes.end(),
                                    A.begin(), A.end(),
                                    std::inserter(difference, difference.begin()));
                if (difference.size() == 0) {
                    e->n1 = A_X;
                    A_X->addEdge(e);
                } else {
                    e->n1 = B_X;
                    B_X->addEdge(e);
                }
            }
        }
        
        E_T.push_back(new EdgeT(A_X, B_X, minCapa));
        V_T.erase(V_T.begin() + X_i);
        delete X;
    }
    
    double **GomoryHuTree = gen_graphInstance(numNodes);
    for (EdgeT *e: E_T) {
        GomoryHuTree[e->nid0][e->nid1] = e->c;
        delete e;
    }
    E_T.clear();
    for (NodeT *n: V_T) {
        delete n;
    }
    V_T.clear();
    
    return GomoryHuTree;
}


void test_GomoryHuTree() {

    //    const int V = 6;
    //    double graph[V][V] = {
    //        {0.0, 1.0, 7.0, 0.0, 0.0, 0.0},
    //        {1.0, 0.0, 1.0, 3.0, 2.0, 0.0},
    //        {7.0, 1.0, 0.0, 0.0, 4.0, 0.0},
    //        {0.0, 3.0, 0.0, 0.0, 1.0, 6.0},
    //        {0.0, 2.0, 4.0, 1.0, 0.0, 2.0},
    //        {0.0, 0.0, 0.0, 6.0, 2.0, 0.0}
    //    };
    
    const int V = 6;
    double graph[V][V] = {
        {0.0,  10.0, 0.0, 0.0, 0.0, 8.0},
        {10.0,  0.0, 4.0, 0.0, 2.0, 3.0},
        {0.0,   4.0, 0.0, 5.0, 4.0, 2.0},
        {0.0,   0.0, 5.0, 0.0, 7.0, 2.0},
        {0.0,   2.0, 4.0, 7.0, 0.0, 3.0},
        {8.0,   3.0, 2.0, 2.0, 3.0, 0.0}
    };
    
    double **G = gen_graphInstance(graph);
    
    double **GomoryHuTree = get_GomoryHuTree(G, V);
    
    for (int i = 0; i < V; i++)
        delete [] GomoryHuTree[i];
    delete GomoryHuTree;
}
