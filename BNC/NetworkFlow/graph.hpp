//
//  graph.hpp
//  BNC
//
//  Created by Chung-Kyun HAN on 20/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef graph_hpp
#define graph_hpp

#include <iostream>
#include <queue>
#include <array>
//



template <size_t numNodes>
double** gen_graphInstance(double (&graph)[numNodes][numNodes]) {
    double **G;
    G = new double*[numNodes];
    for (int i = 0; i < numNodes; i++) {
        G[i] = new double[numNodes];
        for (int j = 0; j < numNodes; j++) {
            G[i][j] = graph[i][j];
        }
    }
    return G;
}
double** gen_graphInstance(int numNodes);
double** gen_graphInstance(double **graph, int numNodes);
double** gen_graphInstance(int numNodes,
                           std::vector< std::array<int, 2> > &edges, std::vector<double> &weights);
//

bool reachability_bfs(double **rGraph, int V, int s, int t, int parent[]);
void reachability_dfs(double **graph, int V, int s, bool visited[]);

#endif /* graph_hpp */
