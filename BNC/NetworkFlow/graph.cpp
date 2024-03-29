//
//  graph.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 20/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "graph.hpp"


double** gen_graphInstance(int numNodes) {
    double **G;
    G = new double*[numNodes];
    for (int i = 0; i < numNodes; i++) {
        G[i] = new double[numNodes];
        for (int j = 0; j < numNodes; j++) {
            G[i][j] = 0.0;
        }
    }
    return G;
}

double** gen_graphInstance(double **graph, int numNodes) {
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

double** gen_graphInstance(int numNodes, std::vector< std::array<int, 2> > &edges, std::vector<double> &weights) {
    assert(edges.size() == weights.size());
    double **G = gen_graphInstance(numNodes);
    for (int i = 0; i < edges.size(); i++) {
        G[edges[i][0]][edges[i][1]] = weights[i];
    }
    
    return G;
}


bool reachability_bfs(double **rGraph, int V, int s, int t, int parent[]) {
    // Create a visited array and mark all vertices as not visited
    bool visited[V];
    std::memset(visited, 0, sizeof(visited));
    
    // Create a queue, enqueue source vertex and mark source vertex
    // as visited
    std::queue <int> q;
    q.push(s);
    visited[s] = true;
    parent[s] = -1;
    
    // Standard BFS Loop
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        
        for (int v = 0; v < V; v++) {
            if (visited[v] == false && rGraph[u][v] > 0.0) {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }
    
    // If we reached sink in BFS starting from source, then return
    // true, else false
    return (visited[t] == true);
}

void reachability_dfs(double **graph, int V, int s, bool visited[]) {
    visited[s] = true;
    for (int i = 0; i < V; i++) {
        if (graph[s][i] && !visited[i])
            reachability_dfs(graph, V, i, visited);
    }
}
