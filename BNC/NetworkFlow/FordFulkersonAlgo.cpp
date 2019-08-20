//
//  FordFulkersonAlgo.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 16/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "FordFulkersonAlgo.hpp"


void sendMaxFlow(double **rGraph, int numNodes, int source_id, int terminal_id, int parent[],
                 double *max_flow) {
    int u, v;
    while (reachability_bfs(rGraph, numNodes, source_id, terminal_id, parent)) {
        // Find minimum residual capacity of the edges along the
        // path filled by BFS. Or we can say find the maximum flow
        // through the path found.
        double path_flow = DBL_MAX;
        for (v = terminal_id; v != source_id; v = parent[v]) {
            u = parent[v];
            if (rGraph[u][v] < path_flow) {
                path_flow = rGraph[u][v];
            }
        }
        
        // update residual capacities of the edges and reverse edges
        // along the path
        for (v = terminal_id; v != source_id; v = parent[v]) {
            u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }
        
        // Add path flow to overall flow
        *max_flow += path_flow;
    }
}

double fordFulkerson(double **graph, int numNodes, int source_id, int terminal_id) {
    
    // Create a residual graph and fill the residual graph with
    // given capacities in the original graph as residual capacities
    // in residual graph
    double **rGraph = gen_graphInstance(graph, numNodes);
    
    int parent[numNodes]; // This array is filled by BFS and to store path
    double max_flow = 0.0; // There is no flow initially
    sendMaxFlow(rGraph, numNodes, source_id, terminal_id, parent,
                &max_flow);
    //
    for (int i = 0; i < numNodes; i++)
        delete [] rGraph[i];
    delete rGraph;
    //
    return max_flow;
}

std::vector<std::array<int, 2>> get_minCut(double **graph, int numNodes, int source_id, int terminal_id) {
    std::vector<std::array<int, 2>> cut;
    
    double **rGraph = gen_graphInstance(graph, numNodes);
    
    int parent[numNodes];
    double max_flow = 0.0;
    sendMaxFlow(rGraph, numNodes, source_id, terminal_id, parent, &max_flow);
    
    bool visited[numNodes];
    std::memset(visited, false, sizeof(visited));
    reachability_dfs(rGraph, numNodes, source_id, visited);
    
    // Print all edges that are from a reachable vertex to
    // non-reachable vertex in the original graph
    for (int u = 0; u < numNodes; u++) {
        for (int v = 0; v < numNodes; v++) {
            if (visited[u] && !visited[v] && graph[u][v]) {
                cut.push_back(std::array<int, 2> { {u, v} });
            }
        }
    }
    
    for (int i = 0; i < numNodes; i++) {
        delete [] rGraph[i];
    }
    delete[] rGraph;
    
    return cut;
}

void test_FordFulkersonAlgo() {
    const int V = 6;
    double graph[V][V] = {
        {0.0, 16.0, 13.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 10.0, 12.0, 0.0, 0.0},
        {0.0, 4.0, 0.0, 0.0, 14.0, 0.0},
        {0.0, 0.0, 9.0, 0.0, 0.0, 20.0},
        {0.0, 0.0, 0.0, 7.0, 0.0, 4.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };

    double **G = gen_graphInstance(graph);
    std::cout << "The maximum possible flow is " << fordFulkerson(G, V, 0, 5);

    std::vector<std::array<int, 2>> cut = get_minCut(G, V, 0, 5);
}
