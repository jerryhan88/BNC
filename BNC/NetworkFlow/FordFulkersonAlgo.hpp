//
//  FordFulkersonAlgo.hpp
//  BNC
//
//  Created by Chung-Kyun HAN on 16/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef FordFulkersonAlgo_hpp
#define FordFulkersonAlgo_hpp

#include <stdio.h>
#include <float.h>
#include <queue>
#include <vector>
#include <set>
#include <array>
#include <map>
#include <algorithm>
#include <iostream>
//
#include "graph.hpp"

std::vector<std::array<int, 2>> get_minCut(double **graph, int V, int s, int t);

void test_GomoryHuTree();

#endif /* FordFulkersonAlgo_hpp */
