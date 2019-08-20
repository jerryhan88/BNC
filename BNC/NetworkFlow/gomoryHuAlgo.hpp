//
//  gomoryHuAlgo.hpp
//  BNC
//
//  Created by Chung-Kyun HAN on 20/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef gomoryHuAlgo_hpp
#define gomoryHuAlgo_hpp

#include <vector>
#include <array>
#include <set>
#include <map>
#include <algorithm>

#include "graph.hpp"
#include "FordFulkersonAlgo.hpp"

double** get_GomoryHuTree(double **graph, int numNodes);

void test_GomoryHuTree();

#endif /* gomoryHuAlgo_hpp */


