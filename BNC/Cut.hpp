//
//  Cut.hpp
//  BNC
//
//  Created by Chung-Kyun HAN on 1/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Cut_hpp
#define Cut_hpp

#include <stdio.h>
#include <set>

#include "Problem.hpp"

typedef std::pair<int, int> arc;

std::set<arc> get_delta_S(std::set<int> S, const Problem &prmt);
std::set<int> get_sigma_S(std::set<int> S, const Problem &prmt);
std::set<int> get_pi_S(std::set<int> S, const Problem &prmt);
std::set<int> get_bar_S(std::set<int> S, const Problem &prmt);


template<typename T1, typename T2>
void get_SEC1_LHS(T1 *lhs, T2 **x_ij,
                  const std::set<int> &S,
                  const Problem &prmt) {
    std::set<arc> delta_S = get_delta_S(S, prmt);
    std::set<int> bar_S = get_bar_S(S, prmt);
    std::set<int> sigma_S = get_sigma_S(S, prmt);
    /*
     * \sum_{ i, j \in \delta(S)} x_{i, j}
     */
    for (arc a: delta_S) {
        *lhs += x_ij[a.first][a.second];
    }
    /*
     * -2 * \sum_{i \in \bar{S} \cap \sigma(S)} \sum_{j \in S} x_{i, j}
     */
    for (int i: bar_S) {
        if (sigma_S.find(i) == sigma_S.end())
            continue;
        for (int j: S) {
            *lhs -= 2 * x_ij[i][j];
        }
    }
    /*
     * -2 * \sum_{i \in \bar{S} \setminus \sigma(S)} \sum_{j \in S \cap \sigma(S)} x_{i, j}
     */
    for (int i: bar_S) {
        if (sigma_S.find(i) != sigma_S.end())
            continue;
        for (int j: S) {
            if (sigma_S.find(j) == sigma_S.end())
                continue;
            *lhs -= 2 * x_ij[i][j];
        }
    }
}


template<typename T1, typename T2>
void get_SEC2_LHS(T1 *lhs, T2 **x_ij,
                  const std::set<int> &S,
                  const Problem &prmt) {
    std::set<arc> delta_S = get_delta_S(S, prmt);
    std::set<int> bar_S = get_bar_S(S, prmt);
    std::set<int> pi_S = get_pi_S(S, prmt);
    /*
     * \sum_{ i, j \in \delta(S)} x_{i, j}
     */
    for (arc a: delta_S) {
        *lhs += x_ij[a.first][a.second];
    }
    /*
     * - 2 * \sum_{i \in S} \sum_{j \in \bar{S} \cap \pi(S)} x_{i, j}
     */
    for (int i: S) {
        for (int j: bar_S) {
            if (pi_S.find(j) == pi_S.end())
                continue;
            *lhs -= 2 * x_ij[i][j];
        }
    }
    /*
     * - 2 * \sum_{i \in S \cap \pi(S)} \sum_{j \in \bar{S} \setminus \pi(S)} x_{i, j}
     */
    for (int i: S) {
        if (pi_S.find(i) == pi_S.end())
            continue;
        for (int j: bar_S) {
            if (pi_S.find(j) != pi_S.end())
                continue;
            *lhs -= 2 * x_ij[i][j];
        }
    }
}


std::vector<std::set<int>> get_SEC(int SEC_no, double **_x_ij, const Problem &prmt);



#endif /* Cut_hpp */
