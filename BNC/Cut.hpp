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
#include "gurobi_c++.h"

typedef std::pair<int, int> edge;

std::set<edge> get_delta_S(const std::set<int> &S, const Problem &prmt);
std::set<int> get_sigma_S(const std::set<int> &S, const Problem &prmt);
std::set<int> get_pi_S(const std::set<int> &S, const Problem &prmt);
std::set<int> get_bar_S(const std::set<int> &S, const Problem &prmt);


template<typename T1, typename T2>
void get_SEC1_LHS(T1 *lhs, T2 **x_ij,
                  const std::set<int> &S,
                  const Problem &prmt) {
    std::set<edge> delta_S = get_delta_S(S, prmt);
    std::set<int> bar_S = get_bar_S(S, prmt);
    std::set<int> sigma_S = get_sigma_S(S, prmt);
    /*
     * \sum_{ i, j \in \delta(S)} x_{i, j}
     */
    for (edge a: delta_S) {
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
    std::set<edge> delta_S = get_delta_S(S, prmt);
    std::set<int> bar_S = get_bar_S(S, prmt);
    std::set<int> pi_S = get_pi_S(S, prmt);
    /*
     * \sum_{ i, j \in \delta(S)} x_{i, j}
     */
    for (edge a: delta_S) {
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


std::set<std::set<int>> get_SEC(int SEC_no,
                                   const std::set<std::set<int>> &ES_SEC,
                                   double **_x_ij, const Problem &prmt);



class cut_handler {
public:
    std::string ch_name;
    Problem *prob;
    //
    cut_handler(std::string ch_name, Problem *prob){
        this->ch_name = ch_name;
        this->prob = prob;
    }
    //
    void indentify_violatedConstrains(double **x_ij) {
        throw "Should override indentify_violatedConstrains()";
    }
    void generate_cuts() {
        throw "Should override generate_cuts()";
    }
};


class SubtourEliminationCut_heuristic : public cut_handler {
public:
    std::set<std::set<int>> existing_sets;
    std::set<std::set<int>> new_sets;
    //
    SubtourEliminationCut_heuristic(std::string ch_name, Problem *prob): cut_handler(ch_name, prob) {}
    //
    void indentify_violatedConstrains(double **x_ij);
private:
    void find_new_nid_wMinLHS_wAddition(double **x_ij, const std::set<int> &S0,
                                        std::set<std::set<int>> &tabu_S,
                                        std::set<int> &tabu_N,
                                        int *min_nid, double *min_lhs);
    void find_included_nid_wMinLHS_wDeletion(double **x_ij, const std::set<int> &S0,
                                             std::set<std::set<int>> &tabu_S,
                                             std::set<int> &tabu_N,
                                             int *min_nid, double *min_lhs);
protected:
    double get_LHS_givenS(double **x_ij, const std::set<int> &S) {
        throw "Should override get_LHS_givenS()";
    }
};

class SE1_cut : public SubtourEliminationCut_heuristic {
public:
    //
    SE1_cut(std::string ch_name, Problem *prob) : SubtourEliminationCut_heuristic(ch_name, prob) {}
private:
    double get_LHS_givenS(double **x_ij, const std::set<int> &S);
};

class SE2_cut : public SubtourEliminationCut_heuristic {
public:
    //
    SE2_cut(std::string ch_name, Problem *prob) : SubtourEliminationCut_heuristic(ch_name, prob) {}
private:
    double get_LHS_givenS(double **x_ij, const std::set<int> &S);
};


#endif /* Cut_hpp */
