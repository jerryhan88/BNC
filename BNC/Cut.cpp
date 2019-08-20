//
//  Cut.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 1/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Cut.hpp"

#define NUM_ITER 25
#define MIN_CAR_SUBSET 3
#define PROBABILITY_ADD_NODE 0.5

std::set<edge> get_delta_S(const std::set<int> &S, const Problem &prmt) {
    std::set<edge> delta_S_P;
    for (int i: S) {
        for (int j: prmt.N) {
            if (S.find(j) != S.end())
                continue;
            edge ij{i, j};
            delta_S_P.insert(ij);
        }
    }
    std::set<edge> delta_S_M;
    for (int i: prmt.N) {
        if (S.find(i) != S.end())
            continue;
        for (int j: S) {
            edge ij{i, j};
            delta_S_M.insert(ij);
        }
    }
    //
    std::set<edge> delta_S(delta_S_P);
    delta_S.insert(delta_S_M.begin(), delta_S_M.end());
    return delta_S;
}

std::set<int> get_sigma_S(const std::set<int> &S, const Problem &prmt) {
    std::set<int> sigma_S;
    for (int i: S) {
        if (prmt.P.find(i) == prmt.P.end())
            continue;
        sigma_S.insert(i + (int) prmt.K.size());
    }
    return sigma_S;
}

std::set<int> get_pi_S(const std::set<int> &S, const Problem &prmt) {
    std::set<int> pi_S;
    for (int i: S) {
        if (prmt.D.find(i) == prmt.D.end())
            continue;
        pi_S.insert(i - (int) prmt.K.size());
    }
    return pi_S;
}

std::set<int> get_bar_S(const std::set<int> &S, const Problem &prmt) {
    std::set<int> bar_S;
    for (int i: prmt.N) {
        if (S.find(i) != S.end())
            continue;
        bar_S.insert(i);
    }
    return bar_S;
}

void get_SEC_aS(int SEC_no,
                const std::set<std::set<int>> &tabu_S,
                const std::set<int> &tabu_N,
                const std::set<int> &S0,
                int *min_i, double *min_lhs,
                double **x_ij, const Problem &prmt) {
    // Get a new subset by adding a new node
    for(int i: prmt.PD) {
        if(tabu_N.find(i) != tabu_N.end())
            continue;
        if (S0.find(i) != S0.end())
            continue;
        std::set<int> S1(S0);
        S1.insert(i);
        if (tabu_S.find(S1) != tabu_S.end())
            continue;
        double lhs = 0.0;
        if (SEC_no == 1) {
            get_SEC1_LHS(&lhs, x_ij, S1, prmt);
        } else{
            assert(SEC_no == 2);
            get_SEC2_LHS(&lhs, x_ij, S1, prmt);
        }
        if (lhs < *min_lhs) {
            *min_lhs = lhs;
            *min_i = i;
        }
        S1.clear();
    }
}

void get_SEC_dS(int SEC_no,
                const std::set<std::set<int>> &tabu_S,
                const std::set<int> &tabu_N,
                const std::set<int> &S0,
                int *min_i, double *min_lhs,
                double **x_ij, const Problem &prmt) {
    // Get a new subset by removing the existing node
    for(int i: S0) {
        std::set<int> S1(S0);
        S1.erase(i);
        if (tabu_S.find(S1) != tabu_S.end())
            continue;
        double lhs = 0.0;
        if (SEC_no == 1) {
            get_SEC1_LHS(&lhs, x_ij, S1, prmt);
        } else{
            assert(SEC_no == 2);
            get_SEC2_LHS(&lhs, x_ij, S1, prmt);
        }
        if (lhs < *min_lhs) {
            *min_lhs = lhs;
            *min_i = i;
        }
        S1.clear();
    }
}

std::set<std::set<int>> get_SEC(int SEC_no,
                                   const std::set<std::set<int>> &ES_SEC,
                                   double **_x_ij, const Problem &prmt) {
    /*
     * Valid subtour elimination constraint 1 (SEC_no = 1)
     *  \sum_{ i, j \in S} x_{i, j}
     *   + \sum_{i \in \bar{S} \cap \sigma(S)} \sum_{j \in S} x_{i, j}
     *   + \sum_{i \in \bar{S} \setminus \sigma(S)} \sum_{j \in S \cap \sigma(S)} x_{i, j}
     *  <= |S| - 1
     *
     * Seperation inequality; 2 * x(\bar{S}) + x(\delta(S)) = 2|S|
     *  \sum_{ i, j \in \delta(S)} x_{i, j}
     *   - 2 * \sum_{i \in \bar{S} \cap \sigma(S)} \sum_{j \in S} x_{i, j}
     *   - 2 * \sum_{i \in \bar{S} \setminus \sigma(S)} \sum_{j \in S \cap \sigma(S)} x_{i, j}
     *  < 2
     *
     * Valid subtour elimination constraint 2 (SEC_no = 2)
     *  \sum_{ i, j \in S} x_{i, j}
     *   + \sum_{i \in S} \sum_{j \in \bar{S} \cap \pi(S)} x_{i, j}
     *   + \sum_{i \in S \cap \pi(S)} \sum_{j \in \bar{S} \setminus \pi(S)} x_{i, j}
     *  <= |S| - 1
     *
     * Seperation inequality; 2 * x(\bar{S}) + x(\delta(S)) = 2|S|
     *  \sum_{ i, j \in \delta(S)} x_{i, j}
     *   - 2 * \sum_{i \in S} \sum_{j \in \bar{S} \cap \pi(S)} x_{i, j}
     *   - 2 * \sum_{i \in S \cap \pi(S)} \sum_{j \in \bar{S} \setminus \pi(S)} x_{i, j}
     *  < 2
     */
    
    std::set<std::set<int>> violated_subsets;
    //
    std::set<int> S0;
    std::set<int> tabu_N;
    std::set<std::set<int>> tabu_S;
    tabu_S.insert(ES_SEC.begin(), ES_SEC.end());
    //
    auto randIt = prmt.PD.begin();
    advance(randIt, rand() % prmt.PD.size());
    S0.insert(*randIt);
    int min_i;
    double min_lhs;
    for (int i = 0; i < NUM_ITER; i++) {
        min_i = -1;
        min_lhs = DBL_MAX;
        std::set<int> S1(S0);
        if (S0.size() <= MIN_CAR_SUBSET) {
            get_SEC_aS(SEC_no,
                       tabu_S,
                       tabu_N, S0, &min_i, &min_lhs, _x_ij, prmt);
            S1.insert(min_i);
        } else if (get_random_number() < PROBABILITY_ADD_NODE) {
            get_SEC_aS(SEC_no,
                       tabu_S,
                       tabu_N, S0, &min_i, &min_lhs, _x_ij, prmt);
            S1.insert(min_i);
        } else {
            get_SEC_dS(SEC_no,
                       tabu_S,
                       tabu_N, S0, &min_i, &min_lhs, _x_ij, prmt);
            S1.erase(min_i);
            tabu_N.insert(min_i);
        }
        if (min_i == -1) {
            continue;
        }
//        printf("%d: size %d, min_i: %d, min_lhs: %.2f\n", i, (int) S1.size(), min_i, min_lhs);
        if (min_lhs < 2.0) {
            violated_subsets.insert(S1);
            tabu_S.insert(S1);
        }
        S0.clear();
        S0 = S1;
    }
    return violated_subsets;
}


void SubtourEliminationCut_heuristic::indentify_violatedConstrains(double **x_ij) {
    assert(new_sets.size() == 0);
    //
    std::set<int> S0;
    std::set<int> tabu_N;
    std::set<std::set<int>> tabu_S;
    tabu_S.insert(existing_sets.begin(), existing_sets.end());
    //
    auto randIt = (*prob).PD.begin();
    advance(randIt, rand() % (*prob).PD.size());
    S0.insert(*randIt);
    int min_nid;
    double min_lhs;
    for (int i = 0; i < NUM_ITER; i++) {
        min_nid = -1;
        min_lhs = DBL_MAX;
        std::set<int> S1(S0);
        if (S0.size() <= MIN_CAR_SUBSET) {
            find_new_nid_wMinLHS_wAddition(x_ij, S0, tabu_S, tabu_N,
                                           &min_nid, &min_lhs);
            S1.insert(min_nid);
        } else if (get_random_number() < PROBABILITY_ADD_NODE) {
            find_new_nid_wMinLHS_wAddition(x_ij, S0, tabu_S, tabu_N,
                                           &min_nid, &min_lhs);
            S1.insert(min_nid);
        } else {
            find_included_nid_wMinLHS_wDeletion(x_ij, S0, tabu_S, tabu_N,
                                                          &min_nid, &min_lhs);
            S1.erase(min_nid);
            tabu_N.insert(min_nid);
        }
        if (min_nid == -1) {
            continue;
        }
        if (min_lhs < 2.0) {
            new_sets.insert(S1);
            tabu_S.insert(S1);
        }
        S0.clear();
        S0 = S1;
    }
}


void SubtourEliminationCut_heuristic::find_new_nid_wMinLHS_wAddition(double **x_ij, const std::set<int> &S0,
                                std::set<std::set<int>> &tabu_S,
                                            std::set<int> &tabu_N,
                                            int *min_nid, double *min_lhs) {
    // Get a new subset by adding a new node
    for(int i: (*prob).PD) {
        if(tabu_N.find(i) != tabu_N.end())
            continue;
        if (S0.find(i) != S0.end())
            continue;
        std::set<int> S1(S0);
        S1.insert(i);
        if (tabu_S.find(S1) != tabu_S.end())
            continue;
        double lhs = get_LHS_givenS(x_ij, S1);
        if (lhs < *min_lhs) {
            *min_lhs = lhs;
            *min_nid = i;
        }
        S1.clear();
    }
}

void SubtourEliminationCut_heuristic::find_included_nid_wMinLHS_wDeletion(double **x_ij, const std::set<int> &S0,
                                            std::set<std::set<int>> &tabu_S,
                                                 std::set<int> &tabu_N,
                                                 int *min_nid, double *min_lhs) {
    // Get a new subset by removing the existing node
    for(int i: S0) {
        std::set<int> S1(S0);
        S1.erase(i);
        if (tabu_S.find(S1) != tabu_S.end())
            continue;
        double lhs = get_LHS_givenS(x_ij, S1);
        if (lhs < *min_lhs) {
            *min_lhs = lhs;
            *min_nid = i;
        }
        S1.clear();
    }
}


double SE1_cut::get_LHS_givenS(double **x_ij, const std::set<int> &S) {
    /*
     * Valid subtour elimination constraint
     *  \sum_{ i, j \in S} x_{i, j}
     *   + \sum_{i \in \bar{S} \cap \sigma(S)} \sum_{j \in S} x_{i, j}
     *   + \sum_{i \in \bar{S} \setminus \sigma(S)} \sum_{j \in S \cap \sigma(S)} x_{i, j}
     *  <= |S| - 1
     *
     * Seperation inequality; 2 * x(\bar{S}) + x(\delta(S)) = 2|S|
     *  \sum_{ i, j \in \delta(S)} x_{i, j}
     *   - 2 * \sum_{i \in \bar{S} \cap \sigma(S)} \sum_{j \in S} x_{i, j}
     *   - 2 * \sum_{i \in \bar{S} \setminus \sigma(S)} \sum_{j \in S \cap \sigma(S)} x_{i, j}
     *  < 2
     */
    //
    
    double lhs = 0.0;
    //
    std::set<edge> delta_S = get_delta_S(S, *prob);
    std::set<int> bar_S = get_bar_S(S, *prob);
    std::set<int> sigma_S = get_sigma_S(S, *prob);
    /*
     * \sum_{ i, j \in \delta(S)} x_{i, j}
     */
    for (edge a: delta_S) {
        lhs += x_ij[a.first][a.second];
    }
    /*
     * -2 * \sum_{i \in \bar{S} \cap \sigma(S)} \sum_{j \in S} x_{i, j}
     */
    for (int i: bar_S) {
        if (sigma_S.find(i) == sigma_S.end())
            continue;
        for (int j: S) {
            lhs -= 2 * x_ij[i][j];
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
            lhs -= 2 * x_ij[i][j];
        }
    }
    //
    return lhs;
}


double SE2_cut::get_LHS_givenS(double **x_ij, const std::set<int> &S) {
    /*
     * Valid subtour elimination constraint
     *  \sum_{ i, j \in S} x_{i, j}
     *   + \sum_{i \in S} \sum_{j \in \bar{S} \cap \pi(S)} x_{i, j}
     *   + \sum_{i \in S \cap \pi(S)} \sum_{j \in \bar{S} \setminus \pi(S)} x_{i, j}
     *  <= |S| - 1
     *
     * Seperation inequality; 2 * x(\bar{S}) + x(\delta(S)) = 2|S|
     *  \sum_{ i, j \in \delta(S)} x_{i, j}
     *   - 2 * \sum_{i \in S} \sum_{j \in \bar{S} \cap \pi(S)} x_{i, j}
     *   - 2 * \sum_{i \in S \cap \pi(S)} \sum_{j \in \bar{S} \setminus \pi(S)} x_{i, j}
     *  < 2
     */
    //
    
    double lhs = 0.0;
    //
    std::set<edge> delta_S = get_delta_S(S, *prob);
    std::set<int> bar_S = get_bar_S(S, *prob);
    std::set<int> pi_S = get_pi_S(S, *prob);
    /*
     * \sum_{ i, j \in \delta(S)} x_{i, j}
     */
    for (edge a: delta_S) {
        lhs += x_ij[a.first][a.second];
    }
    /*
     * - 2 * \sum_{i \in S} \sum_{j \in \bar{S} \cap \pi(S)} x_{i, j}
     */
    for (int i: S) {
        for (int j: bar_S) {
            if (pi_S.find(j) == pi_S.end())
                continue;
            lhs -= 2 * x_ij[i][j];
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
            lhs -= 2 * x_ij[i][j];
        }
    }
    //
    return lhs;
}
