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

std::set<arc> get_delta_S(std::set<int> S, const Problem &prmt) {
    std::set<arc> delta_S_P;
    for (int i: S) {
        for (int j: prmt.N) {
            if (S.find(j) != S.end())
                continue;
            arc ij{i, j};
            delta_S_P.insert(ij);
        }
    }
    std::set<arc> delta_S_M;
    for (int i: prmt.N) {
        if (S.find(i) != S.end())
            continue;
        for (int j: S) {
            arc ij{i, j};
            delta_S_M.insert(ij);
        }
    }
    //
    std::set<arc> delta_S(delta_S_P);
    delta_S.insert(delta_S_M.begin(), delta_S_M.end());
    return delta_S;
}

std::set<int> get_sigma_S(std::set<int> S, const Problem &prmt) {
    std::set<int> sigma_S;
    for (int i: S) {
        if (prmt.P.find(i) == prmt.P.end())
            continue;
        sigma_S.insert(i + (int) prmt.K.size());
    }
    return sigma_S;
}

std::set<int> get_pi_S(std::set<int> S, const Problem &prmt) {
    std::set<int> pi_S;
    for (int i: S) {
        if (prmt.D.find(i) == prmt.D.end())
            continue;
        pi_S.insert(i - (int) prmt.K.size());
    }
    return pi_S;
}

std::set<int> get_bar_S(std::set<int> S, const Problem &prmt) {
    std::set<int> bar_S;
    for (int i: prmt.N) {
        if (S.find(i) != S.end())
            continue;
        bar_S.insert(i);
    }
    return bar_S;
}

void get_SEC_aS(int SEC_no,
                const std::set<int> &tabuS,
                const std::set<int> &S0,
                int *min_i, double *min_lhs,
                double **x_ij, const Problem &prmt) {
    // Get a new subset by adding a new node
    for(int i: prmt.PD) {
        if(tabuS.find(i) != tabuS.end())
            continue;
        if (S0.find(i) != S0.end())
            continue;
        std::set<int> S1(S0);
        S1.insert(i);
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
                const std::set<int> &tabuS,
                const std::set<int> &S0,
                int *min_i, double *min_lhs,
                double **x_ij, const Problem &prmt) {
    // Get a new subset by removing the existing node
    for(int i: S0) {
        std::set<int> S1(S0);
        S1.erase(i);
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

std::vector<std::set<int>> get_SEC(int SEC_no, double **_x_ij, const Problem &prmt) {
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
    
    std::vector<std::set<int>> violated_subsets;
    //
    std::set<int> S0;
    std::set<int> tabuS;
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
            get_SEC_aS(SEC_no, tabuS, S0, &min_i, &min_lhs, _x_ij, prmt);
            S1.insert(min_i);
        } else if (get_random_number() < PROBABILITY_ADD_NODE) {
            get_SEC_aS(SEC_no, tabuS, S0, &min_i, &min_lhs, _x_ij, prmt);
            S1.insert(min_i);
        } else {
            get_SEC_dS(SEC_no, tabuS, S0, &min_i, &min_lhs, _x_ij, prmt);
            S1.erase(min_i);
            tabuS.insert(min_i);
        }
        if (min_i == -1) {
            continue;
        }
        printf("%d: size %d, min_i: %d, min_lhs: %.2f\n", i, (int) S1.size(), min_i, min_lhs);
        if (min_lhs < 2.0) {
            violated_subsets.push_back(S1);
        }
        S0.clear();
        S0 = S1;
    }
    return violated_subsets;
}
