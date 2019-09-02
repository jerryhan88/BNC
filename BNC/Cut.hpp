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
#include <string>


#include "Problem.hpp"
#include "Etc.hpp"
#include "gurobi_c++.h"
#include "NetworkFlow/graph.hpp"
#include "NetworkFlow/FordFulkersonAlgo.hpp"


typedef std::pair<int, int> edge;
#define NUM_ITER 25
#define MIN_CAR_SUBSET 3
#define PROBABILITY_ADD_NODE 0.5
#define DEFAULT_BUFFER_SIZE 2048

std::set<edge> get_delta_S(const std::set<int> &S, const Problem &prmt);
std::set<int> get_sigma_S(const std::set<int> &S, const Problem &prmt);
std::set<int> get_pi_S(const std::set<int> &S, const Problem &prmt);
std::set<int> get_bar_S(const std::set<int> &S, const Problem &prmt);


class cut_handler {
public:
    std::string ch_name;
    Problem *prob;
    signed int num_generatedCuts;
    //
    cut_handler(std::string ch_name, Problem *prob){
        this->ch_name = ch_name;
        this->prob = prob;
        num_generatedCuts = 0;
    }
    //
    virtual void indentify_violatedConstrains(double **x_ij) {
        throw "Should override indentify_violatedConstrains()";
    }
    virtual void add_cuts(GRBModel *grbModel, GRBVar **x_ij) {
        throw "Should override add_cuts()";
    }
    virtual cut_handler* clone() {
        throw "Should override clone()";
    }
    virtual signed int get_numIdentifiedConstraints() {
        throw "Should override get_numIdentifiedConstraints()";
    }
};


class SubtourElimination_cut : public cut_handler {
public:
    std::set< std::set<int> > existing_sets;
    std::set< std::set<int> > new_sets;
    //
    SubtourElimination_cut(std::string ch_name, Problem *prob): cut_handler(ch_name, prob) {
    }
    ~SubtourElimination_cut() {
        existing_sets.clear();
        new_sets.clear();
    }
    //
    signed int get_numIdentifiedConstraints() {
        return (int) new_sets.size();
    }
};

class SE0_cut : public SubtourElimination_cut {
public:
    //
    SE0_cut(std::string ch_name, Problem *prob) : SubtourElimination_cut(ch_name, prob) {}
    SE0_cut* clone(){
        SE0_cut *sec = new SE0_cut(ch_name, prob);
        sec->existing_sets.insert(existing_sets.begin(), existing_sets.end());
        sec->num_generatedCuts = num_generatedCuts;
        return sec;
    }
    void indentify_violatedConstrains(double **x_ij);
    void add_cuts(GRBModel *grbModel, GRBVar **x_ij);
};

class SE1_cut : public SubtourElimination_cut {
public:
    //
    SE1_cut(std::string ch_name, Problem *prob) : SubtourElimination_cut(ch_name, prob) {}
    SE1_cut* clone(){
        SE1_cut *sec = new SE1_cut(ch_name, prob);
        sec->existing_sets.insert(existing_sets.begin(), existing_sets.end());
        sec->num_generatedCuts = num_generatedCuts;
        return sec;
    }
    void indentify_violatedConstrains(double **x_ij);
    void add_cuts(GRBModel *grbModel, GRBVar **x_ij);
    //
    template<typename T1, typename T2>
    static void get_LHS_givenS(T1 **x_ij, const std::set<int> &S, T2 *lhs,
                               Problem *prob) {
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
        
        std::set<edge> delta_S = get_delta_S(S, *prob);
        std::set<int> bar_S = get_bar_S(S, *prob);
        std::set<int> sigma_S = get_sigma_S(S, *prob);
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
};

class SE2_cut : public SubtourElimination_cut {
public:
    //
    SE2_cut(std::string ch_name, Problem *prob) : SubtourElimination_cut(ch_name, prob) {}
    SE2_cut* clone(){
        SE2_cut *sec = new SE2_cut(ch_name, prob);
        sec->existing_sets.insert(existing_sets.begin(), existing_sets.end());
        sec->num_generatedCuts = num_generatedCuts;
        return sec;
    }
    void indentify_violatedConstrains(double **x_ij);
    void add_cuts(GRBModel *grbModel, GRBVar **x_ij);
    template<typename T1, typename T2>
    static void get_LHS_givenS(T1 **x_ij, const std::set<int> &S, T2 *lhs,
                               Problem *prob) {
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
        
        std::set<edge> delta_S = get_delta_S(S, *prob);
        std::set<int> bar_S = get_bar_S(S, *prob);
        std::set<int> pi_S = get_pi_S(S, *prob);
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
};

template <typename F>
void find_new_nid_wMinLHS_wAddition(double **x_ij, const std::set<int> &S0,
                                    std::set<std::set<int>> &tabu_S,
                                    std::set<int> &tabu_N,
                                    int *min_nid, double *min_lhs,
                                    Problem *prob,
                                    F get_LHS_givenS) {
    // Get a new subset by adding a new node
    for(int i: prob->PD) {
        if(tabu_N.find(i) != tabu_N.end())
            continue;
        if (S0.find(i) != S0.end())
            continue;
        std::set<int> S1(S0);
        S1.insert(i);
        if (tabu_S.find(S1) != tabu_S.end())
            continue;
        double lhs = 0.0;
        get_LHS_givenS(x_ij, S1, &lhs, prob);
        if (lhs < *min_lhs) {
            *min_lhs = lhs;
            *min_nid = i;
        }
        S1.clear();
    }
}

template <typename F>
void find_included_nid_wMinLHS_wDeletion(double **x_ij, const std::set<int> &S0,
                                         std::set<std::set<int>> &tabu_S,
                                         std::set<int> &tabu_N,
                                         int *min_nid, double *min_lhs,
                                         Problem *prob,
                                         F get_LHS_givenS) {
    // Get a new subset by removing the existing node
    for(int i: S0) {
        std::set<int> S1(S0);
        S1.erase(i);
        if (tabu_S.find(S1) != tabu_S.end())
            continue;
        double lhs = 0.0;
        get_LHS_givenS(x_ij, S1, &lhs, prob);
        if (lhs < *min_lhs) {
            *min_lhs = lhs;
            *min_nid = i;
        }
        S1.clear();
    }
}

template <typename F>
void indentify_seViolatedConstrains(double **x_ij,
                        std::set< std::set<int> > &new_sets,
                        std::set< std::set<int> > &existing_sets,
                        Problem *prob,
                        F get_LHS_givenS) {
    assert(new_sets.size() == 0);
    //
    std::set<int> S0;
    std::set<int> tabu_N;
    std::set<std::set<int>> tabu_S;
    tabu_S.insert(existing_sets.begin(), existing_sets.end());
    //
    
//    auto randIt = (*prob).PD.begin();
//    advance(randIt, rand() % (*prob).PD.size());
//    S0.insert(*randIt);
    
    S0.insert((*prob).S[1]);
    int min_nid;
    double min_lhs;
    for (int i = 0; i < NUM_ITER; i++) {
        min_nid = -1;
        min_lhs = DBL_MAX;
        std::set<int> S1(S0);
        if (S0.size() <= MIN_CAR_SUBSET) {
            find_new_nid_wMinLHS_wAddition(x_ij, S0, tabu_S, tabu_N,
                                           &min_nid, &min_lhs,
                                           prob,
                                           get_LHS_givenS);
            S1.insert(min_nid);
        } else if (get_random_number() < PROBABILITY_ADD_NODE) {
            find_new_nid_wMinLHS_wAddition(x_ij, S0, tabu_S, tabu_N,
                                           &min_nid, &min_lhs,
                                           prob,
                                           get_LHS_givenS);
            S1.insert(min_nid);
        } else {
//            find_included_nid_wMinLHS_wDeletion(x_ij, S0, tabu_S, tabu_N,
//                                                &min_nid, &min_lhs,
//                                                prob,
//                                                get_LHS_givenS);
//            S1.erase(min_nid);
//            tabu_N.insert(min_nid);
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

template <typename F>
void add_SE_cuts(GRBModel *grbModel, GRBVar **x_ij,
                 std::set< std::set<int> > &new_sets,
                 signed int *num_generatedCuts,
                 std::string ch_name,
                 std::set< std::set<int> > &existing_sets,
                 Problem *prob,
                 F get_LHS_givenS) {
    char buf[DEFAULT_BUFFER_SIZE];
    for (std::set<int> S: new_sets) {
        GRBLinExpr lhs = 0.0;
        get_LHS_givenS(x_ij, S, &lhs, prob);
        *num_generatedCuts += 1;
        sprintf(buf, "%s[%d]", ch_name.c_str(), *num_generatedCuts);
        grbModel->addConstr(lhs >= 2, buf);
    }
    existing_sets.insert(new_sets.begin(), new_sets.end());
    new_sets.clear();
}


class cut_composer {
public:
    std::vector<cut_handler*> chs;
    //
    cut_composer() {};
    cut_composer(Problem *prob, std::vector< std::string > cut_names);
    ~cut_composer() {
        for (cut_handler *ch : chs) {
            delete ch;
        }
        chs.clear();
    }
    //
    cut_composer* clone();
    void solve_seperation_problem(double **x_ij);
    void add_cuts(GRBModel *grbModel, GRBVar **x_ij);
    std::vector<int> get_numViolatedCnsts();
};


#endif /* Cut_hpp */
