//
//  Problem.hpp
//  BNC
//
//  Created by Chung-Kyun HAN on 1/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Problem_hpp
#define Problem_hpp

#include <stdio.h>
#include <vector>
#include <set>
#include <cmath>
#include <float.h>

#include "nlohmann/json.hpp"
#include "Etc.hpp"


class Problem {
public:
    std::string problemName;
    double bv, bw, bu;
    //
    std::vector<int> K, PD, S, N;
    std::set<int> P, D;
    int o, d;
    int *h_k, *n_k, **c_ij;
    double *r_k, *v_k, *w_k;
    double *al_i, *be_i;
    double **t_ij;
    double M;
    //
    Problem () {};
    Problem(int,
            double *, double *, double *,
            int,
            double **, double **,
            double, double, double);
    ~Problem();
    //
    void write_json(std::string ofpath);
    static Problem* read_json(std::string ifpath) {
        std::ifstream is(ifpath);
        nlohmann::json prob_json;
        is >> prob_json;
        
        Problem *prob = new Problem();
        prob->problemName = prob_json["problemName"];
        for (int k: prob_json["K"]) {
            prob->K.push_back(k);
        }
        for (int i: prob_json["PD"]) {
            prob->PD.push_back(i);
        }
        for (int i: prob_json["S"]) {
            prob->S.push_back(i);
        }
        for (int i: prob_json["N"]) {
            prob->N.push_back(i);
        }
        for (int i: prob_json["P"]) {
            prob->P.insert(i);
        }
        for (int i: prob_json["D"]) {
            prob->D.insert(i);
        }
        prob->o = prob_json["o"];
        prob->d = prob_json["d"];
        prob->r_k = new double[prob->K.size()];
        prob->v_k = new double[prob->K.size()];
        prob->w_k = new double[prob->K.size()];
        prob->h_k = new int[prob->K.size()];
        prob->n_k = new int[prob->K.size()];
        for (int k: prob->K) {
            prob->r_k[k] = prob_json["r_k"][k];
            prob->v_k[k] = prob_json["v_k"][k];
            prob->w_k[k] = prob_json["w_k"][k];
            prob->h_k[k] = prob_json["h_k"][k];
            prob->n_k[k] = prob_json["n_k"][k];
        }
        prob->t_ij = new double*[prob->N.size()];
        prob->c_ij = new int*[prob->N.size()];
        for (int i = 0; i < prob->N.size(); i++) {
            prob->t_ij[i] = new double[prob->N.size()];
            prob->c_ij[i] = new int[prob->N.size()];
        }
        prob->al_i = new double[prob->N.size()];
        prob->be_i = new double[prob->N.size()];
        for (int i: prob->N) {
            prob->al_i[i] = prob_json["al_i"][i];
            prob->be_i[i] = prob_json["be_i"][i];
            for (int j: prob->N) {
                prob->t_ij[i][j] = prob_json["t_ij"][i][j];
                prob->c_ij[i][j] = prob_json["c_ij"][i][j];
            }
        }
        prob->M = prob_json["M"];
        prob->bv = prob_json["bv"];
        prob->bw = prob_json["bw"];
        prob->bu = prob_json["bu"];
        
        return prob;
    }
};

Problem gen_problemInstance(int numTasks, int num_rrPoints, int maxReward,
                            double bv, double bw, double bu);


#endif /* Problem_hpp */

