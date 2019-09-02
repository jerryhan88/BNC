//
//  Problem.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 1/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Problem.hpp"

#define DEFAULT_REWARD 1.0
#define DEFAULT_VOLUME 1.0
#define DEFAULT_WEIGHT 1.0
#define DEAULT_TW_BEGIN 0.0
#define DEAULT_TW_END DBL_MAX


double get_euclideanDistance(int i, int j, double *X, double *Y) {
    double dx = X[i] - X[j];
    double dy = Y[i] - Y[j];
    //
    return sqrt(dx * dx + dy * dy);
}

Problem gen_problemInstance(int numTasks, int num_rrPoints, int maxReward,
                              double bv, double bw, double bu) {
    
    double reward[numTasks], volume[numTasks], weight[numTasks];
    for (int k = 0; k < numTasks; k++) {
        reward[k] = DEFAULT_REWARD;
        volume[k] = DEFAULT_VOLUME;
        weight[k] = DEFAULT_WEIGHT;
    }
    int numNodes = (numTasks + 1) * 2 + num_rrPoints;
    double X[numNodes], Y[numNodes];
    for (int i = 0; i < numNodes; i++) {
        X[i] = get_random_number();
        Y[i] = get_random_number();
    }
    double **distance = new double*[numNodes];
    for (int i = 0; i < numNodes; i++) {
        distance[i] = new double[numNodes];
        for (int j = 0; j < numNodes; j++) {
//            distance[i][j] = get_euclideanDistance(i, j, X, Y);
            distance[i][j] = i != j ? 0.01 : 0.00;
        }
    }
    double **timeWindow = new double*[numNodes];
    for (int i = 0; i < numNodes; i++) {
        timeWindow[i] = new double[2];
        timeWindow[i][0] = DEAULT_TW_BEGIN;
        timeWindow[i][1] = DEAULT_TW_END;
    }
    //
    Problem pi(numTasks,
                    reward, volume, weight,
                 numNodes,
                    distance, timeWindow,
                 bv, bw, bu);
    
    for (int i = 0; i < numNodes; i++) {
        delete [] distance[i];
        delete [] timeWindow[i];
    }
    delete [] distance;
    delete [] timeWindow;
    //
    return pi;
}

Problem::Problem(int numTasks,
                        double *reward, double *volume, double *weight,
                     int numNodes,
                        double **distance, double **timeWindow,
                     double bv, double bw, double bu){
    K.reserve(numTasks);
    PD.reserve(2 * numTasks);
    for (int i = 0; i < numTasks; i++) {
        K.push_back(i);
    }
    r_k = new double[K.size()];
    v_k = new double[K.size()];
    w_k = new double[K.size()];
    h_k = new int[K.size()];
    n_k = new int[K.size()];
    for (int k: K) {
        r_k[k] = reward[k];
        v_k[k] = volume[k];
        w_k[k] = weight[k];
        h_k[k] = k + 1;
        n_k[k] = (int) K.size() + k + 1;
        P.insert(h_k[k]);
        PD.push_back(h_k[k]);
        D.insert(n_k[k]);
        PD.push_back(n_k[k]);
    }
    //
    o = 0; d = (int) PD.size() + 1;
    S.push_back(o);
    for (int i = d + 1; i < numNodes; i++) {
        S.push_back(i);
    }
    S.push_back(d);
    N.insert(N.end(), PD.begin(), PD.end());
    N.insert(N.end(), S.begin(), S.end());
    sort(N.begin(),N.end());
    assert(N.size() == numNodes);
    //
    t_ij = new double*[N.size()];
    for (int i = 0; i < N.size(); i++)
        t_ij[i] = new double[N.size()];
    al_i = new double[N.size()];
    be_i = new double[N.size()];
    double max_dist = -1.0;
    for (int i: N) {
        for (int j: N) {
            t_ij[i][j] = distance[i][j];
        }
        double i_max_dist = *std::max_element(t_ij[i], t_ij[i] + N.size());
        if (max_dist < i_max_dist) {
            max_dist = i_max_dist;
        }
        al_i[i] = timeWindow[i][0];
        be_i[i] = timeWindow[i][1];
    }
    c_ij = new int*[N.size()];
    for (int i = 0; i < N.size(); i++)
        c_ij[i] = new int[N.size()];
    for (int i = 0; i < S.size(); i++) {
        for (int j = 0; j < S.size(); j++) {
            if (j < i)
                continue;
            c_ij[S[i]][S[j]] = 1;
            c_ij[S[j]][S[i]] = 0;
        }
    }
    M = N.size() * N.size() * max_dist;
    if (M < 1.0) {
        M = 1.0;
    }
    this->bv = bv; this->bw = bw; this->bu = bu;
}

Problem::~Problem() {
    std::vector<int*> ipV = {h_k, n_k};
    for(auto p: ipV) {
        delete p;
    }
    std::vector<double*> dpV = {r_k, v_k, w_k,
        al_i, be_i};
    for(auto p: dpV) {
        delete p;
    }
    for (int i = 0; i < N.size(); i++) {
        delete [] c_ij[i];
        delete [] t_ij[i];
    }
    delete [] c_ij;
    delete [] t_ij;
    //
    K.clear(); PD.clear(); S.clear(); N.clear();
    P.clear(); D.clear();
}


void Problem::write_json(std::string ofpath) {
    nlohmann::json prob_json;
    //
    prob_json["bv"] = bv;
    prob_json["bw"] = bw;
    prob_json["bu"] = bu;
    //
    prob_json["K"] = nlohmann::json (K);
    prob_json["PD"] = nlohmann::json (PD);
    prob_json["S"] = nlohmann::json (S);
    prob_json["N"] = nlohmann::json (N);
    prob_json["P"] = nlohmann::json (P);
    prob_json["D"] = nlohmann::json (D);
    prob_json["o"] = o;
    prob_json["d"] = d;
    std::vector<int> _r_k(K.size()), _v_k(K.size()), _w_k(K.size()), _h_k(K.size()), _n_k(K.size());
    for (int k: K) {
        _r_k[k] = r_k[k];
        _v_k[k] = v_k[k];
        _w_k[k] = w_k[k];
        _h_k[k] = h_k[k];
        _n_k[k] = n_k[k];
    }
    prob_json["r_k"] = nlohmann::json (_r_k);
    prob_json["v_k"] = nlohmann::json (_v_k);
    prob_json["w_k"] = nlohmann::json (_w_k);
    prob_json["h_k"] = nlohmann::json (_h_k);
    prob_json["n_k"] = nlohmann::json (_n_k);
    //
    std::vector<double> _al_i(N.size()), _be_i(N.size());
    std::vector<std::vector<double>> _t_ij(N.size());
    std::vector<std::vector<int>> _c_ij(N.size());
    for (int i: N) {
        _al_i[i] = al_i[i];
        _be_i[i] = be_i[i];
        std::vector<double> _t_i(N.size());
        std::vector<int> _c_i(N.size());
        for (int j: N) {
            _t_i[j] = t_ij[i][j];
            _c_i[j] = c_ij[i][j];
        }
        _t_ij[i] = _t_i;
        _c_ij[i] = _c_i;
    }
    prob_json["al_i"] = nlohmann::json (_al_i);
    prob_json["be_i"] = nlohmann::json (_be_i);
    prob_json["t_ij"] = nlohmann::json (_t_ij);
    prob_json["c_ij"] = nlohmann::json (_c_ij);
    //
    prob_json["M"] = M;
    
    std::ofstream os(ofpath);
    os << prob_json <<std::endl;
    
}
