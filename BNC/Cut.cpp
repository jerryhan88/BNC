//
//  Cut.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 1/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Cut.hpp"


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

void SE0_cut::indentify_violatedConstrains(double **x_ij) {
    int V = (int) prob->N.size();
    double **rGraph = gen_graphInstance(x_ij, V);
    int parent[V];
    double max_flow;
    //
    bool isFirst = true;
    for (int i: prob->S) {
        if (i == prob->o || i == prob->d)
            continue;
        for (int j: prob->S) {
            if (j == prob->o || j == prob->d)
                continue;
            if (i == j)
                continue;
            if (!isFirst) {
                for (int i = 0; i < V; i++) {
                    for (int j = 0; j < V; j++) {
                        rGraph[i][j] = x_ij[i][j];
                    }
                }
            } else {
                isFirst = false;
            }
            for (int i = 0; i < V; i++) {
                parent[i] = -1;
            }
            max_flow = 0.0;
            //
            sendMaxFlow(rGraph, V, i, j, parent, &max_flow);
            bool visited[V];
            std::memset(visited, false, sizeof(visited));
            reachability_dfs(rGraph, V, i, visited);
            std::set<int> S1;
            for (int i = 0; i < V; i++) {
                if(visited[i]) {
                    S1.insert(i);
                }
            }
            if (S1.size() <= 2)
                continue;
            if (existing_sets.find(S1) == existing_sets.end()) {
                new_sets.insert(S1);
            }
        }
    }
    for (int i = 0; i < V; i++) {
        delete [] rGraph[i];
    }
    delete [] rGraph;
}

void SE0_cut::add_cuts(GRBModel *grbModel, GRBVar **x_ij) {
    char buf[DEFAULT_BUFFER_SIZE];
    for (std::set<int> S: new_sets) {
        GRBLinExpr lhs = 0.0;
        for (int i: S) {
            for (int j: S) {
                lhs += x_ij[i][j];
            }
        }
//        std::set<edge> delta_S = get_delta_S(S, *prob);
//        for (edge a: delta_S) {
//            lhs += x_ij[a.first][a.second];
//        }
        num_generatedCuts += 1;
        sprintf(buf, "%s[%d]", ch_name.c_str(), num_generatedCuts);
//        grbModel->addConstr(lhs >= 2, buf);
        
        grbModel->addConstr(lhs <= (int) S.size() - 2, buf);
    }
    grbModel->update();
    existing_sets.insert(new_sets.begin(), new_sets.end());
    new_sets.clear();
}

void SE1_cut::indentify_violatedConstrains(double **x_ij) {
    indentify_seViolatedConstrains(x_ij,
                                   new_sets,
                                   existing_sets,
                                   prob,
                                   get_LHS_givenS<double, double>);
}

void SE2_cut::indentify_violatedConstrains(double **x_ij) {
    indentify_seViolatedConstrains(x_ij,
                                   new_sets,
                                   existing_sets,
                                   prob,
                                   get_LHS_givenS<double, double>);
}

void SE1_cut::add_cuts(GRBModel *grbModel, GRBVar **x_ij) {
    add_SE_cuts(grbModel, x_ij,
                new_sets,
                &num_generatedCuts,
                ch_name,
                existing_sets,
                prob,
                get_LHS_givenS<GRBVar, GRBLinExpr>);
}

void SE2_cut::add_cuts(GRBModel *grbModel, GRBVar **x_ij) {
    add_SE_cuts(grbModel, x_ij,
                new_sets,
                &num_generatedCuts,
                ch_name,
                existing_sets,
                prob,
                get_LHS_givenS<GRBVar, GRBLinExpr>);
}


cut_composer::cut_composer(Problem *prob, std::vector<std::string> cut_names) {
    for (std::string cn: cut_names) {
        if (cn == "SE0") {
            chs.push_back(new SE0_cut("SE0", prob));
        } else if (cn == "SE1") {
            chs.push_back(new SE1_cut("SE1", prob));
        } else if (cn == "SE2") {
            chs.push_back(new SE2_cut("SE2", prob));
        } else {
            assert(false);
        }
    }
}

cut_composer* cut_composer::clone() {
    cut_composer *cc = new cut_composer();
    for (cut_handler *ch : chs) {
        cc->chs.push_back(ch->clone());
    }
    return cc;
}

void cut_composer::solve_seperation_problem(double **x_ij) {
    for (cut_handler* ch: chs) {
        ch->indentify_violatedConstrains(x_ij);
    }
}

void cut_composer::add_cuts(GRBModel *grbModel, GRBVar **x_ij) {
    for (cut_handler* ch: chs) {
        ch->add_cuts(grbModel, x_ij);
    }
}

std::vector<int> cut_composer::get_numViolatedCnsts() {
    std::vector<int> vec;
    for (cut_handler* ch: chs) {
        vec.push_back(ch->get_numIdentifiedConstraints());
    }
    return vec;
}
