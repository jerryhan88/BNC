//
//  PDPTW_ILP.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 22/7/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "PDPTW_ILP.hpp"

#define DEFAULT_BUFFER_SIZE 20

typedef pair<int, int> arc;


double get_random_number() {
    return ((double) rand()) / RAND_MAX;
}

double distance(double *x, double *y, int i, int j) {
    double dx = x[i] - x[j];
    double dy = y[i] - y[j];
    
    return sqrt(dx * dx + dy * dy);
}


void ex1(int *num_tasks, int *num_rr_points, int *max_reward,
         double *bv, double *bw, double *ul) {
    srand(1);
    *num_tasks = 5; *num_rr_points = 2; *max_reward = 3;
    *bv = 3.0; *bw = 3.0; *ul = 1.0 * 10;
}

class Parameter {
public:
    double bv, bw, ul;
    //
    vector<int> K, PD, S, N;
    set<int> F, P, D;
    int o, d, *h_k, *n_k, **c_ij;
    double *r_k, *v_k, *w_k;
    double *al_i, *be_i, *ga_i;
    double **t_ij;
    double M;
    //
    Parameter(int, int, int,
              double, double, double);
    ~Parameter();
};
Parameter::Parameter(int num_tasks, int num_rr_points, int max_reward,
                     double bv, double bw, double ul) {
    this->bv = bv;
    this->bw = bw;
    this->ul = ul;
    K.reserve(num_tasks);
    PD.reserve(2 * num_tasks);
    for (int i = 0; i < num_tasks; i++) {
        K.push_back(i);
    }
    r_k = new double[K.size()];
    v_k = new double[K.size()];
    w_k = new double[K.size()];
    h_k = new int[K.size()];
    n_k = new int[K.size()];
    for (int k: K) {
        r_k[k] = (int) (get_random_number() * max_reward);
        v_k[k] = 1.0;
        w_k[k] = 1.0;
        h_k[k] = k + 1;
        n_k[k] = (int) K.size() + k + 1;
        P.insert(h_k[k]);
        PD.push_back(h_k[k]);
        D.insert(n_k[k]);
        PD.push_back(n_k[k]);
    }
    for (int k: K) {
        if (get_random_number() < 0.5)
            continue;
        F.insert(k);
        if (F.size() == 2)
            break;
    }
    //
    o = 0; d = (int) PD.size() + 1;
    S.push_back(o);
    S.push_back(d);
    for (int i = 0; i < num_rr_points; i++) {
        S.push_back(d + (i + 1));
    }
    N.insert(N.end(), PD.begin(), PD.end());
    N.insert(N.end(), S.begin(), S.end());
    sort(N.begin(),N.end());
    double Px[N.size()], Py[N.size()];
    for (int i: PD) {
        Px[i] = get_random_number();
        Py[i] = get_random_number();
    }
    t_ij = new double*[N.size()];
    for (int i = 0; i < N.size(); i++)
        t_ij[i] = new double[N.size()];
    al_i = new double[N.size()];
    be_i = new double[N.size()];
    ga_i = new double[N.size()];
    double max_dist = -1.0;
    for (int i: N) {
        for (int j: N) {
            t_ij[i][j] = distance(Px, Py, i, j);
        }
        double i_max_dist = *max_element(t_ij[i], t_ij[i] + N.size());
        if (max_dist < i_max_dist) {
            max_dist = i_max_dist;
        }
        al_i[i] = 0.0;
        be_i[i] = DBL_MAX;
        ga_i[i] = 1.0;
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
}
Parameter::~Parameter() {
    vector<int*> ipV = {h_k, n_k};
    for(auto p: ipV) {
        delete p;
    }
    vector<double*> dpV = {r_k, v_k, w_k,
                            al_i, be_i, ga_i};
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
    F.clear();
}


void def_dvs(GRBModel *model,
                GRBVar **x_ij, GRBVar *u_i,
                const Parameter &prmt){
    char buf[DEFAULT_BUFFER_SIZE];
    for (int i: prmt.N) {
        for (int j: prmt.N) {
            sprintf(buf, "x[%d][%d]", i, j);
            x_ij[i][j] = (*model).addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, buf);
        }
        sprintf(buf, "u[%d]", i);
        u_i[i] = (*model).addVar(0.0, DBL_MAX, 0.0, GRB_CONTINUOUS, buf);
    }
}

void def_FC_cnsts(GRBModel *model, GRBVar **x_ij,
                  const Parameter &prmt) {
    char buf[DEFAULT_BUFFER_SIZE];
    GRBLinExpr linExpr;
    linExpr = 0;
    sprintf(buf, "iF1");
    for (int j: prmt.N) {
        linExpr += x_ij[prmt.o][j];
    }
    (*model).addConstr(linExpr == 1, buf);
    linExpr = 0;
    sprintf(buf, "iF2");
    for (int j: prmt.N) {
        linExpr += x_ij[j][prmt.d];
    }
    (*model).addConstr(linExpr == 1, buf);
    //
    for (int i: prmt.S) {
        if (i == prmt.o || i == prmt.d)
            continue;
        linExpr = 0;
        sprintf(buf, "iFS1[%d]", i);
        for (int j: prmt.N) {
            linExpr += x_ij[i][j];
        }
        (*model).addConstr(linExpr == 1, buf);
        linExpr = 0;
        sprintf(buf, "iFS2[%d]", i);
        for (int j: prmt.N) {
            linExpr += x_ij[j][i];
        }
        (*model).addConstr(linExpr == 1, buf);
    }
    //
    linExpr = 0;
    sprintf(buf, "xF1");
    for (int j: prmt.N) {
        linExpr += x_ij[j][prmt.o];
    }
    (*model).addConstr(linExpr == 0, buf);
    linExpr = 0;
    sprintf(buf, "xF2");
    for (int j: prmt.N) {
        linExpr += x_ij[prmt.d][j];
    }
    (*model).addConstr(linExpr == 0, buf);
    //
    for (int k: prmt.K) {
        if (prmt.F.find(k) != prmt.F.end())
            continue;
        linExpr = 0;
        sprintf(buf, "xFN[%d]", k);
        for (int j: prmt.N) {
            linExpr += x_ij[j][prmt.n_k[k]];
        }
        (*model).addConstr(linExpr == 0, buf);
    }
    //
    for (int k: prmt.F) {
        linExpr = 0;
        sprintf(buf, "tFC[%d]", k);
        for (int j: prmt.N) {
            linExpr += x_ij[prmt.n_k[k]][j];
        }
        for (int j: prmt.N) {
            linExpr -= x_ij[j][prmt.h_k[k]];
        }
        (*model).addConstr(linExpr <= 0, buf);
    }
    //
    for (int i: prmt.PD) {
        linExpr = 0;
        sprintf(buf, "FC[%d]", i);
        for (int j: prmt.N) {
            linExpr += x_ij[i][j];
        }
        for (int j: prmt.N) {
            linExpr -= x_ij[j][i];
        }
        (*model).addConstr(linExpr == 0, buf);
    }
}

void def_AT_cnsts(GRBModel *model, GRBVar **x_ij, GRBVar *u_i,
                  const Parameter &prmt) {
    char buf[DEFAULT_BUFFER_SIZE];
    for (int i: prmt.N) {
        sprintf(buf, "TW1[%d]", i);
        (*model).addConstr(prmt.al_i[i] <= u_i[i], buf);
        sprintf(buf, "TW2[%d]", i);
        (*model).addConstr(u_i[i] <= prmt.be_i[i], buf);
    }
    //
    for (int k: prmt.F) {
        sprintf(buf, "WD_S[%d]", k);
        (*model).addConstr(u_i[prmt.h_k[k]] <= u_i[prmt.n_k[k]], buf);
    }
    //
    for (int i: prmt.S) {
        for (int j: prmt.S) {
            sprintf(buf, "RR_P[%d,%d]", i, j);
            (*model).addConstr(prmt.c_ij[i][j] * u_i[i] <= u_i[j], buf);
        }
    }
    //
    for (int i: prmt.N) {
        for (int j: prmt.N) {
            sprintf(buf, "AT[%d,%d]", i, j);
            (*model).addConstr(u_i[i] + prmt.ga_i[i] + prmt.t_ij[i][j] <= u_i[j] + prmt.M * (1 - x_ij[i][j]), buf);
        }
    }
//    //
    GRBLinExpr linExpr = 0;
    sprintf(buf, "DT");
    for (int i: prmt.N) {
        for (int j: prmt.N) {
            linExpr += prmt.t_ij[i][j] * x_ij[i][j];
        }
    }
    (*model).addConstr(linExpr <= prmt.ul, buf);
}

void def_CP_cnsts(GRBModel *model, GRBVar **x_ij,
                  const Parameter &prmt) {
    char buf[DEFAULT_BUFFER_SIZE];
    GRBLinExpr linExpr;
    linExpr = 0;
    sprintf(buf, "bv");
    for (int k: prmt.K) {
        for (int j: prmt.PD) {
            linExpr += prmt.v_k[k] * x_ij[j][prmt.n_k[k]];
        }
    }
    (*model).addConstr(linExpr <= prmt.bv, buf);
    //
    linExpr = 0;
    sprintf(buf, "bw");
    for (int k: prmt.K) {
        for (int j: prmt.PD) {
            linExpr += prmt.w_k[k] * x_ij[j][prmt.n_k[k]];
        }
    }
    (*model).addConstr(linExpr <= prmt.bw, buf);
}

void get_x_ij(double **_x_ij, GRBVar **x_ij,
                  const Parameter &prmt) {
    for (int i: prmt.N) {
        for (int j: prmt.N) {
            _x_ij[i][j] = x_ij[i][j].get(GRB_DoubleAttr_X);
        }
    }
}

set<arc> get_delta_S(set<int> S, const Parameter &prmt) {
    set<arc> delta_S_P;
    for (int i: S) {
        for (int j: prmt.N) {
            if (S.find(j) != S.end())
                continue;
            arc ij{i, j};
            delta_S_P.insert(ij);
        }
    }
    set<arc> delta_S_M;
    for (int i: prmt.N) {
        if (S.find(i) != S.end())
            continue;
        for (int j: S) {
            arc ij{i, j};
            delta_S_M.insert(ij);
        }
    }
    //
    set<arc> delta_S(delta_S_P);
    delta_S.insert(delta_S_M.begin(), delta_S_M.end());
    return delta_S;
}

set<int> get_sigma_S(set<int> S, const Parameter &prmt) {
    set<int> sigma_S;
    for (int i: S) {
        if (prmt.P.find(i) == prmt.P.end())
            continue;
        sigma_S.insert(i + (int) prmt.K.size());
    }
    return sigma_S;
}

set<int> get_pi_S(set<int> S, const Parameter &prmt) {
    set<int> pi_S;
    for (int i: S) {
        if (prmt.D.find(i) == prmt.D.end())
            continue;
        pi_S.insert(i - (int) prmt.K.size());
    }
    return pi_S;
}

set<int> get_bar_S(set<int> S, const Parameter &prmt) {
    set<int> bar_S;
    for (int i: prmt.N) {
        if (S.find(i) != S.end())
            continue;
        bar_S.insert(i);
    }
    return bar_S;
}

template<typename T1, typename T2>
void get_SEC1_LHS(T1 *lhs, T2 **x_ij,
                  const set<int> &S,
                  const Parameter &prmt) {
    set<arc> delta_S = get_delta_S(S, prmt);
    set<int> bar_S = get_bar_S(S, prmt);
    set<int> sigma_S = get_sigma_S(S, prmt);
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
     * -2 * \sum_{i \in \bar{S} \setminus sigma_S} \sum_{j \in S \cap \sigma(S)} x_{i, j}
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


void get_SEC1_aS(const set<int> &tabuS,
                     const set<int> &S0,
                     int *min_i, double *min_lhs,
                     double **_x_ij, const Parameter &prmt) {
    // Get a new subset by adding a new node
    for(int i: prmt.PD) {
        if(tabuS.find(i) != tabuS.end())
            continue;
        if (S0.find(i) != S0.end())
            continue;
        set<int> S1(S0);
        S1.insert(i);
        double lhs = 0.0;
        get_SEC1_LHS(&lhs, _x_ij, S1, prmt);
        if (lhs < *min_lhs) {
            *min_lhs = lhs;
            *min_i = i;
        }
        S1.clear();
    }
}

void get_SEC1_dS(const set<int> &tabuS,
                         const set<int> &S0,
                         int *min_i, double *min_lhs,
                         double **_x_ij, const Parameter &prmt) {
    // Get a new subset by removing the existing node
    for(int i: S0) {
        set<int> S1(S0);
        S1.erase(i);
        double lhs = 0.0;
        get_SEC1_LHS(&lhs, _x_ij, S1, prmt);
        if (lhs < *min_lhs) {
            *min_lhs = lhs;
            *min_i = i;
        }
        S1.clear();
    }
}

vector<set<int>> get_SEC1(double **_x_ij, const Parameter &prmt) {
    vector<set<int>> violated_subsets;
    //
    set<int> S0;
    set<int> tabuS;
    auto randIt = prmt.PD.begin();
    advance(randIt, rand() % prmt.PD.size());
    S0.insert(*randIt);
    int min_i;
    double min_lhs;
    for (int i = 0; i < 25; i++) {
        min_i = -1;
        min_lhs = DBL_MAX;
        set<int> S1(S0);
        if (S0.size() < 2) {
            get_SEC1_aS(tabuS, S0, &min_i, &min_lhs, _x_ij, prmt);
            S1.insert(min_i);
        } else if (get_random_number() < 0.5) {
            get_SEC1_aS(tabuS, S0, &min_i, &min_lhs, _x_ij, prmt);
            S1.insert(min_i);
        } else {
            get_SEC1_dS(tabuS, S0, &min_i, &min_lhs, _x_ij, prmt);
            S1.erase(min_i);
            tabuS.insert(min_i);
        }
        if (min_i == -1) {
            continue;
        }
        printf("%d: size %d, min_i: %d, min_lhs: %.2f\n", i, S1.size(), min_i, min_lhs);
        if (min_lhs < 2.0) {
            violated_subsets.push_back(S1);
        }
        S0.clear();
        S0 = S1;
    }
    return violated_subsets;
}

void run_example_PDPTW() {
    int num_tasks, num_rr_points, max_reward;
    double bv, bw, ul;
    ex1(&num_tasks, &num_rr_points, &max_reward, &bv, &bw, &ul);
    const Parameter prmt(num_tasks, num_rr_points, max_reward,
                         bv, bw, ul);
    /*
     * Declare env. and decision variables
     */
    GRBEnv *env = NULL;
    GRBVar **x_ij = NULL, *u_i = NULL;
    double **_x_ij;
    x_ij = new GRBVar *[prmt.N.size()];
    _x_ij = new double *[prmt.N.size()];
    for (int i: prmt.N) {
        x_ij[i] = new GRBVar[prmt.N.size()];
        _x_ij[i] = new double[prmt.N.size()];
    }
    u_i = new GRBVar[prmt.N.size()];
    try {
        env = new GRBEnv();
        GRBModel model = GRBModel(*env);
        def_dvs(&model, x_ij, u_i, prmt);
        /*
         * Set a objective function
         */
        GRBLinExpr objF = 0;
        for (int k: prmt.K) {
            for (int j: prmt.PD) {
                objF += prmt.r_k[k] * x_ij[j][prmt.n_k[k]];
            }
        }
        model.setObjective(objF, GRB_MAXIMIZE);
        //
        def_FC_cnsts(&model, x_ij, prmt);
        def_AT_cnsts(&model, x_ij, u_i, prmt);
        def_CP_cnsts(&model, x_ij, prmt);
        model.update();
        model.write("/Users/ckhan/workspace/BNC/BNC/test_LP.lp");
        //
        model.optimize();
        get_x_ij(_x_ij, x_ij, prmt);
        for (int i: prmt.N) {
            for (int j: prmt.N) {
                printf("x[%d][%d]: %.2f\n", i, j, _x_ij[i][j]);
            }
        }
        vector<set<int>> SEC1 = get_SEC1(_x_ij, prmt);
        for (set<int> S: SEC1) {
            GRBLinExpr lhs = 0.0;
            get_SEC1_LHS(&lhs, x_ij, SEC1[0], prmt);
            model.addConstr(lhs >= 2);
        }
        model.write("/Users/ckhan/workspace/BNC/BNC/test_LP_cut.lp");
        model.optimize();
        
        printf("test");
        
//        cout << "!!!!!Obj: " << rmodel.get(GRB_DoubleAttr_ObjVal) << endl;
//        cout << x_ij[0][0].get(GRB_StringAttr_VarName) << " " <<
//        x_ij[0][0].get(GRB_DoubleAttr_X) << endl;
//        GRBModel rmodel = model.relax();
//        rmodel.write("/Users/ckhan/workspace/BNC/BNC/test_LP1.lp");
//        rmodel.optimize();
//        if (rmodel.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
//            rmodel.computeIIS();
//            rmodel.write("/Users/ckhan/workspace/BNC/BNC/test.ilp");
//        }
//        printf("Solve a new model!!\n");
//        rmodel.addConstr(x_ij[0][1] == 1, "testCnst");
//        rmodel.write("/Users/ckhan/workspace/BNC/BNC/test_LP2.lp");
//        rmodel.optimize();
//        GRBConstr cnst = rmodel.getConstrByName("bv");
//        double a = cnst.get(GRB_DoubleAttr_Pi);
//        printf("%.2f", a);
        
    } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
    }

    for (int i = 0; i < prmt.N.size(); i++)
        delete [] x_ij[i];
    delete [] x_ij;
    delete [] u_i;
    delete env;

    printf("ended!");
}
