//
//  MathematicalModel.cpp
//  BNC
//
//  Created by Chung-Kyun HAN on 1/8/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "MathematicalModel.hpp"


MathematicalModel::MathematicalModel(char xType, Problem *prob) {
    this->prob = prob;
    x_ij = new GRBVar *[(*prob).N.size()];
    for (int i: (*prob).N) {
        x_ij[i] = new GRBVar[(*prob).N.size()];
    }
    u_i = new GRBVar[(*prob).N.size()];
    env = new GRBEnv("");
    grbModel = new GRBModel(*env);
    def_dvs(xType);
    def_FC_cnsts();
    def_AT_cnsts();
    def_CP_cnsts();
    def_objF();
    grbModel->set(GRB_IntParam_LogToConsole, 0);
}

MathematicalModel* MathematicalModel::clone() {
    char buf[DEFAULT_BUFFER_SIZE];
    MathematicalModel *mm = new MathematicalModel();
    mm->prob = prob;
    (*grbModel).update();
    mm->grbModel = new GRBModel(*grbModel);
    mm->x_ij = new GRBVar *[(*prob).N.size()];
    mm->u_i = new GRBVar[(*prob).N.size()];
    for (int i: (*prob).N) {
        mm->x_ij[i] = new GRBVar[(*prob).N.size()];
        for (int j: (*prob).N) {
            sprintf(buf, "x[%d][%d]", i, j);
            mm->x_ij[i][j] = (*(*mm).grbModel).getVarByName(buf);
        }
        sprintf(buf, "u[%d]", i);
        mm->u_i[i] = (*(*mm).grbModel).getVarByName(buf);
    }
    return mm;
}

void MathematicalModel::add_intConstr(int i, int j, int rhs) {
    char buf[DEFAULT_BUFFER_SIZE];
    sprintf(buf, "intConst[%d][%d]", i, j);
    if (rhs == 0) {
        (*grbModel).addConstr(x_ij[i][j] <= rhs, buf);
    } else {
        assert(rhs == 1);
        (*grbModel).addConstr(x_ij[i][j] >= rhs, buf);
    }
    (*grbModel).update();
}

void MathematicalModel::def_dvs(char modType){
    char buf[DEFAULT_BUFFER_SIZE];
    bool isIntegerModel = modType == 'I' ? true : false;
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            sprintf(buf, "x[%d][%d]", i, j);
            if (isIntegerModel) {
                x_ij[i][j] = (*grbModel).addVar(0.0, 1.0, 0.0, GRB_BINARY, buf);
            } else {
                x_ij[i][j] = (*grbModel).addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, buf);
            }
        }
        sprintf(buf, "u[%d]", i);
        u_i[i] = (*grbModel).addVar(0.0, DBL_MAX, 0.0, GRB_CONTINUOUS, buf);
    }
}

void MathematicalModel::def_FC_cnsts() {
    char buf[DEFAULT_BUFFER_SIZE];
    GRBLinExpr linExpr;
    /*
     *  \sum_{ j \in N } x_{j, o} = \sum_{ j \in N } x_{d, j} = 0 \label{eq:xF}
     */
    linExpr = 0;
    sprintf(buf, "xF0");
    for (int i: (*prob).N) {
        linExpr += x_ij[i][i];
    }
    (*grbModel).addConstr(linExpr == 0, buf);
    linExpr = 0;
    sprintf(buf, "xF1");
    for (int j: (*prob).N) {
        linExpr += x_ij[j][(*prob).o];
    }
    (*grbModel).addConstr(linExpr == 0, buf);
    linExpr = 0;
    sprintf(buf, "xF2");
    for (int j: (*prob).N) {
        linExpr += x_ij[(*prob).d][j];
    }
    (*grbModel).addConstr(linExpr == 0, buf);
    /*
     *  \sum_{ j \in N } x_{o, j} = \sum_{ j \in N } x_{j, d} = 1 \label{eq:iF}
     */
    linExpr = 0;
    sprintf(buf, "iF1");
    for (int j: (*prob).N) {
        linExpr += x_ij[(*prob).o][j];
    }
    (*grbModel).addConstr(linExpr == 1, buf);
    linExpr = 0;
    sprintf(buf, "iF2");
    for (int j: (*prob).N) {
        linExpr += x_ij[j][(*prob).d];
    }
    (*grbModel).addConstr(linExpr == 1, buf);
    /*
     *  \sum_{ i \neq j \in N } x_{i, j} = \sum_{ i \neq j \in N } x_{j, i} = 1, \forall i \in S \setminus \{ o, d \} \label{eq:iFS}
     */
    for (int i: (*prob).S) {
        if (i == (*prob).o || i == (*prob).d)
            continue;
        linExpr = 0;
        sprintf(buf, "iFS1[%d]", i);
        for (int j: (*prob).N) {
            if (i == j)
                continue;
            linExpr += x_ij[i][j];
        }
        (*grbModel).addConstr(linExpr == 1, buf);
        linExpr = 0;
        sprintf(buf, "iFS2[%d]", i);
        for (int j: (*prob).N) {
            if (i == j)
                continue;
            linExpr += x_ij[j] [i];
        }
        (*grbModel).addConstr(linExpr == 1, buf);
    }
    /*
     *   \sum_{ j \in N } x_{n_{k}, j} \le \sum_{ j \in N } x_{j, h_{k}}, \forall k \in K \label{eq:wV}
     */
    for (int k: (*prob).K) {
        linExpr = 0;
        sprintf(buf, "tFC[%d]", k);
        for (int j: (*prob).N) {
            linExpr += x_ij[(*prob).n_k[k]][j];
        }
        for (int j: (*prob).N) {
            linExpr -= x_ij[j][(*prob).h_k[k]];
        }
        (*grbModel).addConstr(linExpr <= 0, buf);
    }
    /*
     *  \sum_{ j \in N } x_{i, j} = \sum_{ j \in N } x_{j, i},~~~~\forall i \in P \cup D \label{eq:FC}
     */
    for (int i: (*prob).PD) {
        linExpr = 0;
        for (int j: (*prob).N) {
            linExpr += x_ij[i][j];
        }
        sprintf(buf, "FC_1[%d]", i);
        (*grbModel).addConstr(linExpr <= 1, buf);
        for (int j: (*prob).N) {
            linExpr -= x_ij[j][i];
        }
        sprintf(buf, "FC[%d]", i);
        (*grbModel).addConstr(linExpr == 0, buf);
    }
}

void MathematicalModel::def_AT_cnsts() {
    char buf[DEFAULT_BUFFER_SIZE];
    /*
     *  \alpha_{i} \le \mu_{i} \le \beta_{i},~~~~\forall i \in N \label{eq:TW}
     */
    for (int i: (*prob).N) {
        sprintf(buf, "TW1[%d]", i);
        (*grbModel).addConstr((*prob).al_i[i] <= u_i[i], buf);
        sprintf(buf, "TW2[%d]", i);
        (*grbModel).addConstr(u_i[i] <= (*prob).be_i[i], buf);
    }
    /*
     *  \mu_{h_{k}} \le \mu_{n_{k}},~~~~\forall k \in K, \label{eq:WD_S}
     */
    for (int k: (*prob).K) {
        sprintf(buf, "WD_S[%d]", k);
        (*grbModel).addConstr(u_i[(*prob).h_k[k]] <= u_i[(*prob).n_k[k]], buf);
    }
    /*
     *  c_{i, j} \mu_{i} \le \mu_{j},~~~~\forall i, j \in S, \label{eq:RR_P}
     */
    for (int i: (*prob).S) {
        for (int j: (*prob).S) {
            sprintf(buf, "RR_P[%d,%d]", i, j);
            (*grbModel).addConstr((*prob).c_ij[i][j] * u_i[i] <= u_i[j], buf);
        }
    }
    /*
     *  \mu_{i} + t_{i, j} \le \mu_{j} + M (1 - x_{i, j}),~~~~\forall i, j \in N \label{eq:AT}
     */
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            sprintf(buf, "AT[%d,%d]", i, j);
            (*grbModel).addConstr(u_i[i] + (*prob).t_ij[i][j] <= u_i[j] + (*prob).M * (1 - x_ij[i][j]), buf);
        }
    }
    /*
     *  \sum_{ i, j \in N} t_{i, j} x_{i, j}  \le \textbf{u} \label{eq:DT}
     */
    GRBLinExpr linExpr = 0;
    sprintf(buf, "DT");
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            linExpr += (*prob).t_ij[i][j] * x_ij[i][j];
        }
    }
    (*grbModel).addConstr(linExpr <= (*prob).bu, buf);
}


void MathematicalModel::def_CP_cnsts() {
    char buf[DEFAULT_BUFFER_SIZE];
    GRBLinExpr linExpr;
    linExpr = 0;
    sprintf(buf, "bv");
    for (int k: (*prob).K) {
        for (int j: (*prob).PD) {
            linExpr += (*prob).v_k[k] * x_ij[j][(*prob).n_k[k]];
        }
    }
    (*grbModel).addConstr(linExpr <= (*prob).bv, buf);
    //
    linExpr = 0;
    sprintf(buf, "bw");
    for (int k: (*prob).K) {
        for (int j: (*prob).PD) {
            linExpr += (*prob).w_k[k] * x_ij[j][(*prob).n_k[k]];
        }
    }
    (*grbModel).addConstr(linExpr <= (*prob).bw, buf);
}

void MathematicalModel::def_objF() {
    GRBLinExpr objF = 0;
    for (int k: (*prob).K) {
        for (int j: (*prob).PD) {
            objF += (*prob).r_k[k] * x_ij[j][(*prob).n_k[k]];
        }
    }
    (*grbModel).setObjective(objF, GRB_MAXIMIZE);
}

void MathematicalModel::get_x_ij(double **_x_ij) {
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            _x_ij[i][j] = x_ij[i][j].get(GRB_DoubleAttr_X);
        }
    }
}

void MathematicalModel::get_u_i(double *_u_i) {
    for (int i: (*prob).N) {
        _u_i[i] = u_i[i].get(GRB_DoubleAttr_X);
    }
}

MathematicalModel::~MathematicalModel() {
    for (int i = 0; i < (*prob).N.size(); i++)
        delete [] x_ij[i];
    delete [] x_ij;
    delete [] u_i;
    delete grbModel;
}

void subtourelim::callback() {
    try {
        if (where == GRB_CB_MIP) {
            // General MIP callback
            double nodecnt = getDoubleInfo(GRB_CB_MIP_NODCNT);
            double objbst = getDoubleInfo(GRB_CB_MIP_OBJBST);
            double objbnd = getDoubleInfo(GRB_CB_MIP_OBJBND);
            int cutcnt = getIntInfo(GRB_CB_MIP_CUTCNT);
            int solcnt = getIntInfo(GRB_CB_MIP_SOLCNT);
//            if (nodecnt - lastnode >= 100) {
//                lastnode = nodecnt;
//                int actnodes = (int) getDoubleInfo(GRB_CB_MIP_NODLFT);
//                int itcnt = (int) getDoubleInfo(GRB_CB_MIP_ITRCNT);
//                int cutcnt = getIntInfo(GRB_CB_MIP_CUTCNT);
//                cout << nodecnt << " " << actnodes << " " << itcnt
//                << " " << objbst << " " << objbnd << " "
//                << solcnt << " " << cutcnt << endl;
            }
        if (where == GRB_CB_MIPNODE) {
            // Found an integer feasible solution - does it visit every node?
            double **x = new double*[n];
            int *tour = new int[n];
            int i, j, len;
            for (i = 0; i < n; i++)
                x[i] = getSolution(vars[i], n);
            
            //                findsubtour(n, x, &len, tour);
            if (len < n) {
                // Add subtour elimination constraint
                GRBLinExpr expr = 0;
                for (i = 0; i < len; i++)
                    for (j = i+1; j < len; j++)
                        expr += vars[tour[i]][tour[j]];
                addLazy(expr <= len-1);
            }
            
            for (i = 0; i < n; i++)
                delete[] x[i];
            delete[] x;
            delete[] tour;
        }
    } catch (GRBException e) {
        //            cout << "Error number: " << e.getErrorCode() << endl;
        //            cout << e.getMessage() << endl;
    } catch (...) {
        //            cout << "Error during callback" << endl;
    }
}
