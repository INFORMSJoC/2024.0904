#pragma once
#include<vector>
#include<tuple>
#include "ReadIns.h"
#include<chrono>
#include "gurobi_c++.h"

// -------------- the information of subproblem
class InformationSolveSP
{
public:
    int SP_status;
    double SP_gap;
    double SP_total_time;
    double SP_OptObjVal;
    std::vector<int> O;
    std::vector<int> E;
    std::vector<int> A;

    InformationSolveSP() : SP_status(-1), SP_gap(-1), SP_total_time(0.0), SP_OptObjVal(-1) {}
};



// -------------- the information of the final result
class LBBD_BS_Info
{
public:
    int LBBD_iter;
    int MP_status;
    int SP_status;
    int OP_status;
    double UB;
    double LB;
    double obj_MP;
    double obj_SP;
    double obj_mp_theta;
    double obj_mp_f1f2;

    double total_time_SP;
    double total_time_MP;
    double Total_time;

    InformationSolveSP* lastSP;

    double LB_delta_obj;    // get after obtain delta
    double LB_delta_time;
    // std::vector<int> LB_Config_Seq;

    // final optimal solution
    std::vector<int> lB;
    std::vector<int> lC;
    std::vector<int> lO;
    std::vector<int> lE;
    std::vector<int> lA;

    int lbbd_cut_opt;
    int lbbd_cut_nogood;

    double lbbd_gap;
    double sp_gap;
    double mp_gap;


    LBBD_BS_Info() : LBBD_iter(0), MP_status(-1), SP_status(-1), OP_status(-1), UB(100000), LB(1e-3), obj_MP(-1), obj_SP(-1),
        obj_mp_theta(-1), obj_mp_f1f2(-1),
        total_time_SP(0.0),
        total_time_MP(0.0), Total_time(0.0),
        lbbd_cut_opt(-1), lbbd_cut_nogood(0), lbbd_gap(-1), sp_gap(-1), mp_gap(-1) {}
};



void Solve_SP(InformationSolveSP& spSol_Info, Instance ins, std::map< std::tuple<int, int>,
    std::vector<int> > Vbm, GRBEnv& env, double time_limit)
{
    auto sp_t1 = std::chrono::high_resolution_clock::now();

    // build SP model
    GRBModel SP_model = GRBModel(env);
    //SP_model.set(GRB_IntParam_OutputFlag, 0);
    SP_model.set(GRB_DoubleParam_TimeLimit, time_limit);

    // add variables
    std::vector< std::vector<GRBVar> > U(ins.d_T);     //U[p][w]: for linearization
    for (int p = 0; p < ins.V.size(); p++)
    {
        std::vector<GRBVar> row(ins.W.size());
        for (int w = 0; w < ins.W.size(); w++)
        {
            std::string var_name = "u[" + std::to_string(ins.V[p]) + "][" + std::to_string((int)ins.W[w]) + "]";
            GRBVar u = SP_model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, var_name);
            row[w] = u;
        }
        U[p] = row;
    }

    std::map<std::tuple<int, int, int>, GRBVar> Y;
    for (int i = 1; i <= ins.d_T; i++)
        for (int p = 1; p <= ins.d_T; p++)
            for (int j = 1; j <= ins.J.size(); j++)
            {
                std::string var_name = "y[" + std::to_string(i) + "][" + std::to_string(p) + "][" + std::to_string(j) + "]";
                GRBVar y = SP_model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                Y[std::make_tuple(i, p, j)] = y;
            }


    std::map<std::tuple<int, int>, GRBVar> H;	// h
    for (int i = 0; i < ins.d_T; i++)
        for (int l = 0; l < ins.L.size(); l++)
        {
            std::string var_name = "h[" + std::to_string(ins.V[i]) + "][" + std::to_string(ins.L[l]) + "]";
            GRBVar h = SP_model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
            H[std::make_tuple(ins.V[i], ins.L[l])] = h;
        }



    std::map<std::tuple<int, int, int>, GRBVar> Z;
    for (int i = 1; i <= ins.d_T; i++)
        for (int v = 1; v <= i; v++)
            for (int l = 1; l <= ins.L.size(); l++)
            {
                std::string var_name = "z[" + std::to_string(i) + "][" + std::to_string(v) + "][" + std::to_string(l) + "]";
                GRBVar z = SP_model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                Z[std::make_tuple(i, v, l)] = z;
            }

    // add constraints
    for (int p = 1; p <= ins.d_T; p++)
        for (int w = 1; w <= ins.W.size(); w++)
        {
            GRBLinExpr con1 = 0;
            for (int j = 1; j <= ins.J.size(); j++)
                for (int p_ = 1; p_ <= p; p_++)
                    for (int i = 1; i <= ins.d_T; i++)
                        con1 += ins.r_j_w[j][w] * Y[std::make_tuple(i, p_, j)];

            SP_model.addConstr(U[p - 1][w - 1] - con1 + p * ins.sigma_w[w], GRB_GREATER_EQUAL, 0.0);
            SP_model.addConstr(U[p - 1][w - 1] + con1 - p * ins.sigma_w[w], GRB_GREATER_EQUAL, 0.0);
        }

    for (int i = 1; i <= ins.d_T; i++)
    {
        GRBLinExpr con2 = 0;
        for (int l = 1; l <= ins.L.size(); l++)
            con2 += H[std::make_tuple(i, l)];
        SP_model.addConstr(con2, GRB_EQUAL, 1.0);
    }

    for (int v = 2; v <= ins.d_T; v++)
        for (int l = 1; l <= ins.L.size(); l++)
            for (int k = 1; k < v; k++)
                for (int i = v; i <= ins.d_T; i++)
                {
                    GRBLinExpr con3 = H[std::make_tuple(k, l)] + H[std::make_tuple(v, l)] - 2;
                    for (int l_ = 1; l_ <= ins.L.size(); l_++)
                        con3 += Z[std::make_tuple(i, v, l_)];

                    for (int i_ = k; i_ < i; i_++)
                        for (int l_ = 1; l_ <= ins.L.size(); l_++)
                            con3 -= Z[std::make_tuple(i_, k, l_)];

                    SP_model.addConstr(con3, GRB_LESS_EQUAL, 0.0);
                }

    for (int i = 1; i <= ins.d_T; i++)
    {
        GRBLinExpr con4 = 0;
        for (int v = 1; v <= i; v++)
            for (int l = 1; l <= ins.L.size(); l++)
                con4 += Z[std::make_tuple(i, v, l)];
        SP_model.addConstr(con4, GRB_LESS_EQUAL, 1.0);
    }

    for (int v = 1; v <= ins.d_T; v++)
    {
        GRBLinExpr con5 = 0;
        for (int i = v; i <= ins.d_T; i++)
            for (int l = 1; l <= ins.L.size(); l++)
                con5 += Z[std::make_tuple(i, v, l)];
        SP_model.addConstr(con5, GRB_LESS_EQUAL, 1.0);
    }


    for (int i = 1; i <= ins.d_T; i++)
        for (int l = 1; l <= ins.L.size(); l++)
        {
            GRBLinExpr con6 = 0;
            for (int i_ = 1; i_ <= i; i_++)
                con6 += H[std::make_tuple(i_, l)];
            for (int i_ = 1; i_ <= i; i_++)
                for (int v = 1; v <= i_; v++)
                    con6 -= Z[std::make_tuple(i_, v, l)];
            SP_model.addConstr(con6, GRB_LESS_EQUAL, ins.q0);
        }

    for (int v = 1; v <= ins.d_T; v++)
        for (int l = 1; l <= ins.L.size(); l++)
        {
            GRBLinExpr con7 = 0;
            for (int i = v; i <= ins.d_T; i++)
                con7 += Z[std::make_tuple(i, v, l)];
            con7 -= H[std::make_tuple(v, l)];
            SP_model.addConstr(con7, GRB_LESS_EQUAL, 0.0);
        }


    for (int i = 1; i <= ins.d_T; i++)
    {
        GRBLinExpr con8 = 0;
        for (int j = 1; j <= ins.J.size(); j++)
            for (int p = 1; p <= ins.d_T; p++)
                con8 += Y[std::make_tuple(i, p, j)];
        SP_model.addConstr(con8, GRB_EQUAL, 1.0);
    }

    for (int p = 1; p <= ins.d_T; p++)
    {
        GRBLinExpr con9 = 0;
        for (int i = 1; i <= ins.d_T; i++)
            for (int j = 1; j <= ins.J.size(); j++)
                con9 += Y[std::make_tuple(i, p, j)];
        SP_model.addConstr(con9, GRB_EQUAL, 1.0);
    }


    for (int i = 0; i < ins.d_bmj.size(); i++)
    {
        GRBLinExpr con10 = 0;
        for (int v = 0; v < Vbm[std::make_tuple(ins.d_bmj[i].body, ins.d_bmj[i].color)].size(); v++)
            for (int p = 1; p <= ins.d_T; p++)
                con10 += Y[std::make_tuple(Vbm[std::make_tuple(ins.d_bmj[i].body, ins.d_bmj[i].color)][v], p, ins.d_bmj[i].config)];

        SP_model.addConstr(con10, GRB_EQUAL, ins.d_bmj[i].demand);
    }

    for (int i1 = 1; i1 < ins.d_T; i1++)
        for (int p1 = 1; p1 <= ins.d_T; p1++)
            for (int i2 = i1 + 1; i2 <= ins.d_T; i2++)
                for (int l = 1; l <= ins.L.size(); l++)
                {
                    GRBLinExpr con11 = H[std::make_tuple(i1, l)] + H[std::make_tuple(i2, l)] - 3;
                    for (int p = p1; p <= ins.d_T; p++)
                        for (int j = 1; j <= ins.J.size(); j++)
                            con11 += Y[std::make_tuple(i1, p, j)];
                    for (int p = 1; p <= p1; p++)
                        for (int j = 1; j <= ins.J.size(); j++)
                            con11 += Y[std::make_tuple(i2, p, j)];
                    SP_model.addConstr(con11, GRB_LESS_EQUAL, 0.0);
                }


    for (int i = 1; i <= ins.d_T; i++)
    {
        GRBLinExpr con12 = 0;
        for (int p = 1; p <= ins.d_T; p++)
            for (int j = 1; j <= ins.J.size(); j++)
                con12 += p * Y[std::make_tuple(i, p, j)];
        con12 -= i;
        for (int l = 1; l <= ins.L.size(); l++)
            con12 += i * Z[std::make_tuple(i, i, l)];
        con12 -= 1;
        SP_model.addConstr(con12, GRB_GREATER_EQUAL, 0.0);
    }


    for (int i = 2; i <= ins.d_T; i++)
    {
        GRBLinExpr con13 = 0;
        GRBLinExpr con14 = 0;
        for (int p = 1; p <= ins.d_T; p++)
            for (int j = 1; j <= ins.J.size(); j++)
            {
                con13 += p * Y[std::make_tuple(i, p, j)];
                con14 += p * Y[std::make_tuple(i, p, j)];
            }
        con13 -= 1;
        con14 -= 1;

        for (int i_ = 1; i_ < i; i_++)
            for (int v = 1; v <= i_; v++)
                for (int l = 1; l <= ins.L.size(); l++)
                {
                    con13 -= Z[std::make_tuple(i_, v, l)];
                    con14 -= Z[std::make_tuple(i_, v, l)];
                }
        con13 += ins.d_T;
        con14 -= ins.d_T;

        for (int l = 1; l <= ins.L.size(); l++)
        {
            con13 -= ins.d_T * Z[std::make_tuple(i, i, l)];
            con14 += ins.d_T * Z[std::make_tuple(i, i, l)];
        }

        SP_model.addConstr(con13, GRB_GREATER_EQUAL, 0.0);
        SP_model.addConstr(con14, GRB_LESS_EQUAL, 0.0);
    }

    // symmetry breaking
    for (int i = 0; i <= ins.cars.size(); i++)
    {
        int tc = Vbm[std::make_tuple(ins.cars[i].body, ins.cars[i].color)].size();
        if (tc >= 2)
        {
            for (int i1 = 0; i1 < tc - 1; i1++)
                for (int i2 = i1 + 1; i2 < tc; i2++)
                {
                    int c1 = Vbm[std::make_tuple(ins.cars[i].body, ins.cars[i].color)][i1];
                    int c2 = Vbm[std::make_tuple(ins.cars[i].body, ins.cars[i].color)][i2];

                    GRBLinExpr con15 = 1;
                    for (int p = 1; p <= ins.d_T; p++)
                        for (int j = 1; j <= ins.J.size(); j++)
                        {
                            con15 += p * Y[std::make_tuple(c1, p, j)];
                            con15 -= p * Y[std::make_tuple(c2, p, j)];
                        }
                    SP_model.addConstr(con15, GRB_LESS_EQUAL, 0.0);
                }
        }



    }


    // Preprocess
    if (ins.d_T > ins.L.size() * ins.q0)
    {
        for (int i = 1 + (ins.L.size() - 1) * ins.q0 + 1; i < 1 + (ins.L.size() - 1) * ins.q0 + ins.q0; i++)
        {
            int p_i = i - (1 + (ins.L.size() - 1) * ins.q0);
            for (int p = 1; p <= p_i; p++)
            {
                GRBLinExpr sp_con5 = 0;
                for (int j = 1; j <= ins.J.size(); j++)
                    sp_con5 += Y[std::make_tuple(i, p, j)];
                SP_model.addConstr(sp_con5, GRB_EQUAL, 0.0);
            }
        }
        for (int i = 1 + ins.q0 * ins.L.size(); i <= ins.d_T; i++)
        {
            int ti = (i - ins.L.size() * ins.q0) / ins.L.size();
            int ri = (i - ins.L.size() * ins.q0) % ins.L.size();
            int pi = ins.q0 + ti * ins.L.size() + ri;
            for (int p = 1; p < pi; p++)
            {
                GRBLinExpr sp_con6 = 0;
                for (int j = 1; j <= ins.J.size(); j++)
                    sp_con6 += Y[std::make_tuple(i, p, j)];
                SP_model.addConstr(sp_con6, GRB_EQUAL, 0.0);
            }
        }
    }


    // set objective
    GRBLinExpr sp_obj = 0;
    for (int p = 1; p <= ins.d_T; p++)
        for (int w = 1; w <= ins.W.size(); w++)
            sp_obj += ins.gamma * U[p - 1][w - 1];

    SP_model.setObjective(sp_obj, GRB_MINIMIZE);

    // optimize
    SP_model.optimize();
    auto sp_t2 = std::chrono::high_resolution_clock::now();

    auto sp_time = std::chrono::duration_cast<std::chrono::milliseconds>(sp_t2 - sp_t1);


    //spSol_Info.SP_total_time = sp_time.count();
    //spSol_Info.SP_gap = SP_model.get(GRB_DoubleAttr_MIPGap);


    spSol_Info.SP_total_time = SP_model.get(GRB_DoubleAttr_Runtime);


    int SP_status;
    double SP_gap;
    double SP_total_time;
    double SP_OptObjVal;
    std::vector<int> O;
    std::vector<int> E;
    std::vector<int> A;


    if (SP_model.get(GRB_IntAttr_Status) == 2)
    {
        spSol_Info.SP_gap = 0.0;
        spSol_Info.SP_status = 1;
        spSol_Info.SP_OptObjVal = SP_model.get(GRB_DoubleAttr_ObjVal);

        for (int i = 1; i <= ins.d_T; i++)
            for (int p = 1; p <= ins.d_T; p++)
                for (int j = 1; j <= ins.J.size(); j++)
                    if (Y[std::make_tuple(i, p, j)].get(GRB_DoubleAttr_X) > 0.05)
                        spSol_Info.O.push_back(j);

        for (int p = 1; p <= ins.d_T; p++)
            for (int i = 1; i <= ins.d_T; i++)
                for (int j = 1; j <= ins.J.size(); j++)
                    if (Y[std::make_tuple(i, p, j)].get(GRB_DoubleAttr_X) > 0.05)
                        spSol_Info.E.push_back(i);
        for (int i = 1; i <= ins.d_T; i++)
            for (int l = 1; l <= ins.L.size(); l++)
                if (H[std::make_tuple(i, l)].get(GRB_DoubleAttr_X) > 0.05)
                    spSol_Info.A.push_back(l);

        std::cout << "-------------------- Get SP optimal solution.--------------" << std::endl;
    }
    else if (SP_model.get(GRB_DoubleAttr_MIPGap) < 10)
    {
        spSol_Info.SP_gap = SP_model.get(GRB_DoubleAttr_MIPGap);
        std::cout << "-------------------- Do not get SP optimal solution, have Gap--------------" << std::endl;
    }



}