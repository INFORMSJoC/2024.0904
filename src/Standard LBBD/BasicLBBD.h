#pragma once
#include "SolveSP.h"
#include "gurobi_c++.h"
#include<chrono>

// -------------- get the lower bound of the SP (delta)
void get_delta(GRBEnv& env, Instance ins, LBBD_BS_Info& lbbd)
{
    //auto t1 = std::chrono::high_resolution_clock::now();

    GRBModel model = GRBModel(env);
    model.set(GRB_IntParam_OutputFlag, 0);

    // add decision variables
    std::vector< std::vector<GRBVar> > N(ins.d_T);   //N[p][j]: the number of cars with configuration j from position 1 to p
    for (int p = 0; p < ins.V.size(); p++)
    {
        std::vector<GRBVar> row(ins.J.size());
        for (int j = 0; j < ins.J.size(); j++)
        {
            std::string var_name = "n[" + std::to_string(ins.V[p]) + "][" + std::to_string((int)ins.J[j]) + "]";
            GRBVar n = model.addVar(0.0, (double)ins.d_T, 0.0, GRB_INTEGER, var_name);
            row[j] = n;
        }
        N[p] = row;
    }

    std::vector< std::vector<GRBVar> > U(ins.d_T);     //U[p][w]: for linearization
    for (int p = 0; p < ins.V.size(); p++)
    {
        std::vector<GRBVar> row(ins.W.size());
        for (int w = 0; w < ins.W.size(); w++)
        {
            std::string var_name = "u[" + std::to_string(ins.V[p]) + "][" + std::to_string((int)ins.W[w]) + "]";
            GRBVar u = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, var_name);
            row[w] = u;
        }
        U[p] = row;
    }


    // add constraints
    for (int p = 0; p < ins.V.size(); p++)              // p_v = ins.V[p] = 1, 2, 3, 4, ...
    {
        for (int w = 0; w < ins.W.size(); w++)          // w_v = ins.W[w] = 1, 2, 3, ...
        {
            GRBLinExpr cons1 = 0;
            for (auto it = ins.r_j_w.begin(); it != ins.r_j_w.end(); it++)
                cons1 += it->second[(int)ins.W[w]] * N[p][it->first - 1];
            model.addConstr(U[p][w] - cons1 + ins.V[p] * ins.sigma_w[(int)ins.W[w]], GRB_GREATER_EQUAL, 0.0);
            model.addConstr(U[p][w] + cons1 - ins.V[p] * ins.sigma_w[(int)ins.W[w]], GRB_GREATER_EQUAL, 0.0);
        }
    }

    for (int p = 0; p < ins.V.size(); p++)
    {
        GRBLinExpr cons2 = 0;
        for (int j = 0; j < ins.J.size(); j++)
            cons2 += N[p][j];
        model.addConstr(cons2 - ins.V[p], GRB_EQUAL, 0.0);
    }

    for (int p = 0; p < ins.V.size() - 1; p++)
    {
        for (int j = 0; j < ins.J.size(); j++)
        {
            model.addConstr(N[p][j] - N[p + 1][j], GRB_LESS_EQUAL, 0.0);
            model.addConstr(N[p + 1][j] - N[p][j] - 1, GRB_LESS_EQUAL, 0.0);
        }
    }

    // set objective function
    GRBLinExpr obj = 0;
    for (int p = 0; p < ins.V.size(); p++)
        for (int w = 0; w < ins.W.size(); w++)
            obj += U[p][w];
    model.setObjective(ins.gamma * obj, GRB_MINIMIZE);

    // optimize model
    model.optimize();
    //auto t2 = std::chrono::high_resolution_clock::now();
    //auto du = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    // get the solution
    //for (int j = 0; j < ins.J.size(); j++)
    //    if (N[0][j].get(GRB_DoubleAttr_X) > 0.5)
    //    {
    //        lbbd.LB_Config_Seq.push_back(ins.J[j]);
    //        break;
    //    }
    //for (int p = 1; p < ins.V.size(); p++)
    //    for (int j = 0; j < ins.J.size(); j++)
    //        if (N[p][j].get(GRB_DoubleAttr_X) - N[p - 1][j].get(GRB_DoubleAttr_X) > 0.5)
    //        {
    //            lbbd.LB_Config_Seq.push_back(ins.J[j]);
    //            break;
    //        }

    std::cout << "LB of SP (delta): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;

    lbbd.LB_delta_obj = model.get(GRB_DoubleAttr_ObjVal);
    lbbd.LB_delta_time = model.get(GRB_DoubleAttr_Runtime);
    //lbbd.LB_delta_time = du.count();

}



// --------------- basic LBBD framework - no B&C to solve SP
void BasicLBBD(Instance ins, std::string const& ins_name, std::string const& txtfile, std::string const& csv1,
    std::string const& csv2, double T_time_limit)
{
    // ins_name: 10_1 for example
    // txtfile: output the final results
    // csv2: record each loop
    // csv1: record the final results
    // T_time_limit: time limit for LBBD loop
    // mp_time_lit: time limit for MP model
    // sp_time_lit: time limit for SP model


    // open all files for writting
    std::ofstream out_txt_f;
    out_txt_f.open(txtfile, std::ios::out);

    std::ofstream out_csv_f1;
    out_csv_f1.open(csv1, std::ios::app);

    std::ofstream out_csv_f2;
    out_csv_f2.open(csv2, std::ios::app);

    GRBEnv env = GRBEnv(true);
    env.start();

    LBBD_BS_Info* lbbdInfo = new LBBD_BS_Info();


    // ********************************************************************************************************************
    // get the lower bound of the SP (delta) and the best sequence of configuration (Best_Con_Seq)
    get_delta(env, ins, (*lbbdInfo));

    std::cout << "The LB of the subproblem (delta) is:  " << lbbdInfo->LB_delta_obj << std::endl;
    std::cout << "Time used to get delta:  " << lbbdInfo->LB_delta_time << " seconds" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    out_txt_f << "The LB of the subproblem (delta) is: " << lbbdInfo->LB_delta_obj << std::endl;
    out_txt_f << "Time used to get delta: " << lbbdInfo->LB_delta_time << std::endl;

    out_txt_f << "-----------------------------------------------------" << std::endl;
    out_csv_f1 << ins_name << ',' << std::to_string(lbbdInfo->LB_delta_obj) << ',' << std::to_string(lbbdInfo->LB_delta_time) << ',';
    out_csv_f2 << ins_name << ',';


    // ********************************************************************************************************************
    double remain_time = T_time_limit;

    // build MP model
    auto mp_t1 = std::chrono::high_resolution_clock::now();
    GRBModel MP_model = GRBModel(env);
    //MP_model.set(GRB_IntParam_OutputFlag, 0);
    MP_model.set(GRB_IntParam_Method, 2);
    MP_model.set(GRB_DoubleParam_TimeLimit, remain_time);   // time limit for solving MP

    // ------------------------ add variables
    std::map<std::tuple<int, int, int>, GRBVar> X;
    for (int i = 0; i < ins.d_T; i++)
        for (int b = 0; b < ins.B.size(); b++)
            for (int m = 0; m < ins.M.size(); m++)
            {
                std::string var_name = "x[" + std::to_string(ins.V[i]) + "]["
                    + std::to_string(ins.B[b]) + "][" + std::to_string(ins.M[m]) + "]";
                GRBVar x = MP_model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                X[std::make_tuple(ins.V[i], ins.B[b], ins.M[m])] = x;
            }

    GRBVar theta = MP_model.addVar(lbbdInfo->LB_delta_obj, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "theta");


    std::map<std::tuple<int, int>, GRBVar> F;
    for (int i = 0; i < ins.d_T - 1; i++)
    {
        std::string var_name = "f[" + std::to_string(ins.V[i]) + "][" + std::to_string(ins.V[i] + 1) + "]";
        GRBVar f = MP_model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, var_name);
        F[std::make_tuple(ins.V[i], ins.V[i] + 1)] = f;
    }


    std::map<std::tuple<int, int>, GRBVar> G;
    for (int i = 0; i < ins.d_T - 1; i++)
    {
        std::string var_name = "g[" + std::to_string(ins.V[i]) + "][" + std::to_string(ins.V[i] + 1) + "]";
        GRBVar g = MP_model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, var_name);
        G[std::make_tuple(ins.V[i], ins.V[i] + 1)] = g;
    }



    // ------------------ add constraints
    for (int i = 1; i <= ins.d_T; i++)
    {
        GRBLinExpr con1 = 0;
        for (int b = 1; b <= ins.B.size(); b++)
            for (int m = 1; m <= ins.M.size(); m++)
                con1 += X[std::make_tuple(i, b, m)];
        MP_model.addConstr(con1, GRB_EQUAL, 1.0);
    }

    for (int i = 0; i < ins.cars.size(); i++)
    {
        GRBLinExpr con2 = 0;
        for (int v = 1; v <= ins.d_T; v++)
            con2 += X[std::make_tuple(v, ins.cars[i].body, ins.cars[i].color)];
        MP_model.addConstr(con2, GRB_EQUAL, ins.cars[i].demand);
    }

    for (int p = 1; p <= ins.d_T - ins.c0; p++)
        for (int m = 1; m <= ins.M.size(); m++)
        {
            GRBLinExpr con3 = 0;
            for (int i = p; i <= p + ins.c0; i++)
                for (int b = 1; b <= ins.B.size(); b++)
                    con3 += X[std::make_tuple(i, b, m)];
            MP_model.addConstr(con3, GRB_LESS_EQUAL, ins.c0);
        }

    for (int i = 1; i < ins.d_T; i++)
        for (int b = 1; b <= ins.B.size(); b++)
        {
            GRBLinExpr con4 = F[std::make_tuple(i, i + 1)];
            for (int m = 1; m <= ins.M.size(); m++)
            {
                con4 -= X[std::make_tuple(i, b, m)];
                con4 += X[std::make_tuple(i + 1, b, m)];
            }
            MP_model.addConstr(con4, GRB_GREATER_EQUAL, 0.0);

        }

    for (int i = 1; i < ins.d_T; i++)
        for (int m = 1; m <= ins.M.size(); m++)
        {
            GRBLinExpr con5 = G[std::make_tuple(i, i + 1)];
            for (int b = 1; b <= ins.B.size(); b++)
            {
                con5 -= X[std::make_tuple(i, b, m)];
                con5 += X[std::make_tuple(i + 1, b, m)];
            }
            MP_model.addConstr(con5, GRB_GREATER_EQUAL, 0.0);

        }
    MP_model.addConstr(theta, GRB_GREATER_EQUAL, lbbdInfo->LB_delta_obj);

    // ------------------------------ set objective -----------------
    GRBLinExpr mp_obj = theta;
    for (int i = 1; i < ins.d_T; i++)
    {
        mp_obj += ins.alpha * F[std::make_tuple(i, i + 1)];
        mp_obj += ins.beta * G[std::make_tuple(i, i + 1)];
    }

    MP_model.setObjective(mp_obj, GRB_MINIMIZE);

    auto mp_t2 = std::chrono::high_resolution_clock::now();
    auto time_build = std::chrono::duration_cast<std::chrono::milliseconds>(mp_t2 - mp_t1);


    double loop_time_mp = time_build.count();
    double loop_time_sp = 0.0;
    //double obj_mp = -1;
    //double obj_f1f2 = -1;
    //double obj_theta = -1;
    //double obj_sp = -1;
    std::map< std::tuple<int, int>, std::vector<int> > V_bm;

    lbbdInfo->Total_time += lbbdInfo->LB_delta_time;
    bool mpsp = true;


    // begin LBBD loop
    while (lbbdInfo->UB - lbbdInfo->LB > 1e-3)
    {
        if (lbbdInfo->Total_time > T_time_limit * 1000)
            break;
        if (remain_time < 2.0)
            break;

        if (!lbbdInfo->lB.empty())
        {
            lbbdInfo->lB.clear();
            lbbdInfo->lC.clear();
            lbbdInfo->lO.clear();
            lbbdInfo->lA.clear();
            lbbdInfo->lE.clear();
            lbbdInfo->MP_status = -1;
            lbbdInfo->SP_status = -1;
            lbbdInfo->OP_status = -1;
            lbbdInfo->mp_gap = -1;
            lbbdInfo->sp_gap = -1;
            V_bm.clear();
            lbbdInfo->MP_status = -1;
            lbbdInfo->SP_status = -1;
        }

        bool add_opt = 0;
        bool add_nogood = 0;
        loop_time_sp = 0.0;

        lbbdInfo->LBBD_iter += 1;
        std::cout << "Iteration: " << lbbdInfo->LBBD_iter << std::endl;

        // solve MP first
        auto mp_tt1 = std::chrono::high_resolution_clock::now();
        MP_model.optimize();
        auto mp_tt2 = std::chrono::high_resolution_clock::now();
        auto mp_dur = std::chrono::duration_cast<std::chrono::milliseconds>(mp_tt2 - mp_tt1);

        remain_time -= MP_model.get(GRB_DoubleAttr_Runtime);

        loop_time_mp = MP_model.get(GRB_DoubleAttr_Runtime);

        //if (lbbdInfo->LBBD_iter == 1)

        //    loop_time_mp += mp_dur.count();
        //else
        //    loop_time_mp = mp_dur.count();


        lbbdInfo->total_time_MP += loop_time_mp;
        lbbdInfo->Total_time += loop_time_mp;

        if (MP_model.get(GRB_IntAttr_Status) != 2)
        {
            if (MP_model.get(GRB_DoubleAttr_MIPGap) < 10)
                lbbdInfo->mp_gap = MP_model.get(GRB_DoubleAttr_MIPGap);

            mpsp = false;

            if (lbbdInfo->LBBD_iter == 1)
                out_csv_f2 << std::to_string(lbbdInfo->LBBD_iter) << ',';
            else
                out_csv_f2 << " " << ',' << std::to_string(lbbdInfo->LBBD_iter) << ',';

            out_csv_f2 << std::to_string(lbbdInfo->UB) << ','
                << std::to_string(lbbdInfo->LB) << ','
                << std::to_string(loop_time_mp * 0.001) << ','
                << std::to_string(lbbdInfo->MP_status) << ',' << "-1" << ',' << "-1" << ','
                << "-1" << ','
                << std::to_string(lbbdInfo->mp_gap) << ','
                << std::to_string(0) << ','
                << std::to_string(-1) << ','
                << std::to_string(-1) << std::endl;
            break;
        }

        lbbdInfo->MP_status = 1;
        lbbdInfo->mp_gap = 0.0;

        std::cout << "----------------Solve SP.--------------------" << std::endl;


        lbbdInfo->obj_MP = MP_model.get(GRB_DoubleAttr_ObjVal);
        lbbdInfo->obj_mp_theta = theta.get(GRB_DoubleAttr_X);
        lbbdInfo->obj_mp_f1f2 = lbbdInfo->obj_MP - lbbdInfo->obj_mp_theta;

        lbbdInfo->LB = lbbdInfo->obj_MP;

        // get the solution from MP
        std::vector < std::tuple<int, int, int> > X_ibm;
        for (int i = 1; i <= ins.d_T; i++)
            for (int b = 1; b <= ins.B.size(); b++)
                for (int m = 1; m <= ins.M.size(); m++)
                    if (X[std::make_tuple(i, b, m)].get(GRB_DoubleAttr_X) > 0.05)
                    {
                        lbbdInfo->lB.push_back(b);
                        lbbdInfo->lC.push_back(m);
                        V_bm[std::make_tuple(b, m)].push_back(i);
                        X_ibm.push_back(std::make_tuple(i, b, m));      
                        std::cout << "(" << i << "," << b << "," << m << ")-";
                    }
        std::cout << std::endl;


        // solve SP
        InformationSolveSP* spInfo = new InformationSolveSP();
        Solve_SP(*spInfo, ins, V_bm, env, remain_time);

        remain_time -= spInfo->SP_total_time;

        lbbdInfo->lastSP = spInfo;
        loop_time_sp += spInfo->SP_total_time;
        lbbdInfo->total_time_SP += loop_time_sp;
        lbbdInfo->Total_time += loop_time_sp;

        lbbdInfo->obj_SP = spInfo->SP_OptObjVal;
        lbbdInfo->SP_status = spInfo->SP_status;
        lbbdInfo->sp_gap = spInfo->SP_gap;

        if (lbbdInfo->SP_status == -1)
        {
            mpsp = false;
            if (lbbdInfo->LBBD_iter == 1)
                out_csv_f2 << std::to_string(lbbdInfo->LBBD_iter) << ',';
            else
                out_csv_f2 << " " << ',' << std::to_string(lbbdInfo->LBBD_iter) << ',';

            out_csv_f2 << std::to_string(lbbdInfo->UB) << ','
                << std::to_string(lbbdInfo->LB) << ','
                << std::to_string(loop_time_mp) << ','
                << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->obj_MP) << ',' << std::to_string(lbbdInfo->obj_mp_theta) << ','
                << std::to_string(lbbdInfo->obj_MP - lbbdInfo->obj_mp_theta) << ','
                << std::to_string(lbbdInfo->mp_gap) << ','
                << std::to_string(loop_time_sp) << ','
                << std::to_string(lbbdInfo->SP_status) << ','
                << std::to_string(lbbdInfo->obj_SP) << ',' << std::to_string(lbbdInfo->sp_gap) << std::endl;
            break;
        }
        else
        {
            for (int dd = 0; dd < ins.d_T; dd++)
            {
                lbbdInfo->lA.push_back(spInfo->A[dd]);
                lbbdInfo->lO.push_back(spInfo->O[dd]);
                lbbdInfo->lE.push_back(spInfo->E[dd]);
            }
        }

        // judge if add lbbd cuts
        if (lbbdInfo->UB - (lbbdInfo->obj_MP - lbbdInfo->obj_mp_theta + lbbdInfo->obj_SP) > 1e-2)
        {
            add_opt = 1;
            lbbdInfo->lbbd_cut_opt += 1;
            lbbdInfo->UB = lbbdInfo->obj_MP - lbbdInfo->obj_mp_theta + lbbdInfo->obj_SP;
            GRBLinExpr opt_cut = 0;
            for (int opt = 0; opt < X_ibm.size(); opt++)
                opt_cut += X[X_ibm[opt]];

            MP_model.addConstr(theta, GRB_GREATER_EQUAL, lbbdInfo->obj_SP * (opt_cut - X_ibm.size() + 1));
        }
        else
        {
            add_nogood = 1;
            lbbdInfo->lbbd_cut_nogood += 1;
            GRBLinExpr opt_cut = 0;
            for (int opt = 0; opt < X_ibm.size(); opt++)
                opt_cut += X[X_ibm[opt]];
            MP_model.addConstr(opt_cut, GRB_LESS_EQUAL, X_ibm.size() - 1);

        }


        // write csv2
        if (lbbdInfo->LBBD_iter == 1)
            out_csv_f2 << std::to_string(lbbdInfo->LBBD_iter) << ',';
        else
            out_csv_f2 << " " << ',' << std::to_string(lbbdInfo->LBBD_iter) << ',';

        out_csv_f2 << std::to_string(lbbdInfo->UB) << ','
            << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(loop_time_mp) << ','
            << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->obj_MP) << ',' << std::to_string(lbbdInfo->obj_mp_theta) << ','
            << std::to_string(lbbdInfo->obj_MP - lbbdInfo->obj_mp_theta) << ','
            << std::to_string(lbbdInfo->mp_gap) << ','
            << std::to_string(loop_time_sp) << ','
            << std::to_string(lbbdInfo->SP_status) << ','
            << std::to_string(lbbdInfo->obj_SP) << ',' << std::to_string(lbbdInfo->sp_gap) << std::endl;



        loop_time_mp = 0.0;
        loop_time_sp = 0.0;
        delete spInfo;
    }


    if (mpsp == false)
    {
        std::cout << "Do not get the optimal solution because MP or SP is not solved optimally." << std::endl;
        if (lbbdInfo->LBBD_iter == 1)
            lbbdInfo->lbbd_gap = -1;
        else
            lbbdInfo->lbbd_gap = (lbbdInfo->UB - lbbdInfo->LB) / lbbdInfo->LB;

    }
    else
    {
        if (lbbdInfo->UB - lbbdInfo->LB > 1e-2)
        {
            std::cout << "Do not get the optimal solution because UB != LB." << std::endl;
            lbbdInfo->lbbd_gap = (lbbdInfo->UB - lbbdInfo->LB) / lbbdInfo->LB;

            out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time) << std::endl;
            out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
            out_txt_f << "MP Status: " << std::to_string(lbbdInfo->MP_status) << std::endl;
            out_txt_f << "SP Status: " << std::to_string(lbbdInfo->SP_status) << std::endl;
            out_txt_f << "OP Status: " << std::to_string(lbbdInfo->OP_status) << std::endl;
            out_txt_f << "Current upper bound: " << std::to_string(lbbdInfo->UB) << std::endl;
            out_txt_f << "Current lower bound: " << std::to_string(lbbdInfo->LB) << std::endl;
            out_txt_f << "Total Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP) << std::endl;
            out_txt_f << "Total Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP) << std::endl;

        }
        else
        {
            std::cout << "Get the optimal solution." << std::endl;
            lbbdInfo->lbbd_gap = 0;
            lbbdInfo->OP_status = 1;

            // write txt
            out_txt_f << "Get the optimal solution!" << std::endl;
            out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time) << std::endl;
            out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
            out_txt_f << "Total Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP) << std::endl;
            out_txt_f << "Total Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP) << std::endl;
            out_txt_f << "Avg.Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP / lbbdInfo->LBBD_iter) << std::endl;
            out_txt_f << "Avg.Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP / lbbdInfo->LBBD_iter) << std::endl;

            out_txt_f << "LBBD opt_cut: " << std::to_string(lbbdInfo->lbbd_cut_opt) << std::endl;
            out_txt_f << "LBBD nogood_cut: " << std::to_string(lbbdInfo->lbbd_cut_nogood) << std::endl;

            out_txt_f << "--------------------------------------------------------" << std::endl;
            out_txt_f << "The body type of each car: " << std::endl;
            for (int bb = 0; bb < lbbdInfo->lB.size(); bb++) out_txt_f << std::to_string(lbbdInfo->lB[bb]) << " - ";
            out_txt_f << std::endl;

            out_txt_f << "The color of each car: " << std::endl;
            for (int cc = 0; cc < lbbdInfo->lC.size(); cc++) out_txt_f << std::to_string(lbbdInfo->lC[cc]) << " - ";
            out_txt_f << std::endl;

            out_txt_f << "The configuration of each car: " << std::endl;
            for (int oo = 0; oo < lbbdInfo->lO.size(); oo++) out_txt_f << std::to_string(lbbdInfo->lO[oo]) << " - ";
            out_txt_f << std::endl;

            out_txt_f << "The lane to which each car is allcoated: " << std::endl;
            for (int ll = 0; ll < lbbdInfo->lA.size(); ll++) out_txt_f << std::to_string(lbbdInfo->lA[ll]) << " - ";
            out_txt_f << std::endl;

            out_txt_f << "The upstream car_id in each downstream position: " << std::endl;
            for (int ee = 0; ee < lbbdInfo->lE.size(); ee++) out_txt_f << std::to_string(lbbdInfo->lE[ee]) << " - ";
            out_txt_f << std::endl;
        }

    }


    // write csv1
    out_csv_f1 << std::to_string(lbbdInfo->OP_status) << ',' << std::to_string(lbbdInfo->MP_status) << ','
        << std::to_string(lbbdInfo->SP_status) << ',' << std::to_string(lbbdInfo->LBBD_iter) << ',';

    if (mpsp == false && lbbdInfo->LBBD_iter == 1)
        out_csv_f1 << "-1" << ',';
    else
    {
        out_csv_f1 << std::to_string(lbbdInfo->lbbd_gap) << ',';
    }

    // UB, LB, objmp, objsp, num.opt.cut, num.nogood.cut
    out_csv_f1 << std::to_string(lbbdInfo->Total_time) << ',' << std::to_string(lbbdInfo->total_time_MP) << ','
        << std::to_string(lbbdInfo->total_time_SP) << ',' << std::to_string(lbbdInfo->total_time_MP / lbbdInfo->LBBD_iter) << ','
        << std::to_string(lbbdInfo->total_time_SP / lbbdInfo->LBBD_iter) << ','
        << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ',' << std::to_string(lbbdInfo->obj_MP) << ','
        << std::to_string(lbbdInfo->obj_mp_f1f2) << ',' << std::to_string(lbbdInfo->obj_mp_theta) << ',' << std::to_string(lbbdInfo->obj_SP) << ','
        << std::to_string(lbbdInfo->lbbd_cut_opt) << ',' << std::to_string(lbbdInfo->lbbd_cut_nogood) << std::endl;

    out_csv_f1.close();
    out_csv_f2.close();
    out_txt_f.close();

}