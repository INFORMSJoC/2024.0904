#pragma once
#include "HeuInitial.h"
#include "BuildG2.h"

using namespace std;

// -------------- get the lower bound of the SP (delta)
void get_delta(GRBEnv& env, Instance ins, LBBD_Info& lbbd)
{
    auto t1 = std::chrono::steady_clock::now();

    GRBModel model = GRBModel(env);
    model.set(GRB_IntParam_OutputFlag, 0);

    // add decision variables
    std::vector< std::vector<GRBVar> > N(ins.d_T);   //N[p][j]: the number of cars with configuration j from position 1 to p
    for (int p = 0; p < ins.d_T; p++)
    {
        std::vector<GRBVar> row(ins.num_j);
        for (int j = 0; j < ins.num_j; j++)
        {
            std::string var_name = "n[" + std::to_string(p + 1) + "][" + std::to_string(j + 1) + "]";
            GRBVar n = model.addVar(0.0, (double)ins.d_T, 0.0, GRB_INTEGER, var_name);
            row[j] = n;
        }
        N[p] = row;
    }

    std::vector< std::vector<GRBVar> > U(ins.d_T);     //U[p][w]: for linearization
    for (int p = 0; p < ins.d_T; p++)
    {
        std::vector<GRBVar> row(ins.num_w);
        for (int w = 0; w < ins.num_w; w++)
        {
            std::string var_name = "u[" + std::to_string(p + 1) + "][" + std::to_string(w + 1) + "]";
            GRBVar u = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, var_name);
            row[w] = u;
        }
        U[p] = row;
    }


    // add constraints
    for (int p = 0; p < ins.d_T; p++)              // p_v = ins.V[p] = 1, 2, 3, 4, ...
    {
        for (int w = 0; w < ins.num_w; w++)          // w_v = ins.W[w] = 1, 2, 3, ...
        {
            GRBLinExpr cons1 = 0;
            for (auto it = ins.r_j_w.begin(); it != ins.r_j_w.end(); it++)
                cons1 += it->second[w + 1] * N[p][it->first - 1];
            model.addConstr(U[p][w] - cons1 + (p + 1) * ins.sigma_w[w + 1], GRB_GREATER_EQUAL, 0.0);
            model.addConstr(U[p][w] + cons1 - (p + 1) * ins.sigma_w[w + 1], GRB_GREATER_EQUAL, 0.0);
        }
    }

    for (int p = 0; p < ins.d_T; p++)
    {
        GRBLinExpr cons2 = 0;
        for (int j = 0; j < ins.num_j; j++)
            cons2 += N[p][j];
        model.addConstr(cons2 - (p + 1), GRB_EQUAL, 0.0);
    }

    for (int p = 0; p < ins.d_T - 1; p++)
    {
        for (int j = 0; j < ins.num_j; j++)
        {
            model.addConstr(N[p][j] - N[p + 1][j], GRB_LESS_EQUAL, 0.0);
            model.addConstr(N[p + 1][j] - N[p][j] - 1, GRB_LESS_EQUAL, 0.0);
        }
    }

    // set objective function
    GRBLinExpr obj = 0;
    for (int p = 0; p < ins.d_T; p++)
        for (int w = 0; w < ins.num_w; w++)
            obj += U[p][w];
    model.setObjective(ins.gamma * obj, GRB_MINIMIZE);

    // optimize model
    model.optimize();
    auto t2 = std::chrono::steady_clock::now();
    auto du = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    // get the solution
    for (int j = 0; j < ins.num_j; j++)
        if (N[0][j].get(GRB_DoubleAttr_X) > 0.5)
        {
            lbbd.LB_Config_Seq.push_back(j + 1);
            break;
        }
    for (int p = 1; p < ins.d_T; p++)
        for (int j = 0; j < ins.num_j; j++)
            if (N[p][j].get(GRB_DoubleAttr_X) - N[p - 1][j].get(GRB_DoubleAttr_X) > 0.5)
            {
                lbbd.LB_Config_Seq.push_back(j + 1);
                break;
            }


    std::cout << "Configuration sequence in LB: [" << lbbd.LB_Config_Seq[0];
    for (int i = 1; i < lbbd.LB_Config_Seq.size(); i++)
        std::cout << "--> " << lbbd.LB_Config_Seq[i];
    std::cout << "]" << std::endl;

    std::cout << "LB of SP (delta): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;

    lbbd.LB_delta_obj = model.get(GRB_DoubleAttr_ObjVal);
    lbbd.LB_delta_time = model.get(GRB_DoubleAttr_Runtime) * 1000;

}



// -------------- define LBBD function on G2
void LBBDg2(Instance ins, std::string const& ins_name, std::string const& txtfile,
    std::string const& csv1, std::string const& csv2, 
    bool& addCut1, bool& addCut2, bool& addCut3, bool& LiftCut2, bool& LiftCut3, int& sp2_method, bool& use_Heu, double time_limit)
{
    // ins_name: 10_1 for example
    // txtfile: output the final result
    // csv1: record the final summarized results
    // csv2: record the information for each iteration / shortest-path
    // use_Heu: whether use heuristic to solve SP after getting a solution to MP
    // time_limit: for the whole LBBD 1hours for example


    // open all files for writting
    std::ofstream out_txt_f;
    out_txt_f.open(txtfile, std::ios::out);

    std::ofstream out_csv_f1;
    out_csv_f1.open(csv1, std::ios::app);

    std::ofstream out_csv_f2;
    out_csv_f2.open(csv2, std::ios::app);


    GRBEnv env = GRBEnv(true);
    env.start();

    LBBD_Info* lbbdInfo = new LBBD_Info();
    optSol* OptimalSol = new optSol();

    // ********************************************************************************************************************
    // get the lower bound of the SP (delta) and the best sequence of configuration (Best_Con_Seq)
    get_delta(env, ins, (*lbbdInfo));

    std::cout << "The LB of the subproblem (delta) is:  " << lbbdInfo->LB_delta_obj << std::endl;
    std::cout << "Time used to get delta:  " << lbbdInfo->LB_delta_time * 0.001 << " seconds" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    out_txt_f << "The LB of the subproblem (delta) is: " << lbbdInfo->LB_delta_obj << std::endl;
    out_txt_f << "Time used to get delta: " << lbbdInfo->LB_delta_time * 0.001 << std::endl;
    out_txt_f << "-----------------------------------------------------" << std::endl;
    out_csv_f1 << ins_name << ',' << std::to_string(lbbdInfo->LB_delta_obj) << ',' << std::to_string(lbbdInfo->LB_delta_time * 0.001) << ',';
    out_csv_f2 << ins_name << ',';



    // ********************************************************************************************************************
    // build graph2
    Graph2 g2;
    BuildGraph2(g2, ins);

    lbbdInfo->graph2 = &g2;

    out_txt_f << "Num.Nodes.G2: " << g2.nodes.size() << std::endl;
    out_txt_f << "Num.Arcs.G2: " << g2.arcs.size() << std::endl;
    out_txt_f << "Time.Build.G2: " << g2.tt1 * 0.001 << " seconds" << std::endl;
    out_txt_f << "-----------------------------------------------------" << std::endl;
    out_csv_f1 << std::to_string(g2.nodes.size()) << ',' << std::to_string(g2.arcs.size()) << ',';
    out_csv_f1 << std::to_string(g2.tt1 * 0.001) << ',';


    // ********************************************************************************************************************
    // *************** Beigin LBBD on G2 **********************************************************************************
    // build the master problem MP
    auto mp_t1 = std::chrono::steady_clock::now();
    GRBModel MP_model = GRBModel(env);

    // ----------------------add variables
    std::map< std::tuple<int, int>, GRBVar > X;

    for (int i = 0; i < g2.arcs.size(); i++)
    {
        string var_name = "x[" + std::to_string(g2.arcs[i]->from) + "][" + std::to_string(g2.arcs[i]->to) + "]";
        GRBVar x = MP_model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        std::tuple<int, int> lab(g2.arcs[i]->from, g2.arcs[i]->to);
        X[lab] = x;
    }


    GRBVar theta = MP_model.addVar(lbbdInfo->LB_delta_obj, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "theta");

    // ---------------------add MP constraints-type2----------------
    GRBLinExpr mp_con1 = 0;
    for (int i = 0; i < g2.arcs_0.size(); i++)
        mp_con1 += X[std::make_tuple(0, g2.arcs_0[i]->to)];
    MP_model.addConstr(mp_con1, GRB_EQUAL, 1.0);


    GRBLinExpr mp_con2 = 0;
    for (int i = 0; i < g2.nodes_terminal.size(); i++)
    {
        int a_to = g2.nodes_terminal[i]->id;
        for (int a = 0; a < g2.in_A[a_to].size(); a++)
            mp_con2 += X[std::make_tuple(g2.in_A[a_to][a]->from, a_to)];
    }
    MP_model.addConstr(mp_con2, GRB_EQUAL, 1.0);


    for (int n11 = 0; n11 < g2.nodes_inter.size(); n11++)
    {
        int n1 = g2.nodes_inter[n11]->id;
        GRBLinExpr item1 = 0;
        for (int a_to = 0; a_to < g2.out_A[n1].size(); a_to++)
        {
            int n2 = g2.out_A[n1][a_to]->to;
            item1 += X[std::make_tuple(n1, n2)];
        }
        GRBLinExpr item2 = 0;
        for (int a_from = 0; a_from < g2.in_A[n1].size(); a_from++)
        {
            int n2 = g2.in_A[n1][a_from]->from;
            item2 += X[std::make_tuple(n2, n1)];
        }
        MP_model.addConstr(item1, GRB_EQUAL, item2);

    }


    for (int i = 0; i < ins.cars.size(); i++)
    {
        int b = ins.cars[i].body;
        int m = ins.cars[i].color;
        GRBLinExpr mp_con3 = 0;
        for (int a = 0; a < g2.arcs.size(); a++)
        {
            if (g2.arcs[a]->body == b && g2.arcs[a]->color == m)
                mp_con3 += X[make_tuple(g2.arcs[a]->from, g2.arcs[a]->to)];
        }
        MP_model.addConstr(mp_con3, GRB_GREATER_EQUAL, ins.cars[i].demand);
    }


    MP_model.addConstr(theta, GRB_GREATER_EQUAL, lbbdInfo->LB_delta_obj);



    // ------------------------ set objective -----------------------
    GRBLinExpr mp_obj = 0;
    for (int i = 0; i < g2.arcs.size(); i++)
        mp_obj += g2.arcs[i]->cost * X[std::make_tuple(g2.arcs[i]->from, g2.arcs[i]->to)];
    MP_model.setObjective(mp_obj + theta, GRB_MINIMIZE);

    auto mp_t2 = std::chrono::steady_clock::now();
    auto time_build = std::chrono::duration_cast<std::chrono::milliseconds>(mp_t2 - mp_t1);

    // lbbdInfo->Total_time += time_build.count();
    // lbbdInfo->total_time_MP += time_build.count();


    double loop_time_mp = time_build.count();
    double loop_time_sp = 0.0;
    lbbdInfo->Total_time += lbbdInfo->LB_delta_time + g2.tt1;

    double remain_time = time_limit - lbbdInfo->Total_time * 0.001;

    // MP_model.set(GRB_IntParam_OutputFlag, 0);
    MP_model.set(GRB_IntParam_Method, 2);
    MP_model.set(GRB_DoubleParam_TimeLimit, remain_time);   // time limit for solving MP


    bool mpsp_status = true;
    bool is_heu_opt = false;


    std::vector<int> Body;
    std::vector<int> Color;
    std::map<std::tuple<int, int>, std::vector<int>> V_bm;
    std::vector<int> EE;
    std::vector<int> OO;
    std::vector<int> AA;

    std::map<int, std::vector<std::tuple<int, int>>> Config_BodyColor;
    StatConBodyColor(ins, Config_BodyColor);
    InformationSolveSP* spInfo = new InformationSolveSP();
    lbbdInfo->lastSP = spInfo;
    std::map<double, std::vector<SP1_Node_Sol*>> SP_feasibleSolList;

    // begin LBBD loop
    while (lbbdInfo->UB - lbbdInfo->LB > 1e-3)
    {
        if (remain_time < 1.0)
            break;

        if (!Body.empty())
        {
            Body.clear();
            Color.clear();
            V_bm.clear();
            lbbdInfo->lE.clear();
            lbbdInfo->lA.clear();
            lbbdInfo->lO.clear();
            lbbdInfo->mp_gap = -1;
            lbbdInfo->sp_gap = -1;
            lbbdInfo->MP_status = -1;
            lbbdInfo->SP_status = -1;

            for (auto it = SP_feasibleSolList.begin(); it != SP_feasibleSolList.end(); it++)
                for (int inode = 0; inode < it->second.size(); inode++)
                    delete it->second[inode];
            SP_feasibleSolList.clear();

            spInfo->spInfoUpdate();
        }

        int add_opt = 0;
        int add_nogood = 0;
        is_heu_opt = false;
        loop_time_sp = 0.0;
        mpsp_status = true;

        lbbdInfo->LBBD_iter += 1;
        std::cout << "Iteration: " << lbbdInfo->LBBD_iter << std::endl;
        // solve MP first 
        MP_model.optimize();

        if (lbbdInfo->LBBD_iter == 1)
            loop_time_mp += MP_model.get(GRB_DoubleAttr_Runtime) * 1000;
        else
            loop_time_mp = MP_model.get(GRB_DoubleAttr_Runtime) * 1000;

        remain_time -= loop_time_mp * 0.001;
        if (MP_model.get(GRB_IntAttr_Status) != 2)
        {
            if (MP_model.get(GRB_DoubleAttr_MIPGap) > 10)
                lbbdInfo->mp_gap = -1;
            else
                lbbdInfo->mp_gap = MP_model.get(GRB_DoubleAttr_MIPGap);

            lbbdInfo->total_time_MP += loop_time_mp;
            lbbdInfo->Total_time += loop_time_mp;
            mpsp_status = false;
            break;
        }

        lbbdInfo->MP_status = 1;
        lbbdInfo->mp_gap = 0.0;

        std::cout << "----------- Finish solving MP." << std::endl;

        lbbdInfo->obj_MP = MP_model.get(GRB_DoubleAttr_ObjVal);
        lbbdInfo->obj_mp_theta = theta.get(GRB_DoubleAttr_X);
        lbbdInfo->obj_mp_f1f2 = lbbdInfo->obj_MP - lbbdInfo->obj_mp_theta;


        lbbdInfo->LB = lbbdInfo->obj_MP;
        lbbdInfo->total_time_MP += loop_time_mp;
        lbbdInfo->Total_time += loop_time_mp;


        // get the solution from MP
        int car_n = 0;
        for (int a = 0; a < g2.arcs.size(); a++)
        {
            Arc_* aa = g2.arcs[a];
            if (X[std::make_tuple(aa->from, aa->to)].get(GRB_DoubleAttr_X) > 0.05)
            {
                Body.push_back(aa->body);
                Color.push_back(aa->color);
                car_n += 1;
                V_bm[make_tuple(aa->body, aa->color)].push_back(car_n);
                std::cout << "car: " << car_n << " body: " << aa->body << " color: " << aa->color << std::endl;
            }
        }


        // solve SP
        if (use_Heu == true && fabs(lbbdInfo->obj_mp_theta - lbbdInfo->LB_delta_obj) < 1e-3)
        {
            // solve SP by heursitic first, if it is infeasible, then solve by B&C
            auto sp_h_t1 = std::chrono::steady_clock::now();
            bool is_heu_feasible = HeuDone(env, ins, Config_BodyColor, Body, Color, *lbbdInfo, sp2_method);
            auto sp_h_t2 = std::chrono::steady_clock::now();
            auto sp_h_dur = std::chrono::duration_cast<std::chrono::milliseconds>(sp_h_t2 - sp_h_t1);

            loop_time_sp = sp_h_dur.count();

            if (is_heu_feasible == true)
            {
                //lbbdInfo->total_SP_solved += spInfo->num_SP_solved;
                lbbdInfo->SP_status = 1;
                lbbdInfo->sp_gap = 0.0;
                lbbdInfo->obj_SP = lbbdInfo->LB_delta_obj;

                std::cout << "---------- Finsh solving SP by Heu." << std::endl;
                is_heu_opt = true;

            }
            else
            {
                // solve SP by B&B
                Solve_SP(*spInfo, SP_feasibleSolList, ins, V_bm, env, addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, remain_time);
                std::cout << "---------- Finsh solving SP by B&C." << std::endl;

                lbbdInfo->total_SP_solved += spInfo->num_SP_solved;
                loop_time_sp += spInfo->SP_total_time;

                lbbdInfo->total_lazy1 += spInfo->num_lazy1;
                lbbdInfo->total_lazy2 += spInfo->num_lazy2;
                lbbdInfo->total_lazy3 += spInfo->num_lazy3;

                lbbdInfo->total_num_SP2_NogoodCut += spInfo->num_nogood_cut;
                lbbdInfo->total_num_SP2_solved += spInfo->num_sp2_solved;
                lbbdInfo->total_time_SP2 += spInfo->time_solve_sp2;

                lbbdInfo->SP_status = spInfo->SP_status;
                lbbdInfo->sp_gap = spInfo->SP_gap;


                if (lbbdInfo->SP_status == -1)
                {
                    // do not find the optimal solution to SP
                    mpsp_status = false;
                    lbbdInfo->obj_SP = -1;
                    lbbdInfo->Total_time += loop_time_sp;
                    lbbdInfo->total_time_SP += loop_time_sp;
                    remain_time -= loop_time_sp * 0.001;
                    break;
                }
                else
                {
                    lbbdInfo->obj_SP = spInfo->SP_OptObjVal;
                    for (int dd = 0; dd < ins.d_T; dd++)
                    {
                        lbbdInfo->lE.push_back(spInfo->E[dd]);
                        lbbdInfo->lO.push_back(spInfo->O[dd]);
                        lbbdInfo->lA.push_back(spInfo->A[dd]);

                    }
                }
            }
        }
        else
        {
            // solve SP by B&B
            Solve_SP(*spInfo, SP_feasibleSolList, ins, V_bm, env, addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, remain_time);
            std::cout << "---------- Finsh solving SP by B&C." << std::endl;

            loop_time_sp += spInfo->SP_total_time;
            lbbdInfo->total_SP_solved += spInfo->num_SP_solved;

            lbbdInfo->total_lazy1 += spInfo->num_lazy1;
            lbbdInfo->total_lazy2 += spInfo->num_lazy2;
            lbbdInfo->total_lazy3 += spInfo->num_lazy3;

            lbbdInfo->total_num_SP2_NogoodCut += spInfo->num_nogood_cut;
            lbbdInfo->total_num_SP2_solved += spInfo->num_sp2_solved;
            lbbdInfo->total_time_SP2 += spInfo->time_solve_sp2;

            lbbdInfo->SP_status = spInfo->SP_status;
            lbbdInfo->sp_gap = spInfo->SP_gap;


            if (lbbdInfo->SP_status == -1)
            {
                // do not find the optimal solution to SP
                mpsp_status = false;
                lbbdInfo->obj_SP = -1;
                lbbdInfo->Total_time += loop_time_sp;
                lbbdInfo->total_time_SP += loop_time_sp;
                remain_time -= loop_time_sp * 0.001;
                break;
            }
            else
            {
                lbbdInfo->obj_SP = spInfo->SP_OptObjVal;
                for (int dd = 0; dd < ins.d_T; dd++)
                {
                    lbbdInfo->lE.push_back(spInfo->E[dd]);
                    lbbdInfo->lO.push_back(spInfo->O[dd]);
                    lbbdInfo->lA.push_back(spInfo->A[dd]);
                }
            }
        }


        lbbdInfo->total_time_SP += loop_time_sp;
        lbbdInfo->Total_time += loop_time_sp;
        remain_time -= loop_time_sp * 0.001;

        lbbdInfo->obj_SP = spInfo->SP_OptObjVal;

        // judge if add lbbd cuts

        std::vector<std::tuple<int, int>> X_pos;
        for (auto it = X.begin(); it != X.end(); it++)
            if (it->second.get(GRB_DoubleAttr_X) > 0.5)
            {
                X_pos.push_back(it->first);
            }

        if (lbbdInfo->UB - (lbbdInfo->obj_MP - lbbdInfo->obj_mp_theta + lbbdInfo->obj_SP) > 1e-2)
        {
            lbbdInfo->lbbd_cut_opt += 1;
            lbbdInfo->UB = lbbdInfo->obj_MP - lbbdInfo->obj_mp_theta + lbbdInfo->obj_SP;
            GRBLinExpr opt_cut = 0;
            for (int opt = 0; opt < X_pos.size(); opt++)
                opt_cut += X[X_pos[opt]];

            MP_model.addConstr(theta, GRB_GREATER_EQUAL, lbbdInfo->obj_SP * (opt_cut - X_pos.size() + 1));
        }
        else
        {
            lbbdInfo->lbbd_cut_nogood += 1;
            GRBLinExpr opt_cut = 0;
            for (int opt = 0; opt < X_pos.size(); opt++)
                opt_cut += X[X_pos[opt]];

            MP_model.addConstr(opt_cut, GRB_LESS_EQUAL, X_pos.size() - 1);
        }


        // write csv2
        // f2: Iter th - UB - LB - MP.time - mp.status - mp.obj - mp.gap - SP.time - sp.status - sp.obj - sp.gap - num.lazy1 - num.lazy2 
        // - num.lazy3 - num.sp.nogood
        // num.sp2.solved - time.sp2.solved - avg.time.sp2 - num.integer.found - num.feasible.solution.found - 
        // num.lazy3.added - num.reduced.lazy3.added
        if (lbbdInfo->LBBD_iter == 1)
            out_csv_f2 << std::to_string(lbbdInfo->LBBD_iter) << ',';
        else
            out_csv_f2 << " " << ',' << std::to_string(lbbdInfo->LBBD_iter) << ',';

        out_csv_f2 << std::to_string(lbbdInfo->UB) << ','
            << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(loop_time_mp * 0.001) << ','
            << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->obj_MP) << ',' << std::to_string(lbbdInfo->obj_mp_theta) << ','
            << std::to_string(lbbdInfo->obj_MP - lbbdInfo->obj_mp_theta) << ','
            << std::to_string(lbbdInfo->mp_gap) << ','
            << std::to_string(loop_time_sp * 0.001) << ','
            << std::to_string(lbbdInfo->SP_status) << ','
            << std::to_string(lbbdInfo->obj_SP) << ',';

        if (is_heu_opt == true)
        {
            // solve by heuristic
            out_csv_f2 << std::to_string(0) << ',' << std::to_string(0) << ','
                << std::to_string(0) << ',' << std::to_string(0) << ',' << std::to_string(0) << ','
                << std::to_string(lbbdInfo->lastSP->num_sp2_solved) << ',' << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001) << ','
                << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001) << ','
                << std::to_string(1) << ',' << std::to_string(1) << ','
                << std::to_string(0) << ',' << std::to_string(0) << std::endl;
        }
        else
        {
            // solve by B&C
            out_csv_f2 << std::to_string(spInfo->SP_gap) << ',' << std::to_string(spInfo->num_lazy1) << ','
                << std::to_string(spInfo->num_lazy2) << ',' << std::to_string(spInfo->num_lazy3) << ','
                << std::to_string(spInfo->num_nogood_cut) << ','
                << std::to_string(spInfo->num_sp2_solved) << ',' << std::to_string(spInfo->time_solve_sp2 * 0.001) << ',';
            if (spInfo->num_sp2_solved == 0)
                out_csv_f2 << std::to_string(0) << ',';
            else

                out_csv_f2 << std::to_string(spInfo->time_solve_sp2 * 0.001 / spInfo->num_sp2_solved) << ',';

            out_csv_f2 << std::to_string(spInfo->num_integer_sol) << ',' << std::to_string(spInfo->num_feasible_sol) << ',';

            //if (spInfo->numLazy3List.size() == 0)
            //    out_csv_f2 << std::to_string(0) << ',' << std::to_string(0);
            //else
            //{
            //    std::string ss = "[";
            //    for (int nn = 0; nn < spInfo->numLazy3List.size() - 1; nn++)
            //    {
            //        ss += std::to_string(spInfo->numLazy3List[nn]);
            //        ss += "-";
            //    }
            //    ss += std::to_string(spInfo->numLazy3List.back());
            //    ss += "]";
            //    out_csv_f2 << ss << ',';

            //    std::string ss3 = "[";
            //    for (int nn = 0; nn < spInfo->reduced_numlazy3List.size() - 1; nn++)
            //    {
            //        ss3 += std::to_string(spInfo->reduced_numlazy3List[nn]);
            //        ss3 += "-";
            //    }
            //    ss3 += std::to_string(spInfo->reduced_numlazy3List.back());
            //    ss3 += "]";
            //    out_csv_f2 << ss3;
            //}


            out_csv_f2 << std::endl;

        }


        loop_time_mp = 0.0;
        loop_time_sp = 0.0;

    }


    if (is_heu_opt == true)     // get the optimal solution by heuristic
        lbbdInfo->Heu_sp = 1;
    else
        lbbdInfo->Heu_sp = 0;


    // wirte txt file and csv1
    if (mpsp_status == false)
    {
        std::cout << "Do not find the optimal solution because either MP or SP is infeasible." << std::endl;

        out_txt_f << "Do not find the optimal solution because either MP or SP is infeasible." << std::endl;
        out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
        out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "MP Status: " << std::to_string(lbbdInfo->MP_status) << std::endl;
        out_txt_f << "SP Status: " << std::to_string(lbbdInfo->SP_status) << std::endl;
        out_txt_f << "OP Status: " << std::to_string(lbbdInfo->OP_status) << std::endl;
        out_txt_f << "Current upper bound: " << std::to_string(lbbdInfo->UB) << std::endl;
        out_txt_f << "Current lower bound: " << std::to_string(lbbdInfo->LB) << std::endl;
        out_txt_f << "Total Time used to solve MP: " << to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
        out_txt_f << "Total Time used to solve SP: " << to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;



        if (lbbdInfo->LBBD_iter == 1)
            out_csv_f2 << std::to_string(lbbdInfo->LBBD_iter) << ',';
        else
            out_csv_f2 << " " << ',' << std::to_string(lbbdInfo->LBBD_iter) << ',';

        out_csv_f2 << std::to_string(lbbdInfo->UB) << ','
            << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(loop_time_mp * 0.001) << ','
            << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->obj_MP) << ',' << std::to_string(lbbdInfo->obj_mp_theta) << ','
            << std::to_string(lbbdInfo->obj_MP - lbbdInfo->obj_mp_theta) << ','
            << std::to_string(lbbdInfo->mp_gap) << ','
            << std::to_string(loop_time_sp * 0.001) << ','
            << std::to_string(lbbdInfo->SP_status) << ','
            << std::to_string(lbbdInfo->obj_SP) << ',';

        out_csv_f2 << std::to_string(spInfo->SP_gap) << ',' << std::to_string(spInfo->num_lazy1) << ','
            << std::to_string(spInfo->num_lazy2) << ',' << std::to_string(spInfo->num_lazy3) << ','
            << std::to_string(spInfo->num_nogood_cut) << ','
            << std::to_string(spInfo->num_sp2_solved) << ',' << std::to_string(spInfo->time_solve_sp2 * 0.001) << ',';
        if (spInfo->num_sp2_solved == 0)
            out_csv_f2 << std::to_string(0) << ',';
        else
        {
            out_csv_f2 << std::to_string(spInfo->time_solve_sp2 * 0.001 / spInfo->num_sp2_solved) << ',';
        }

        out_csv_f2 << std::to_string(spInfo->num_integer_sol) << ',' << std::to_string(spInfo->num_feasible_sol) << ',';
        //if (spInfo->numLazy3List.size() == 0)
        //    out_csv_f2 << std::to_string(0) << ',' << std::to_string(0);
        //else
        //{
        //    std::string ss = "[";
        //    for (int nn = 0; nn < spInfo->numLazy3List.size() - 1; nn++)
        //    {
        //        ss += std::to_string(spInfo->numLazy3List[nn]);
        //        ss += "-";
        //    }
        //    ss += std::to_string(spInfo->numLazy3List.back());
        //    ss += "]";
        //    out_csv_f2 << ss << ',';

        //    std::string ss3 = "[";
        //    for (int nn = 0; nn < spInfo->reduced_numlazy3List.size() - 1; nn++)
        //    {
        //        ss3 += std::to_string(spInfo->reduced_numlazy3List[nn]);
        //        ss3 += "-";
        //    }
        //    ss3 += std::to_string(spInfo->reduced_numlazy3List.back());
        //    ss3 += "]";
        //    out_csv_f2 << ss3;
        //}


        out_csv_f2 << std::endl;
    }
    else
    {
        if (lbbdInfo->UB - lbbdInfo->LB < 1e-3)
        {
            // find the optimal solution
            std::cout << "Get the optimal solution." << std::endl;
            lbbdInfo->lbbd_gap = 0.0;
            lbbdInfo->OP_status = 1;

            out_txt_f << "Get the optimal solution!" << std::endl;
            out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
            out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
            out_txt_f << "Total Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
            out_txt_f << "Total Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;
            out_txt_f << "Avg.Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP * 0.001 / lbbdInfo->LBBD_iter) << std::endl;
            out_txt_f << "Avg.Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->total_SP_solved) << std::endl;

            out_txt_f << "LBBD opt_cut: " << std::to_string(lbbdInfo->lbbd_cut_opt) << std::endl;
            out_txt_f << "LBBD nogood_cut: " << std::to_string(lbbdInfo->lbbd_cut_nogood) << std::endl;

            out_txt_f << "Total Num.SP-Lazy1 cons: " << std::to_string(lbbdInfo->total_lazy1) << std::endl;
            out_txt_f << "Total Num.SP-Lazy2 cons: " << std::to_string(lbbdInfo->total_lazy2) << std::endl;
            out_txt_f << "Total Num.SP-Clique cons: " << std::to_string(lbbdInfo->total_lazy3) << std::endl;
            out_txt_f << "Total Num.SP-No-good cons: " << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << std::endl;

            out_txt_f << "Total Num.SP-S2 solved: " << std::to_string(lbbdInfo->total_num_SP2_solved) << std::endl;
            out_txt_f << "Total Time used to solve SP-S2: " << std::to_string(lbbdInfo->total_time_SP2 * 0.001) << std::endl;
            out_txt_f << "Avg.Time used to solve SP-S2: " << std::to_string(lbbdInfo->total_time_SP2 * 0.001 / lbbdInfo->total_num_SP2_solved) << std::endl;

            out_txt_f << "Avg Num.SP-Lazy1 cons: " << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->LBBD_iter) << std::endl;
            out_txt_f << "Avg Num.SP-Lazy2 cons: " << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->LBBD_iter) << std::endl;
            out_txt_f << "Avg Num.SP-Clique cons: " << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->LBBD_iter) << std::endl;
            out_txt_f << "Avg Num.SP-No-good cons: " << std::to_string(lbbdInfo->total_num_SP2_NogoodCut / lbbdInfo->LBBD_iter) << std::endl;

            out_txt_f << "--------------------------------------------------------" << std::endl;
            out_txt_f << "The body type of each car: " << std::endl;
            for (int bb = 0; bb < Body.size(); bb++) out_txt_f << std::to_string(Body[bb]) << " - ";
            out_txt_f << std::endl;

            out_txt_f << "The color of each car: " << std::endl;
            for (int cc = 0; cc < Color.size(); cc++) out_txt_f << std::to_string(Color[cc]) << " - ";
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
        else
        {
            std::cout << "Do not find the optimal solution because LB != UB." << std::endl;
            lbbdInfo->lbbd_gap = (lbbdInfo->UB - lbbdInfo->LB) / lbbdInfo->LB;

            out_txt_f << "Do not find the optimal solution because LB != UB." << std::endl;
            out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
            out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
            out_txt_f << "MP Status: " << std::to_string(lbbdInfo->MP_status) << std::endl;
            out_txt_f << "SP Status: " << std::to_string(lbbdInfo->SP_status) << std::endl;
            out_txt_f << "OP Status: " << std::to_string(lbbdInfo->OP_status) << std::endl;
            out_txt_f << "Current upper bound: " << std::to_string(lbbdInfo->UB) << std::endl;
            out_txt_f << "Current lower bound: " << std::to_string(lbbdInfo->LB) << std::endl;
            out_txt_f << "Total Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
            out_txt_f << "Total Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;

        }

    }

    // f1: solve by Heu - mpObj - spObj - total time - num.LBBD.Iter - UB - LB
    // OP Status - MP Status - SP Status - Gap - mp.gap - sp.gap - num.opt.cut - num.nogood.cut -
    // total.MP.time - total.SP.time - avg.time.MP - ave.time.SP - num.lazy1 - num.lazy2 - num.lazy3 - num.sp.nogood
    // avg.num.lazy1 - ave.num.lazy2 - avg.num.lazy3 - avg.num.sp.nogood - total.time.sp2 - total.num.sp2 - avg.time.sp2 
    // num.lazy3 added - num.reduced.lazy3.added
    out_csv_f1 << std::to_string(lbbdInfo->Heu_sp) << ',' << std::to_string(lbbdInfo->obj_MP) << ',' << std::to_string(lbbdInfo->obj_SP) << ','
        << std::to_string(lbbdInfo->obj_mp_f1f2) << ',' << std::to_string(lbbdInfo->obj_mp_theta) << ','
        << std::to_string(lbbdInfo->Total_time * 0.001) << ',' << std::to_string(lbbdInfo->LBBD_iter) << ','
        << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
        << std::to_string(lbbdInfo->OP_status) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->SP_status) << ','
        << std::to_string(lbbdInfo->lbbd_gap) << ',' << std::to_string(lbbdInfo->mp_gap) << ',' << std::to_string(lbbdInfo->sp_gap) << ','
        << std::to_string(lbbdInfo->lbbd_cut_opt) << ',' << std::to_string(lbbdInfo->lbbd_cut_nogood) << ','
        << std::to_string(lbbdInfo->total_time_MP * 0.001) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001) << ','
        << std::to_string(lbbdInfo->total_SP_solved) << ','
        << std::to_string(lbbdInfo->total_time_MP * 0.001 / lbbdInfo->LBBD_iter) << ',';
    if (lbbdInfo->total_SP_solved == 0)
        out_csv_f1 << std::to_string(0) << ',';
    else
        out_csv_f1 << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->total_SP_solved) << ',';


    out_csv_f1 << std::to_string(lbbdInfo->total_lazy1) << ',' << std::to_string(lbbdInfo->total_lazy2) << ',' << std::to_string(lbbdInfo->total_lazy3) << ','
        << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << ',';

    if (lbbdInfo->total_SP_solved == 0)
        out_csv_f1 << std::to_string(0) << ',' << std::to_string(0) << ',' << std::to_string(0) << ',' << std::to_string(0) << ',';
    else
        out_csv_f1 << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->LBBD_iter) << ','
        << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->LBBD_iter) << ','
        << std::to_string(lbbdInfo->total_num_SP2_NogoodCut / lbbdInfo->LBBD_iter) << ',';


    out_csv_f1 << std::to_string(lbbdInfo->total_time_SP2 * 0.001) << ',' << std::to_string(lbbdInfo->total_num_SP2_solved) << ',';

    if (lbbdInfo->total_num_SP2_solved == 0)
        out_csv_f1 << std::to_string(0) << ',';
    else
        out_csv_f1 << std::to_string(lbbdInfo->total_time_SP2 * 0.001 / lbbdInfo->total_num_SP2_solved) << ',';


    out_csv_f1 << std::endl;

    out_csv_f1.close();
    out_csv_f2.close();
    out_txt_f.close();


    delete spInfo;

    lbbdInfo->Clear();
    delete lbbdInfo;
    delete OptimalSol;

    if (SP_feasibleSolList.size() != 0)
        for (auto it = SP_feasibleSolList.begin(); it != SP_feasibleSolList.end(); it++)
            for (int i = 0; i < it->second.size(); i++)
                delete it->second[i];
    SP_feasibleSolList.clear();

    return;

}