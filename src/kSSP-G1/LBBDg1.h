#pragma once
#include "HeuInitial.h"
#include "BuildG1.h"
#include "include/easy_digraph.h"

using namespace directed_graph;
using namespace std;

// this code must find the Kahn shortest path
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




// -------------- define LBBD function on G1: use an initial solution for beginning (by Kahn) + K-shortest-path
void LBBDg1(Instance ins, std::string const& ins_name, std::string const& txtfile,
    std::string const& csv1, std::string const& csv2, std::string const& KSP_method,
    bool& addCut1, bool& addCut2, bool& addCut3, bool& LiftCut2, bool& LiftCut3, int& sp2_method,
    bool& Kahn_Heu, bool& other_Heu, double timeLimit)
{
    // ins_name: 10_1 for example
    // txtfile: output the final result
    // csv1: record the final summarized results
    // csv2: record the information for each iteration / shortest-path
    // KSP_method: the method used to find the k-shortest-path
    // spCon1Lift: lifiting SP lazy cons1 or not
    // spCon2Lift: lifiting SP lazy cons2 or not
    // spCon3Lift: lifiting SP lazy cons3 or not
    // sp2_method: use which model to solve sp-s2
    // Kahn_Heu: whether use heuristic method to solve SP for the Kahn-SP
    // other_Heu: whether use heuristic method to solve SP for other SP

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
    // build new G1 and find the first solution of the first level
    CAR_COUNT = ins.cars.size();
    solution_newG1 solution_new_g1;
    Solve(ins.cars, ins.c0, ins.alpha, ins.beta, solution_new_g1);
    std::cout << "Solution of the new g1 size: " << std::to_string(solution_new_g1.solution.size()) << std::endl;

    double loop_time_mp = solution_new_g1.time_find_solution;
    lbbdInfo->total_time_MP += loop_time_mp;
    lbbdInfo->Total_time += lbbdInfo->LB_delta_time + loop_time_mp;


    if (solution_new_g1.finish_build == false)
    {
        // do not finish build the new G1 within the timeLimit, return
        
        lbbdInfo->obj_MP = -1;
        out_txt_f << "Do not find the optimal solution - the first SP status is -1." << std::endl;
        out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
        out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "MP Status: " << std::to_string(lbbdInfo->MP_status) << std::endl;
        out_txt_f << "SP Status: " << std::to_string(lbbdInfo->SP_status) << std::endl;
        out_txt_f << "OP Status: " << std::to_string(lbbdInfo->OP_status) << std::endl;
        out_txt_f << "Current upper bound: X" << std::endl;
        out_txt_f << "Current lower bound: " << std::to_string(lbbdInfo->LB) << std::endl;
        out_txt_f << "Total Time used to solve MP: " << to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
        out_txt_f << "Total Time used to solve SP: " << to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;

        // write csv2
        // f2: Iter th - UB - LB - MP.time - mp.status - mp.obj - mp.gap - SP.time - sp.status - sp.obj - sp.gap - num.lazy1 - num.lazy2 
        // - num.lazy3 - num.sp.nogood
        // num.sp2.solved - time.sp2.solved - avg.time.sp2 - num.integer.found - num.feasible.solution.found - 
        // num.lazy3.added - num.reduced.lazy3.added
        out_csv_f2 << std::to_string(lbbdInfo->LBBD_iter) << ',' << std::to_string(-1) << ',' << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(loop_time_mp * 0.001) << ',' << std::to_string(lbbdInfo->MP_status) << ','
            << std::to_string(lbbdInfo->obj_MP) << ','
            << std::to_string(lbbdInfo->mp_gap) << ','
            << std::to_string(0) << ','
            << std::to_string(lbbdInfo->SP_status) << ',' << std::to_string(-1) << ','
            << std::to_string(lbbdInfo->sp_gap) << ','
            << std::to_string(lbbdInfo->lastSP->num_lazy1) << ',' << std::to_string(lbbdInfo->lastSP->num_lazy2) << ','
            << std::to_string(lbbdInfo->lastSP->num_lazy3) << ',' << std::to_string(lbbdInfo->lastSP->num_nogood_cut) << ','
            << std::to_string(lbbdInfo->lastSP->num_sp2_solved) << ',' << std::to_string(0) << ',';
        if (lbbdInfo->lastSP->num_sp2_solved == 0)
            out_csv_f2 << std::to_string(0) << ',';
        else
            out_csv_f2 << std::to_string(0) << ',';

        out_csv_f2 << std::to_string(0) << ',' << std::to_string(0) << ',';
        out_csv_f2 << std::endl;

        // f1: solve by Heu - Kahn_Heu - Kahn+SP2 - mpObj - spObj - total time - num.LBBD.Iter - UB - LB
        // OP Status - MP Status - SP Status - Gap - mp.gap - sp.gap - 
        // total.MP.time - total.SP.time - avg.time.MP - ave.time.SP - num.lazy1 - num.lazy2 - num.lazy3 - num.sp.nogood
        // avg.num.lazy1 - ave.num.lazy2 - avg.num.lazy3 - avg.num.sp.nogood - total.time.sp2 - total.num.sp2 - avg.time.sp2 
        // num.lazy3 added - num.reduced.lazy3.added
        out_csv_f1 << std::to_string(0) << ',' << std::to_string(0) << ',';
        out_csv_f1 << std::to_string(0) << ',' << std::to_string(0) << ',' << std::to_string(0) << ',';

        out_csv_f1 << std::to_string(lbbdInfo->Heu_sp) << ','
            << std::to_string(lbbdInfo->Kahn_Heu) << ',' << std::to_string(lbbdInfo->Kahn_SP2) << ','
            << std::to_string(lbbdInfo->obj_MP) << ',' << std::to_string(-1) << ','
            << std::to_string(lbbdInfo->Total_time * 0.001) << ',' << std::to_string(lbbdInfo->LBBD_iter) << ','
            << std::to_string(-1) << ',' << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(lbbdInfo->OP_status) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->SP_status) << ','
            << std::to_string(lbbdInfo->lbbd_gap) << ',' << std::to_string(lbbdInfo->mp_gap) << ',' << std::to_string(lbbdInfo->sp_gap) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001) << ','
            << std::to_string(lbbdInfo->total_SP_solved) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001 / lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy1) << ',' << std::to_string(lbbdInfo->total_lazy2) << ',' << std::to_string(lbbdInfo->total_lazy3) << ','
            << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << ',' << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_num_SP2_NogoodCut / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_time_SP2 * 0.001) << ',' << std::to_string(lbbdInfo->total_num_SP2_solved) << ',';

        out_csv_f1 << std::to_string(0) << ',' << std::endl;
        //if (lbbdInfo->total_num_SP2_solved == 0)
        //    out_csv_f1 << std::to_string(0) << ',';
        //else
        //    out_csv_f1 << std::to_string(lbbdInfo->total_time_SP2 * 0.001 / lbbdInfo->total_num_SP2_solved) << ',';
        out_csv_f1.close();
        out_csv_f2.close();
        out_txt_f.close();


        lbbdInfo->Clear();

        delete lbbdInfo;
        delete OptimalSol;

        return;
    }
    
    // get the optimal upstream sequence using the new G1
    std::map<int, std::tuple<int, int>> carId_bodyColor;
    for (int i = 0; i < ins.cars.size(); i++)
        carId_bodyColor[ins.cars[i].id] = std::make_tuple(ins.cars[i].body, ins.cars[i].color);

    int car_n = 0;
    std::map< std::tuple<int, int>, std::vector<int> > V_bm;
    for (int i = 0; i < solution_new_g1.solution.size() - 1; i++)
    {
        int car_id = solution_new_g1.solution[i];
        int car_b = std::get<0>(carId_bodyColor[car_id]);
        int car_c = std::get<1>(carId_bodyColor[car_id]);
        lbbdInfo->lB.push_back(car_b);
        lbbdInfo->lC.push_back(car_c);
        car_n += 1;
        V_bm[std::make_tuple(car_b, car_c)].push_back(car_n);
    }

    lbbdInfo->obj_MP = solution_new_g1.obj_shortest_path;
    lbbdInfo->mp_gap = 0.0;
    lbbdInfo->MP_status = 1;

    double loop_time_sp = 0.0;
    bool is_newg1_Heu_feasible = false;
    bool is_heu_opt = false;
    bool mpsp_status = true;
    lbbdInfo->LB = solution_new_g1.obj_shortest_path + lbbdInfo->LB_delta_obj;
    double remain_time = timeLimit - lbbdInfo->Total_time * 0.001;



    // ------------- solve SP for the shortest path of the new G1
    std::map<int, std::vector<std::tuple<int, int>>> Config_BodyColor;
    StatConBodyColor(ins, Config_BodyColor);

    InformationSolveSP* spInfo = new InformationSolveSP();
    std::map<double, std::vector<SP1_Node_Sol*>> SP_feasibleSolList;
    lbbdInfo->lastSP = spInfo;
    bool use_heu_newG1 = Kahn_Heu;

    // whether use the heurtistic method
    if (use_heu_newG1 == 1)
    {
        // first construct an initial solution by the best sequence of the subproblem / configuration
        auto mp_t1 = std::chrono::steady_clock::now();
        is_newg1_Heu_feasible = HeuDone(env, ins, Config_BodyColor, lbbdInfo->lB, lbbdInfo->lC, *lbbdInfo, sp2_method);
        auto mp_t2 = std::chrono::steady_clock::now();
        auto mp_t = std::chrono::duration_cast<std::chrono::milliseconds>(mp_t2 - mp_t2);   // the time for solving SP by heuristic

        loop_time_sp = mp_t.count();

        if (is_newg1_Heu_feasible == true)
        {
            lbbdInfo->obj_SP = lbbdInfo->LB_delta_obj;
            lbbdInfo->SP_status = 1;
            lbbdInfo->sp_gap = 0.0;

            std::cout << "--------------------- Finish solving the SP of new G1 by Heu. -----------" << std::endl;
            is_heu_opt = true;
        }
        else
        {
            // heu fails, solve SP by B&C
            Solve_SP(*spInfo, SP_feasibleSolList, ins, V_bm, env, addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, remain_time);
            std::cout << "---------- Finsh solving the SP of new G1 by by B&C." << std::endl;

            loop_time_sp += spInfo->SP_total_time;
            lbbdInfo->total_SP_solved += spInfo->num_SP_solved;
            // lbbdInfo->total_time_SP += loop_time_sp;

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
            }
            else
            {
                lbbdInfo->obj_SP = spInfo->SP_OptObjVal;
                lbbdInfo->SP_status = 1;
                // lbbdInfo->Kahn_SP2 = 1;
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
        // solve Kahn-SP by B&B directly
        Solve_SP(*spInfo, SP_feasibleSolList, ins, V_bm, env, addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, remain_time);
        std::cout << "---------- Finsh solving the SP of new G1 by B&C." << std::endl;


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
        }
        else
        {
            lbbdInfo->obj_SP = spInfo->SP_OptObjVal;
            lbbdInfo->SP_status = 1;
            lbbdInfo->Kahn_SP2 = 1;
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
    // stop_time += loop_time_sp;

    remain_time -= loop_time_sp * 0.001;

    if (mpsp_status == false)
    {
        // the first solution to the SP is infeasible
        // output the infeasible solution
        lbbdInfo->obj_SP = -1;
        out_txt_f << "Do not find the optimal solution - the first SP status is -1." << std::endl;
        out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
        out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "MP Status: " << std::to_string(lbbdInfo->MP_status) << std::endl;
        out_txt_f << "SP Status: " << std::to_string(lbbdInfo->SP_status) << std::endl;
        out_txt_f << "OP Status: " << std::to_string(lbbdInfo->OP_status) << std::endl;
        out_txt_f << "Current upper bound: X" << std::endl;
        out_txt_f << "Current lower bound: " << std::to_string(lbbdInfo->LB) << std::endl;
        out_txt_f << "Total Time used to solve MP: " << to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
        out_txt_f << "Total Time used to solve SP: " << to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;

        // write csv2
        // f2: Iter th - UB - LB - MP.time - mp.status - mp.obj - mp.gap - SP.time - sp.status - sp.obj - sp.gap - num.lazy1 - num.lazy2 
        // - num.lazy3 - num.sp.nogood
        // num.sp2.solved - time.sp2.solved - avg.time.sp2 - num.integer.found - num.feasible.solution.found - 
        // num.lazy3.added - num.reduced.lazy3.added
        out_csv_f2 << std::to_string(lbbdInfo->LBBD_iter) << ',' << std::to_string(-1) << ',' << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(loop_time_mp * 0.001) << ',' << std::to_string(lbbdInfo->MP_status) << ','
            << std::to_string(lbbdInfo->obj_MP) << ','
            << std::to_string(lbbdInfo->mp_gap) << ','
            << std::to_string(loop_time_sp * 0.001) << ','
            << std::to_string(lbbdInfo->SP_status) << ',' << std::to_string(-1) << ','
            << std::to_string(lbbdInfo->sp_gap) << ','
            << std::to_string(lbbdInfo->lastSP->num_lazy1) << ',' << std::to_string(lbbdInfo->lastSP->num_lazy2) << ','
            << std::to_string(lbbdInfo->lastSP->num_lazy3) << ',' << std::to_string(lbbdInfo->lastSP->num_nogood_cut) << ','
            << std::to_string(lbbdInfo->lastSP->num_sp2_solved) << ',' << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001) << ',';
        if (lbbdInfo->lastSP->num_sp2_solved == 0)
            out_csv_f2 << std::to_string(0) << ',';
        else
            out_csv_f2 << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001 / lbbdInfo->lastSP->num_sp2_solved) << ',';

        out_csv_f2 << std::to_string(lbbdInfo->lastSP->num_integer_sol) << ',' << std::to_string(lbbdInfo->lastSP->num_feasible_sol) << ',';


        out_csv_f2 << std::endl;

        // f1: solve by Heu - Kahn_Heu - Kahn+SP2 - mpObj - spObj - total time - num.LBBD.Iter - UB - LB
        // OP Status - MP Status - SP Status - Gap - mp.gap - sp.gap - 
        // total.MP.time - total.SP.time - avg.time.MP - ave.time.SP - num.lazy1 - num.lazy2 - num.lazy3 - num.sp.nogood
        // avg.num.lazy1 - ave.num.lazy2 - avg.num.lazy3 - avg.num.sp.nogood - total.time.sp2 - total.num.sp2 - avg.time.sp2 
        // num.lazy3 added - num.reduced.lazy3.added

        out_csv_f1 << std::to_string(0) << ',' << std::to_string(0) << ',';
        out_csv_f1 << std::to_string(0) << ',' << std::to_string(0) << ',' << std::to_string(0) << ',';

        out_csv_f1 << std::to_string(lbbdInfo->Heu_sp) << ','
            << std::to_string(lbbdInfo->Kahn_Heu) << ',' << std::to_string(lbbdInfo->Kahn_SP2) << ','
            << std::to_string(lbbdInfo->obj_MP) << ',' << std::to_string(-1) << ','
            << std::to_string(lbbdInfo->Total_time * 0.001) << ',' << std::to_string(lbbdInfo->LBBD_iter) << ','
            << std::to_string(-1) << ',' << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(lbbdInfo->OP_status) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->SP_status) << ','
            << std::to_string(lbbdInfo->lbbd_gap) << ',' << std::to_string(lbbdInfo->mp_gap) << ',' << std::to_string(lbbdInfo->sp_gap) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001) << ','
            << std::to_string(lbbdInfo->total_SP_solved) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001 / lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy1) << ',' << std::to_string(lbbdInfo->total_lazy2) << ',' << std::to_string(lbbdInfo->total_lazy3) << ','
            << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << ',' << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_num_SP2_NogoodCut / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_time_SP2 * 0.001) << ',' << std::to_string(lbbdInfo->total_num_SP2_solved) << ',';

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

        for (auto it = SP_feasibleSolList.begin(); it != SP_feasibleSolList.end(); it++)
            for (int i = 0; i < it->second.size(); i++)
                delete it->second[i];
        SP_feasibleSolList.clear();

        return;
    }

    // get the optimal solution to the SP
    lbbdInfo->UB = lbbdInfo->obj_MP + lbbdInfo->obj_SP;
    // write csv2
    out_csv_f2 << std::to_string(lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
        << std::to_string(loop_time_mp * 0.001) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->obj_MP) << ','
        << std::to_string(lbbdInfo->mp_gap) << ','
        << std::to_string(loop_time_sp * 0.001) << ',' << std::to_string(lbbdInfo->SP_status) << ',' << std::to_string(lbbdInfo->obj_SP) << ','
        << std::to_string(lbbdInfo->lastSP->SP_gap) << ','
        << std::to_string(lbbdInfo->lastSP->num_lazy1) << ',' << std::to_string(lbbdInfo->lastSP->num_lazy2) << ','
        << std::to_string(lbbdInfo->lastSP->num_lazy3) << ',' << std::to_string(lbbdInfo->lastSP->num_nogood_cut) << ','
        << std::to_string(lbbdInfo->lastSP->num_sp2_solved) << ',' << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001) << ',';
    if (lbbdInfo->lastSP->num_sp2_solved == 0)
        out_csv_f2 << std::to_string(0) << ',';
    else
        out_csv_f2 << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001 / lbbdInfo->lastSP->num_sp2_solved) << ',';


    if (is_heu_opt == true)
        out_csv_f2 << std::to_string(1) << ',' << std::to_string(1) << ',';
    //out_csv_f2 << std::to_string(1) << ',' << std::to_string(1) << ','
    //<< std::to_string(0) << ',' << std::to_string(0) << ',';
    else
    {
        out_csv_f2 << std::to_string(lbbdInfo->lastSP->num_integer_sol) << ',' 
            << std::to_string(lbbdInfo->lastSP->num_feasible_sol) << ',';
    }

    out_csv_f2 << std::endl;
    // whether the solution to the SP of the new G1 is optimal
    if (fabs(lbbdInfo->obj_SP - lbbdInfo->LB_delta_obj) < 1e-2)
    {
        // the solution is the global optimal
        lbbdInfo->OP_status = 1;
        if (is_newg1_Heu_feasible == true)
            lbbdInfo->Heu_sp = 1;
        
        std::cout << "Get the global optimal solution by the shorest path of the new G1" << std::endl;
        std::cout << "Finish Solving " << ins_name << std::endl;
        std::cout << "-----------------------------------------------" << std::endl;
        lbbdInfo->lbbd_gap = 0.0;

        out_txt_f << "Get the global optimal solution using the shorest path of the new G1" << std::endl;
        out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
        out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "Total Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
        out_txt_f << "Total Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;
        out_txt_f << "Avg.Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
        out_txt_f << "Avg.Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;

        out_txt_f << "Total Num.SP-Lazy1 cons: " << std::to_string(lbbdInfo->total_lazy1) << std::endl;
        out_txt_f << "Total Num.SP-Lazy2 cons: " << std::to_string(lbbdInfo->total_lazy2) << std::endl;
        out_txt_f << "Total Num.SP-Clique cons: " << std::to_string(lbbdInfo->total_lazy3) << std::endl;
        out_txt_f << "Total Num.SP-No-good cons: " << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << std::endl;

        out_txt_f << "Total Num.SP-S2 solved: " << std::to_string(lbbdInfo->total_num_SP2_solved) << std::endl;
        out_txt_f << "Total Time used to solve SP-S2: " << std::to_string(lbbdInfo->total_time_SP2 * 0.001) << std::endl;

        if (lbbdInfo->total_num_SP2_solved == 0)
            out_txt_f << "Avg.Time used to solve SP-S2: 0" << std::endl;
        else
            out_txt_f << "Avg.Time used to solve SP-S2: " << std::to_string(lbbdInfo->total_time_SP2 * 0.001 / lbbdInfo->total_num_SP2_solved) << std::endl;

        out_txt_f << "Avg Num.SP-Lazy1 cons: " << std::to_string(lbbdInfo->total_lazy1) << std::endl;
        out_txt_f << "Avg Num.SP-Lazy2 cons: " << std::to_string(lbbdInfo->total_lazy2) << std::endl;
        out_txt_f << "Avg Num.SP-Clique cons: " << std::to_string(lbbdInfo->total_lazy3) << std::endl;
        out_txt_f << "Avg Num.SP-No-good cons: " << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << std::endl;

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

        // f1: solve by Heu - Kahn_Heu - Kahn+SP2 - mpObj - spObj - total time - num.LBBD.Iter - UB - LB
        // OP Status - MP Status - SP Status - Gap - mp.gap - sp.gap - 
        // total.MP.time - total.SP.time - avg.time.MP - ave.time.SP - num.lazy1 - num.lazy2 - num.lazy3 - num.sp.nogood
        // avg.num.lazy1 - ave.num.lazy2 - avg.num.lazy3 - avg.num.sp.nogood - total.time.sp2 - total.num.sp2 - avg.time.sp2 
        // num.lazy3 added - num.reduced.lazy3.added

        out_csv_f1 << std::to_string(0) << ',' << std::to_string(0) << ',';
        out_csv_f1 << std::to_string(0) << ',' << std::to_string(0) << ',' << std::to_string(0) << ',';

        out_csv_f1 << std::to_string(lbbdInfo->Heu_sp) << ','
            << std::to_string(lbbdInfo->Kahn_Heu) << ',' << std::to_string(lbbdInfo->Kahn_SP2) << ','
            << std::to_string(lbbdInfo->obj_MP) << ',' << std::to_string(lbbdInfo->obj_SP) << ','
            << std::to_string(lbbdInfo->Total_time * 0.001) << ',' << std::to_string(lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(lbbdInfo->OP_status) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->SP_status) << ','
            << std::to_string(lbbdInfo->lbbd_gap) << ',' << std::to_string(lbbdInfo->mp_gap) << ',' << std::to_string(lbbdInfo->sp_gap) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001) << ','
            << std::to_string(lbbdInfo->total_SP_solved) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001 / lbbdInfo->LBBD_iter) << ',';

        out_csv_f1 << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->LBBD_iter) << ',';

        out_csv_f1 << std::to_string(lbbdInfo->total_lazy1) << ',' << std::to_string(lbbdInfo->total_lazy2) << ','
            << std::to_string(lbbdInfo->total_lazy3) << ','
            << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << ',';


        out_csv_f1 << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->LBBD_iter) << ','
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


        //g1.Clear();
        delete spInfo;

        lbbdInfo->Clear();
        delete lbbdInfo;
        delete OptimalSol;

        for (auto it = SP_feasibleSolList.begin(); it != SP_feasibleSolList.end(); it++)
            for (int i = 0; i < it->second.size(); i++)
                delete it->second[i];
        SP_feasibleSolList.clear();


        return;

    }



    // ******************************************************************************************************************** 
    // the first shortest path of the new G1 is not the global optimal solution
    // build the original G1
    lbbdInfo->LBBD_iter += 1;
    lbbdInfo->SP_status = -1;
    lbbdInfo->OP_status = -1;
    lbbdInfo->obj_SP = -1;
    lbbdInfo->Heu_sp = 0;
    lbbdInfo->sp_gap = -1;
    lbbdInfo->Kahn_Heu = 0;
    lbbdInfo->Kahn_SP2 = 0;

    loop_time_mp = 0.0;
    loop_time_sp = 0.0;

    lbbdInfo->lB.clear();
    lbbdInfo->lC.clear();
    lbbdInfo->lE.clear();
    lbbdInfo->lO.clear();
    lbbdInfo->lA.clear();
    V_bm.clear();

    std::map< std::tuple<int, int>, std::vector<int> > V_bm3;

    Graph1 g1;
    BuildGraph1(g1, ins, ins_name, remain_time);

    lbbdInfo->graph1 = &g1;


    if (g1.is_built == false)
    {
        std::cout << "Remain time is not enough to build the original G1" << std::endl;
        out_txt_f << "Remain time is not enough to build the original G1" << std::endl;

        loop_time_mp = g1.tt1;
        lbbdInfo->total_time_MP += loop_time_mp;
        lbbdInfo->Total_time += loop_time_mp;
        lbbdInfo->MP_status = -1;
        lbbdInfo->OP_status = -1;
        lbbdInfo->SP_status = -1;
        lbbdInfo->LBBD_iter += 1;
        lbbdInfo->obj_MP = -1;
        lbbdInfo->obj_SP = -1;



        // write file 1 and file 2
        out_csv_f1 << std::to_string(0) << ',' << std::to_string(0) << ',';
        out_csv_f1 << std::to_string(0) << ',' << std::to_string(0) << ',' << std::to_string(0) << ',';

        out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
        out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "MP Status: " << std::to_string(lbbdInfo->MP_status) << std::endl;
        out_txt_f << "SP Status: " << std::to_string(lbbdInfo->SP_status) << std::endl;
        out_txt_f << "OP Status: " << std::to_string(lbbdInfo->OP_status) << std::endl;
        out_txt_f << "Current upper bound: " << std::to_string(lbbdInfo->UB) << std::endl;
        out_txt_f << "Current lower bound: " << std::to_string(lbbdInfo->LB) << std::endl;
        out_txt_f << "Total Time used to solve MP: " << to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
        out_txt_f << "Total Time used to solve SP: " << to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;



        out_csv_f2 << std::to_string(lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(loop_time_mp * 0.001) << ',' << std::to_string(lbbdInfo->MP_status) << ','
            << std::to_string(lbbdInfo->obj_MP) << ','
            << std::to_string(lbbdInfo->mp_gap) << ','
            << std::to_string(loop_time_sp * 0.001) << ','
            << std::to_string(lbbdInfo->SP_status) << ',' << std::to_string(-1) << ','
            << std::to_string(-1) << ','
            << std::to_string(-1) << ',' << std::to_string(-1) << ','
            << std::to_string(-1) << ',' << std::to_string(-1) << ','
            << std::to_string(-1) << ',' << std::to_string(-1) << ',';
        out_csv_f2 << std::to_string(-1) << ',';


        out_csv_f2 << std::to_string(-1) << ',' << std::to_string(-1) << ',';
        out_csv_f2 << std::endl;

        // f1: solve by Heu - Kahn_Heu - Kahn+SP2 - mpObj - spObj - total time - num.LBBD.Iter - UB - LB
        // OP Status - MP Status - SP Status - Gap - mp.gap - sp.gap - 
        // total.MP.time - total.SP.time - avg.time.MP - ave.time.SP - num.lazy1 - num.lazy2 - num.lazy3 - num.sp.nogood
        // avg.num.lazy1 - ave.num.lazy2 - avg.num.lazy3 - avg.num.sp.nogood - total.time.sp2 - total.num.sp2 - avg.time.sp2 
        // num.lazy3 added - num.reduced.lazy3.added
        lbbdInfo->lbbd_gap = (lbbdInfo->UB - lbbdInfo->LB) / lbbdInfo->LB;
        lbbdInfo->Kahn_Heu = 0;
        lbbdInfo->Kahn_SP2 = 0;
        out_csv_f1 << std::to_string(lbbdInfo->Heu_sp) << ','
            << std::to_string(lbbdInfo->Kahn_Heu) << ',' << std::to_string(lbbdInfo->Kahn_SP2) << ','
            << std::to_string(lbbdInfo->obj_MP) << ',' << std::to_string(-1) << ','
            << std::to_string(lbbdInfo->Total_time * 0.001) << ',' << std::to_string(lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(lbbdInfo->OP_status) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->SP_status) << ','
            << std::to_string(lbbdInfo->lbbd_gap) << ',' << std::to_string(-1) << ',' << std::to_string(-1) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001) << ','
            << std::to_string(lbbdInfo->total_SP_solved) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001 / lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy1) << ',' << std::to_string(lbbdInfo->total_lazy2) << ',' << std::to_string(lbbdInfo->total_lazy3) << ','
            << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << ',' << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_num_SP2_NogoodCut / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_time_SP2 * 0.001) << ',' << std::to_string(lbbdInfo->total_num_SP2_solved) << ',';

        if (lbbdInfo->total_num_SP2_solved == 0)
            out_csv_f1 << std::to_string(0) << ',';
        else
            out_csv_f1 << std::to_string(lbbdInfo->total_time_SP2 * 0.001 / lbbdInfo->total_num_SP2_solved) << ',';


        out_csv_f1 << std::endl;

        out_csv_f1.close();
        out_csv_f2.close();
        out_txt_f.close();


        //g1.Clear();
        lbbdInfo->Clear();
        delete lbbdInfo;
        delete OptimalSol;

        return;
    }

    out_txt_f << "Num.Nodes.G1: " << g1.nodes.size() << std::endl;
    out_txt_f << "Num.Arcs.G1: " << g1.arcs.size() << std::endl;
    out_txt_f << "Time.Build.G1: " << g1.tt1 * 0.001 << " seconds" << std::endl;
    out_txt_f << "Time.Find.SP: " << g1.tt2 * 0.001 << " seconds" << std::endl;
    out_txt_f << "-----------------------------------------------------" << std::endl;
    out_csv_f1 << std::to_string(g1.nodes.size()) << ',' << std::to_string(g1.arcs.size()) << ',';
    out_csv_f1 << std::to_string(g1.tt1 * 0.001) << ',' << std::to_string(g1.tt2 * 0.001) << ',' << std::to_string(g1.ShortestPathCost) << ',';


    // get the optimal solution to MP by Kahn-SP
    std::reverse(g1.ShortestPathArc.begin(), g1.ShortestPathArc.end());
    int car_n1 = 0;
    for (int n = 0; n < g1.ShortestPathArc.size(); n++)
    {
        Arc1* ar = g1.ShortestPathArc[n];
        lbbdInfo->lB.push_back(ar->body);
        lbbdInfo->lC.push_back(ar->color);
        car_n1 += 1;
        V_bm3[std::make_tuple(ar->body, ar->color)].push_back(car_n1);
    }

    lbbdInfo->obj_MP = g1.ShortestPathCost;
    lbbdInfo->mp_gap = 0.0;
    lbbdInfo->MP_status = 1;

    std::cout << "---------- Get the optimal solution to Kahn shortest path." << std::endl;
    loop_time_mp = g1.tt2;
    bool is_Kahn_Heu_feasible = false;
    is_heu_opt = false;
    mpsp_status = true;

    lbbdInfo->total_time_MP += loop_time_mp;
    lbbdInfo->Total_time += loop_time_mp;

    if (fabs(g1.ShortestPathCost - solution_new_g1.obj_shortest_path) > 1e-2)
        lbbdInfo->LB = g1.ShortestPathCost + lbbdInfo->LB_delta_obj;

    remain_time -= g1.tt1 * 0.001;
    remain_time -= g1.tt2 * 0.001;

    // solve the SP for the Kahn shortest path
    InformationSolveSP* spInfo2 = new InformationSolveSP();
    std::map<double, std::vector<SP1_Node_Sol*>> SP_feasibleSolList2;
    lbbdInfo->lastSP = spInfo2;

    if (Kahn_Heu == 1)
    {
        // use the heuristic method for the Kahn shorest path
        auto mp_t1 = std::chrono::steady_clock::now();
        is_Kahn_Heu_feasible = HeuDone(env, ins, Config_BodyColor, lbbdInfo->lB, lbbdInfo->lC, *lbbdInfo, sp2_method);
        auto mp_t2 = std::chrono::steady_clock::now();
        auto mp_t = std::chrono::duration_cast<std::chrono::milliseconds>(mp_t2 - mp_t2);   // the time for solving SP by heuristic

        loop_time_sp = mp_t.count();

        if (is_Kahn_Heu_feasible == true)
        {
            //lbbdInfo->total_SP_solved += 1;
            //lbbdInfo->total_time_SP += loop_time_sp;
            lbbdInfo->obj_SP = lbbdInfo->LB_delta_obj;
            //lbbdInfo->total_SP_solved += 1;
            lbbdInfo->SP_status = 1;
            lbbdInfo->sp_gap = 0.0;
            //lbbdInfo->Heu_sp = 1;
            //lbbdInfo->Kahn_Heu = 1;

            std::cout << "------------- Finish solving Kahn-SP by Heu." << std::endl;
            is_heu_opt = true;
        }
        else
        {
            // solve by the B&C
            Solve_SP(*spInfo2, SP_feasibleSolList2, ins, V_bm3, env, addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, remain_time);
            std::cout << "---------- Finsh solving Kahn-SP by B&C." << std::endl;


            loop_time_sp += spInfo2->SP_total_time;
            lbbdInfo->total_SP_solved += spInfo2->num_SP_solved;
            // lbbdInfo->total_time_SP += loop_time_sp;

            lbbdInfo->total_lazy1 += spInfo2->num_lazy1;
            lbbdInfo->total_lazy2 += spInfo2->num_lazy2;
            lbbdInfo->total_lazy3 += spInfo2->num_lazy3;

            lbbdInfo->total_num_SP2_NogoodCut += spInfo2->num_nogood_cut;
            lbbdInfo->total_num_SP2_solved += spInfo2->num_sp2_solved;
            lbbdInfo->total_time_SP2 += spInfo2->time_solve_sp2;


            lbbdInfo->SP_status = spInfo2->SP_status;
            lbbdInfo->sp_gap = spInfo2->SP_gap;

            if (lbbdInfo->SP_status == -1)
            {
                // do not find the optimal solution to SP
                mpsp_status = false;
            }
            else
            {
                lbbdInfo->obj_SP = spInfo2->SP_OptObjVal;
                lbbdInfo->SP_status = 1;
                // lbbdInfo->Kahn_SP2 = 1;
                for (int dd = 0; dd < ins.d_T; dd++)
                {
                    lbbdInfo->lE.push_back(spInfo2->E[dd]);
                    lbbdInfo->lO.push_back(spInfo2->O[dd]);
                    lbbdInfo->lA.push_back(spInfo2->A[dd]);

                }
            }

        }

    }
    else
    {
        // solve by the B&C
        Solve_SP(*spInfo2, SP_feasibleSolList2, ins, V_bm3, env, addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, remain_time);
        std::cout << "---------- Finish solving Kahn-SP by B&C." << std::endl;


        loop_time_sp += spInfo2->SP_total_time;
        lbbdInfo->total_SP_solved += spInfo2->num_SP_solved;
        lbbdInfo->total_lazy1 += spInfo2->num_lazy1;
        lbbdInfo->total_lazy2 += spInfo2->num_lazy2;
        lbbdInfo->total_lazy3 += spInfo2->num_lazy3;

        lbbdInfo->total_num_SP2_NogoodCut += spInfo2->num_nogood_cut;
        lbbdInfo->total_num_SP2_solved += spInfo2->num_sp2_solved;
        lbbdInfo->total_time_SP2 += spInfo2->time_solve_sp2;


        lbbdInfo->SP_status = spInfo2->SP_status;
        lbbdInfo->sp_gap = spInfo2->SP_gap;

        if (lbbdInfo->SP_status == -1)
        {
            // do not find the optimal solution to SP
            mpsp_status = false;
        }
        else
        {
            lbbdInfo->obj_SP = spInfo2->SP_OptObjVal;
            lbbdInfo->SP_status = 1;
            lbbdInfo->Kahn_SP2 = 1;
            for (int dd = 0; dd < ins.d_T; dd++)
            {
                lbbdInfo->lE.push_back(spInfo2->E[dd]);
                lbbdInfo->lO.push_back(spInfo2->O[dd]);
                lbbdInfo->lA.push_back(spInfo2->A[dd]);

            }
        }


    }


    lbbdInfo->total_time_SP += loop_time_sp;
    lbbdInfo->Total_time += loop_time_sp;
    // stop_time += loop_time_sp;

    remain_time -= loop_time_sp * 0.001;


    if (mpsp_status == false)
    {
        // the first solution to SP is infeasible
        // output the infeasible solution
        lbbdInfo->obj_SP = -1;
        out_txt_f << "Do not find the optimal solution - the SP status of the Kahn shortest path is -1." << std::endl;
        out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
        out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "MP Status: " << std::to_string(lbbdInfo->MP_status) << std::endl;
        out_txt_f << "SP Status: " << std::to_string(lbbdInfo->SP_status) << std::endl;
        out_txt_f << "OP Status: " << std::to_string(lbbdInfo->OP_status) << std::endl;
        out_txt_f << "Current upper bound: X" << std::endl;
        out_txt_f << "Current lower bound: " << std::to_string(lbbdInfo->LB) << std::endl;
        out_txt_f << "Total Time used to solve MP: " << to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
        out_txt_f << "Total Time used to solve SP: " << to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;

        // write csv2
        // f2: Iter th - UB - LB - MP.time - mp.status - mp.obj - mp.gap - SP.time - sp.status - sp.obj - sp.gap - num.lazy1 - num.lazy2 
        // - num.lazy3 - num.sp.nogood
        // num.sp2.solved - time.sp2.solved - avg.time.sp2 - num.integer.found - num.feasible.solution.found - 
        // num.lazy3.added - num.reduced.lazy3.added
        out_csv_f2 << std::to_string(lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(loop_time_mp * 0.001) << ',' << std::to_string(lbbdInfo->MP_status) << ','
            << std::to_string(lbbdInfo->obj_MP) << ','
            << std::to_string(lbbdInfo->mp_gap) << ','
            << std::to_string(loop_time_sp * 0.001) << ','
            << std::to_string(lbbdInfo->SP_status) << ',' << std::to_string(-1) << ','
            << std::to_string(lbbdInfo->sp_gap) << ','
            << std::to_string(lbbdInfo->lastSP->num_lazy1) << ',' << std::to_string(lbbdInfo->lastSP->num_lazy2) << ','
            << std::to_string(lbbdInfo->lastSP->num_lazy3) << ',' << std::to_string(lbbdInfo->lastSP->num_nogood_cut) << ','
            << std::to_string(lbbdInfo->lastSP->num_sp2_solved) << ',' << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001) << ',';
        if (lbbdInfo->lastSP->num_sp2_solved == 0)
            out_csv_f2 << std::to_string(0) << ',';
        else
            out_csv_f2 << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001 / lbbdInfo->lastSP->num_sp2_solved) << ',';

        out_csv_f2 << std::to_string(lbbdInfo->lastSP->num_integer_sol) << ',' << std::to_string(lbbdInfo->lastSP->num_feasible_sol) << ',';
        out_csv_f2 << std::endl;

        // f1: solve by Heu - Kahn_Heu - Kahn+SP2 - mpObj - spObj - total time - num.LBBD.Iter - UB - LB
        // OP Status - MP Status - SP Status - Gap - mp.gap - sp.gap - 
        // total.MP.time - total.SP.time - avg.time.MP - ave.time.SP - num.lazy1 - num.lazy2 - num.lazy3 - num.sp.nogood
        // avg.num.lazy1 - ave.num.lazy2 - avg.num.lazy3 - avg.num.sp.nogood - total.time.sp2 - total.num.sp2 - avg.time.sp2 
        // num.lazy3 added - num.reduced.lazy3.added
        lbbdInfo->Kahn_Heu = 0;
        lbbdInfo->Kahn_SP2 = 0;
        out_csv_f1 << std::to_string(lbbdInfo->Heu_sp) << ','
            << std::to_string(lbbdInfo->Kahn_Heu) << ',' << std::to_string(lbbdInfo->Kahn_SP2) << ','
            << std::to_string(lbbdInfo->obj_MP) << ',' << std::to_string(-1) << ','
            << std::to_string(lbbdInfo->Total_time * 0.001) << ',' << std::to_string(lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(lbbdInfo->OP_status) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->SP_status) << ','
            << std::to_string(lbbdInfo->lbbd_gap) << ',' << std::to_string(lbbdInfo->mp_gap) << ',' << std::to_string(lbbdInfo->sp_gap) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001) << ','
            << std::to_string(lbbdInfo->total_SP_solved) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001 / lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy1) << ',' << std::to_string(lbbdInfo->total_lazy2) << ',' << std::to_string(lbbdInfo->total_lazy3) << ','
            << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << ',' << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_num_SP2_NogoodCut / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_time_SP2 * 0.001) << ',' << std::to_string(lbbdInfo->total_num_SP2_solved) << ',';

        if (lbbdInfo->total_num_SP2_solved == 0)
            out_csv_f1 << std::to_string(0) << ',';
        else
            out_csv_f1 << std::to_string(lbbdInfo->total_time_SP2 * 0.001 / lbbdInfo->total_num_SP2_solved) << ',';


        out_csv_f1 << std::endl;

        out_csv_f1.close();
        out_csv_f2.close();
        out_txt_f.close();


        //g1.Clear();
        delete spInfo2;
        lbbdInfo->Clear();

        delete lbbdInfo;
        delete OptimalSol;

        for (auto it = SP_feasibleSolList2.begin(); it != SP_feasibleSolList2.end(); it++)
            for (int i = 0; i < it->second.size(); i++)
                delete it->second[i];
        SP_feasibleSolList2.clear();

        return;

    }

    if (lbbdInfo->UB - (lbbdInfo->obj_MP + lbbdInfo->obj_SP) > 1e-2)
        lbbdInfo->UB = lbbdInfo->obj_MP + lbbdInfo->obj_SP;

    // f2: Iter th - UB - LB - MP.time - mp.status - mp.obj - mp.gap - SP.time - sp.status - sp.obj - sp.gap - num.lazy1 - num.lazy2 
        // - num.lazy3 - num.sp.nogood
        // num.sp2.solved - time.sp2.solved - avg.time.sp2 - num.integer.found - num.feasible.solution.found - 
        // num.lazy3.added - num.reduced.lazy3.added
    out_csv_f2 << std::to_string(lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
        << std::to_string(loop_time_mp * 0.001) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->obj_MP) << ','
        << std::to_string(lbbdInfo->mp_gap) << ','
        << std::to_string(loop_time_sp * 0.001) << ',' << std::to_string(lbbdInfo->SP_status) << ',' << std::to_string(lbbdInfo->obj_SP) << ','
        << std::to_string(lbbdInfo->lastSP->SP_gap) << ','
        << std::to_string(lbbdInfo->lastSP->num_lazy1) << ',' << std::to_string(lbbdInfo->lastSP->num_lazy2) << ','
        << std::to_string(lbbdInfo->lastSP->num_lazy3) << ',' << std::to_string(lbbdInfo->lastSP->num_nogood_cut) << ','
        << std::to_string(lbbdInfo->lastSP->num_sp2_solved) << ',' << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001) << ',';
    if (lbbdInfo->lastSP->num_sp2_solved == 0)
        out_csv_f2 << std::to_string(0) << ',';
    else
        out_csv_f2 << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001 / lbbdInfo->lastSP->num_sp2_solved) << ',';


    if (is_heu_opt == true)
        out_csv_f2 << std::to_string(1) << ',' << std::to_string(1) << ',';
    //out_csv_f2 << std::to_string(1) << ',' << std::to_string(1) << ','
    //<< std::to_string(0) << ',' << std::to_string(0) << ',';
    else
    {
        out_csv_f2 << std::to_string(lbbdInfo->lastSP->num_integer_sol) << ',' << std::to_string(lbbdInfo->lastSP->num_feasible_sol) << ',';
    }
    out_csv_f2 << std::endl;


    // whether the solution to Kahn-SP is optimal
    if (fabs(lbbdInfo->obj_SP - lbbdInfo->LB_delta_obj) < 1e-2 || (lbbdInfo->obj_MP + lbbdInfo->LB_delta_obj - lbbdInfo->UB) > 1e-2)
    {
        // optput the optimal solution with the Kahn-Shortest path
        // is_Kahn = true?
        lbbdInfo->OP_status = 1;
        if (lbbdInfo->Kahn_Heu == true || lbbdInfo->Kahn_SP2 == true)
            lbbdInfo->Heu_sp = 1;

        std::cout << "Get the global optimal solution of the Kahn shortest path." << std::endl;
        std::cout << "Finish Solving " << ins_name << std::endl;
        std::cout << "-----------------------------------------------" << std::endl;
        lbbdInfo->lbbd_gap = 0.0;

        out_txt_f << "Get the globaloptimal solution of the Kahn shorest path!" << std::endl;
        out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
        out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "Total Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
        out_txt_f << "Total Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;
        out_txt_f << "Avg.Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
        out_txt_f << "Avg.Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;

        out_txt_f << "Total Num.SP-Lazy1 cons: " << std::to_string(lbbdInfo->total_lazy1) << std::endl;
        out_txt_f << "Total Num.SP-Lazy2 cons: " << std::to_string(lbbdInfo->total_lazy2) << std::endl;
        out_txt_f << "Total Num.SP-Clique cons: " << std::to_string(lbbdInfo->total_lazy3) << std::endl;
        out_txt_f << "Total Num.SP-No-good cons: " << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << std::endl;

        out_txt_f << "Total Num.SP-S2 solved: " << std::to_string(lbbdInfo->total_num_SP2_solved) << std::endl;
        out_txt_f << "Total Time used to solve SP-S2: " << std::to_string(lbbdInfo->total_time_SP2 * 0.001) << std::endl;

        if (lbbdInfo->total_num_SP2_solved == 0)
            out_txt_f << "Avg.Time used to solve SP-S2: 0" << std::endl;
        else
            out_txt_f << "Avg.Time used to solve SP-S2: " << std::to_string(lbbdInfo->total_time_SP2 * 0.001 / lbbdInfo->total_num_SP2_solved) << std::endl;

        out_txt_f << "Avg Num.SP-Lazy1 cons: " << std::to_string(lbbdInfo->total_lazy1) << std::endl;
        out_txt_f << "Avg Num.SP-Lazy2 cons: " << std::to_string(lbbdInfo->total_lazy2) << std::endl;
        out_txt_f << "Avg Num.SP-Clique cons: " << std::to_string(lbbdInfo->total_lazy3) << std::endl;
        out_txt_f << "Avg Num.SP-No-good cons: " << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << std::endl;

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

        // f1: solve by Heu - Kahn_Heu - Kahn+SP2 - mpObj - spObj - total time - num.LBBD.Iter - UB - LB
        // OP Status - MP Status - SP Status - Gap - mp.gap - sp.gap - 
        // total.MP.time - total.SP.time - avg.time.MP - ave.time.SP - num.lazy1 - num.lazy2 - num.lazy3 - num.sp.nogood
        // avg.num.lazy1 - ave.num.lazy2 - avg.num.lazy3 - avg.num.sp.nogood - total.time.sp2 - total.num.sp2 - avg.time.sp2 
        // num.lazy3 added - num.reduced.lazy3.added
        lbbdInfo->Kahn_Heu = 0;
        lbbdInfo->Kahn_SP2 = 0;
        out_csv_f1 << std::to_string(lbbdInfo->Heu_sp) << ','
            << std::to_string(lbbdInfo->Kahn_Heu) << ',' << std::to_string(lbbdInfo->Kahn_SP2) << ','
            << std::to_string(lbbdInfo->obj_MP) << ',' << std::to_string(lbbdInfo->obj_SP) << ','
            << std::to_string(lbbdInfo->Total_time * 0.001) << ',' << std::to_string(lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(lbbdInfo->OP_status) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->SP_status) << ','
            << std::to_string(lbbdInfo->lbbd_gap) << ',' << std::to_string(lbbdInfo->mp_gap) << ',' << std::to_string(lbbdInfo->sp_gap) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001) << ','
            << std::to_string(lbbdInfo->total_SP_solved) << ','
            << std::to_string(lbbdInfo->total_time_MP * 0.001 / lbbdInfo->LBBD_iter) << ',';

        out_csv_f1 << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->LBBD_iter) << ',';

        //if (lbbdInfo->total_SP_solved == 0)
        //    out_csv_f1 << std::to_string(0) << ',';
        //else
        //    out_csv_f1 << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->total_SP_solved) << ',';


        out_csv_f1 << std::to_string(lbbdInfo->total_lazy1) << ',' << std::to_string(lbbdInfo->total_lazy2) << ','
            << std::to_string(lbbdInfo->total_lazy3) << ','
            << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << ',';


        out_csv_f1 << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->LBBD_iter) << ','
            << std::to_string(lbbdInfo->total_num_SP2_NogoodCut / lbbdInfo->LBBD_iter) << ',';

        //if (lbbdInfo->total_SP_solved == 0)
        //    out_csv_f1 << std::to_string(0) << ',' << std::to_string(0) << ',' << std::to_string(0) << ',' << std::to_string(0) << ',';
        //else
        //    out_csv_f1 << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->total_SP_solved) << ','
        //    << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->total_SP_solved) << ','
        //    << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->total_SP_solved) << ','
        //    << std::to_string(lbbdInfo->total_num_SP2_NogoodCut / lbbdInfo->total_SP_solved) << ',';

        out_csv_f1 << std::to_string(lbbdInfo->total_time_SP2 * 0.001) << ',' << std::to_string(lbbdInfo->total_num_SP2_solved) << ',';

        if (lbbdInfo->total_num_SP2_solved == 0)
            out_csv_f1 << std::to_string(0) << ',';
        else
            out_csv_f1 << std::to_string(lbbdInfo->total_time_SP2 * 0.001 / lbbdInfo->total_num_SP2_solved) << ',';

        out_csv_f1 << std::endl;

        out_csv_f1.close();
        out_csv_f2.close();
        out_txt_f.close();


        //g1.Clear();
        delete spInfo;

        lbbdInfo->Clear();
        delete lbbdInfo;
        delete OptimalSol;

        for (auto it = SP_feasibleSolList.begin(); it != SP_feasibleSolList.end(); it++)
            for (int i = 0; i < it->second.size(); i++)
                delete it->second[i];
        SP_feasibleSolList.clear();


        return;
    }

    // ********************************************************************************************************************
    // use LBBD k-shortest path to solve - already have a solution to MP (by Kahn)
    std::cout << "****************Begin KSP: " << std::endl;

    std::vector<int> k_shortest_path_node;
    std::vector<Arc1*> k_shortest_path_arc;
    int c_n = 0;
    std::map<std::tuple<int, int>, std::vector<int>> V_bm1;

    double dif_mp_time = 0;
    double ori_cum_time = 0;
    double cum_mp_time = 0;

    bool is_other_Heu = false;
    bool is_sp_opt = false;


    bool find_opt = false;
    std::string g_read_file = "../GBS/" + ins_name + ".txt";
    //std::string g_read_file = "C://CppProjects//DeFG1//x64//GBS//10_2.txt";
    EasyDirectedGraph<size_t, uint32_t, uint32_t>* g2 = new EasyDirectedGraph<size_t, uint32_t, uint32_t>(g_read_file);

    auto start = std::chrono::steady_clock::now();
    g2->init_kssp(KSP_method, lbbdInfo->graph1->source_id, lbbdInfo->graph1->terminal_id, 1, false);    // 1 is the version pf KSP method
    // std::vector<uint32_t> res;

    OptimalSol->obj_mp = lbbdInfo->obj_MP;
    OptimalSol->obj_sp = lbbdInfo->obj_SP;
    OptimalSol->obj_opt = lbbdInfo->obj_MP + lbbdInfo->obj_SP;

    InformationSolveSP* spInfo1 = new InformationSolveSP();
    std::map<double, std::vector<SP1_Node_Sol*>> SP_feasibleSolList1;


    do
    {
        auto tt1 = std::chrono::steady_clock::now();

        lbbdInfo->lB.clear();
        lbbdInfo->lC.clear();
        lbbdInfo->lE.clear();
        lbbdInfo->lO.clear();
        lbbdInfo->lA.clear();
        loop_time_mp = 0.0;
        loop_time_sp = 0.0;
        lbbdInfo->SP_status = -1;

        if (k_shortest_path_node.size() != 0)
        {
            k_shortest_path_node.clear();
            k_shortest_path_arc.clear();
        }


        auto sol_path = g2->next_path();
        if (sol_path.second == 0)
            break;


        for (size_t nd = 0; nd < sol_path.first.size(); nd++)
            k_shortest_path_node.push_back(sol_path.first[nd]);



        // judge if this path is the same as the shortest path
        bool isequ = is_equal(g1.ShortestPathNode, k_shortest_path_node);
        if (isequ == true) continue;


        auto stop = std::chrono::steady_clock::now();
        auto ksp_dur = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); // 1/1000 seconds

        cum_mp_time = ksp_dur.count();
        dif_mp_time = cum_mp_time - ori_cum_time;
        ori_cum_time = cum_mp_time;

        loop_time_mp = dif_mp_time;
        lbbdInfo->total_time_MP += loop_time_mp;
        lbbdInfo->Total_time += loop_time_mp;

        remain_time -= loop_time_mp * 0.001;

        lbbdInfo->obj_MP = sol_path.second;
        lbbdInfo->LB = lbbdInfo->obj_MP + lbbdInfo->LB_delta_obj;



        for (size_t nd = 0; nd < k_shortest_path_node.size() - 1; nd++)
        {
            int k_from = k_shortest_path_node[nd];
            int k_to = k_shortest_path_node[nd + 1];
            k_shortest_path_arc.push_back(g1.from_to_arc[std::make_tuple(k_from, k_to)]);
        }

        for (size_t aa = 0; aa < k_shortest_path_arc.size(); aa++)
        {
            Arc1* arr = k_shortest_path_arc[aa];
            lbbdInfo->lB.push_back(arr->body);
            lbbdInfo->lC.push_back(arr->color);
            c_n += 1;
            V_bm1[std::make_tuple(arr->body, arr->color)].push_back(c_n);
        }


        // create a new SP object
        // InformationSolveSP* spInfo1 = new InformationSolveSP();
        // std::map<double, std::vector<SP1_Node_Sol*>> SP_feasibleSolList1;
        is_sp_opt = false;
        lbbdInfo->lastSP = spInfo1;


        // res.push_back(sol_path.second);
        lbbdInfo->LBBD_iter += 1;



        is_other_Heu = false;
        // build a new SP object
        if (other_Heu == true)
        {
            // first use Heu to solve SP
            auto ts1 = std::chrono::steady_clock::now();
            is_other_Heu = HeuDone(env, ins, Config_BodyColor, lbbdInfo->lB, lbbdInfo->lC, *lbbdInfo, sp2_method);
            auto te1 = std::chrono::steady_clock::now();
            auto t1_dur = std::chrono::duration_cast<std::chrono::milliseconds>(te1 - ts1);


            // is_sp_opt = true;
            loop_time_sp = t1_dur.count();
            // lbbdInfo->total_SP_solved += 1;

            if (is_other_Heu == true)
            {
                lbbdInfo->obj_SP = lbbdInfo->LB_delta_obj;
                lbbdInfo->SP_status = 1;
                is_sp_opt = true;
                lbbdInfo->sp_gap = 0.0;
                // lbbdInfo->lastSP = spInfo1;
                lbbdInfo->Total_time += loop_time_sp;
                lbbdInfo->total_time_SP += loop_time_sp;
                lbbdInfo->Heu_sp = 1;

                std::cout << "------------- Finish solving SP by heuristic." << std::endl;
                // is_heu_opt = true;
            }
            else
            {
                // solve SP by B&B
                // loop_sp_time is different
                // lbbdInfo->total_SP_solved += 1;
                Solve_SP(*spInfo1, SP_feasibleSolList1, ins, V_bm1, env, addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, remain_time);
                std::cout << "---------- Finish solving SP by B&C." << std::endl;
                // lbbdInfo->lastSP = spInfo1;


                loop_time_sp += spInfo1->SP_total_time;
                lbbdInfo->total_SP_solved += spInfo1->num_SP_solved;

                lbbdInfo->total_lazy1 += spInfo1->num_lazy1;
                lbbdInfo->total_lazy2 += spInfo1->num_lazy2;
                lbbdInfo->total_lazy3 += spInfo1->num_lazy3;

                lbbdInfo->total_num_SP2_NogoodCut += spInfo1->num_nogood_cut;
                lbbdInfo->total_num_SP2_solved += spInfo1->num_sp2_solved;
                lbbdInfo->total_time_SP2 += spInfo1->time_solve_sp2;


                lbbdInfo->SP_status = spInfo1->SP_status;
                lbbdInfo->sp_gap = spInfo1->SP_gap;

                lbbdInfo->Total_time += loop_time_sp;
                lbbdInfo->total_time_SP += loop_time_sp;

                if (lbbdInfo->SP_status == -1)
                {
                    // do not find the optimal solution to SP
                    mpsp_status = false;
                    // write csv2
                    // f2: Iter th - UB - LB - MP.time - mp.status - mp.obj - mp.gap - SP.time - sp.status - sp.obj - sp.gap - num.lazy1 - num.lazy2 
                    // - num.lazy3 - num.sp.nogood
                    // num.sp2.solved - time.sp2.solved - avg.time.sp2 - num.integer.found - num.feasible.solution.found - 
                    // num.lazy3.added - num.reduced.lazy3.added
                    out_csv_f2 << " " << ',' << std::to_string(lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
                        << std::to_string(loop_time_mp * 0.001) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->obj_MP) << ','
                        << std::to_string(lbbdInfo->mp_gap) << ','
                        << std::to_string(loop_time_sp * 0.001) << ',' << std::to_string(lbbdInfo->SP_status) << ',' << "XX" << ','
                        << std::to_string(lbbdInfo->sp_gap) << ','
                        << std::to_string(lbbdInfo->lastSP->num_lazy1) << ',' << std::to_string(lbbdInfo->lastSP->num_lazy2) << ','
                        << std::to_string(lbbdInfo->lastSP->num_lazy3) << ',' << std::to_string(lbbdInfo->lastSP->num_nogood_cut) << ','
                        << std::to_string(lbbdInfo->lastSP->num_sp2_solved) << ',' << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001) << ',';
                    if (lbbdInfo->lastSP->num_sp2_solved == 0)
                        out_csv_f2 << std::to_string(0) << ',';
                    else
                        out_csv_f2 << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001 / lbbdInfo->lastSP->num_sp2_solved) << ',';

                    out_csv_f2 << std::to_string(lbbdInfo->lastSP->num_integer_sol) << ',' << std::to_string(lbbdInfo->lastSP->num_feasible_sol) << ',';     // integer solution and feasible solution

                    /*std::string ss = "[";
                    for (int nn = 0; nn < lbbdInfo->lastSP->numLazy3List.size() - 1; nn++)
                    {
                        ss += std::to_string(lbbdInfo->lastSP->numLazy3List[nn]);
                        ss += "-";
                    }
                    ss += std::to_string(lbbdInfo->lastSP->numLazy3List.back());
                    ss += "]";
                    out_csv_f2 << ss << ',';

                    std::string ss3 = "[";
                    for (int nn = 0; nn < lbbdInfo->lastSP->reduced_numlazy3List.size() - 1; nn++)
                    {
                        ss3 += std::to_string(lbbdInfo->lastSP->reduced_numlazy3List[nn]);
                        ss3 += "-";
                    }
                    ss3 += std::to_string(lbbdInfo->lastSP->reduced_numlazy3List.back());
                    ss3 += "]";
                    out_csv_f2 << ss3;*/

                    out_csv_f2 << std::endl;
                    // record txt
                    out_txt_f << "Do not find the optimal solution - the ksp SP status is -1." << std::endl;
                    out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
                    out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
                    out_txt_f << "MP Status: " << std::to_string(lbbdInfo->MP_status) << std::endl;
                    out_txt_f << "SP Status: " << std::to_string(lbbdInfo->SP_status) << std::endl;
                    out_txt_f << "OP Status: " << std::to_string(lbbdInfo->OP_status) << std::endl;
                    out_txt_f << "Current upper bound: " << std::to_string(lbbdInfo->UB) << std::endl;
                    out_txt_f << "Current lower bound: " << std::to_string(lbbdInfo->LB) << std::endl;
                    out_txt_f << "Total Time used to solve MP: " << to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
                    out_txt_f << "Total Time used to solve SP: " << to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;


                    break;



                }
                else
                {

                    std::cout << "---------- Finish solving SP by B&C and optimally." << std::endl;
                    lbbdInfo->obj_SP = spInfo1->SP_OptObjVal;
                    lbbdInfo->sp_gap = 0.0;
                    is_sp_opt = true;

                    for (int dd = 0; dd < ins.d_T; dd++)
                    {
                        lbbdInfo->lE.push_back(spInfo1->E[dd]);
                        lbbdInfo->lO.push_back(spInfo1->O[dd]);
                        lbbdInfo->lA.push_back(spInfo1->A[dd]);

                    }
                }
            }




        }
        else
        {
            // solve SP by B&B
            // loop_sp_time is different
            // lbbdInfo->total_SP_solved += 1;
            Solve_SP(*spInfo1, SP_feasibleSolList1, ins, V_bm1, env, addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, remain_time);
            std::cout << "---------- Finish solving SP by B&C." << std::endl;
            // lbbdInfo->lastSP = spInfo1;

            loop_time_sp += spInfo1->SP_total_time;
            lbbdInfo->total_SP_solved += spInfo1->num_SP_solved;
            lbbdInfo->total_lazy1 += spInfo1->num_lazy1;
            lbbdInfo->total_lazy2 += spInfo1->num_lazy2;
            lbbdInfo->total_lazy3 += spInfo1->num_lazy3;

            lbbdInfo->total_num_SP2_NogoodCut += spInfo1->num_nogood_cut;
            lbbdInfo->total_num_SP2_solved += spInfo1->num_sp2_solved;
            lbbdInfo->total_time_SP2 += spInfo1->time_solve_sp2;


            lbbdInfo->SP_status = spInfo1->SP_status;
            lbbdInfo->sp_gap = spInfo1->SP_gap;

            lbbdInfo->Total_time += loop_time_sp;
            lbbdInfo->total_time_SP += loop_time_sp;

            if (lbbdInfo->SP_status == -1)
            {
                // do not find the optimal solution to SP
                mpsp_status = false;
                lbbdInfo->total_SP_solved += 1;
                // write csv2
                // f2: Iter th - UB - LB - MP.time - mp.status - mp.obj - mp.gap - SP.time - sp.status - sp.obj - sp.gap - num.lazy1 - num.lazy2 
                // - num.lazy3 - num.sp.nogood
                // num.sp2.solved - time.sp2.solved - avg.time.sp2 - num.integer.found - num.feasible.solution.found - 
                // num.lazy3.added - num.reduced.lazy3.added
                out_csv_f2 << " " << ',' << std::to_string(lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
                    << std::to_string(loop_time_mp * 0.001) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->obj_MP) << ','
                    << std::to_string(lbbdInfo->mp_gap) << ','
                    << std::to_string(loop_time_sp * 0.001) << ',' << std::to_string(lbbdInfo->SP_status) << ',' << std::to_string(-1) << ','
                    << std::to_string(lbbdInfo->sp_gap) << ','
                    << std::to_string(lbbdInfo->lastSP->num_lazy1) << ',' << std::to_string(lbbdInfo->lastSP->num_lazy2) << ','
                    << std::to_string(lbbdInfo->lastSP->num_lazy3) << ',' << std::to_string(lbbdInfo->lastSP->num_nogood_cut) << ','
                    << std::to_string(lbbdInfo->lastSP->num_sp2_solved) << ',' << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001) << ',';
                if (lbbdInfo->lastSP->num_sp2_solved == 0)
                    out_csv_f2 << std::to_string(0) << ',';
                else
                    out_csv_f2 << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001 / lbbdInfo->lastSP->num_sp2_solved) << ',';

                out_csv_f2 << std::to_string(lbbdInfo->lastSP->num_integer_sol) << ',' << std::to_string(lbbdInfo->lastSP->num_feasible_sol) << ',';     // integer solution and feasible solution


                out_csv_f2 << std::endl;
                // record txt
                out_txt_f << "Do not find the optimal solution - the ksp SP status is -1." << std::endl;
                out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
                out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
                out_txt_f << "MP Status: " << std::to_string(lbbdInfo->MP_status) << std::endl;
                out_txt_f << "SP Status: " << std::to_string(lbbdInfo->SP_status) << std::endl;
                out_txt_f << "OP Status: " << std::to_string(lbbdInfo->OP_status) << std::endl;
                out_txt_f << "Current upper bound: " << std::to_string(lbbdInfo->UB) << std::endl;
                out_txt_f << "Current lower bound: " << std::to_string(lbbdInfo->LB) << std::endl;
                out_txt_f << "Total Time used to solve MP: " << to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
                out_txt_f << "Total Time used to solve SP: " << to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;

                break;
            }
            else
            {
                lbbdInfo->obj_SP = spInfo1->SP_OptObjVal;
                lbbdInfo->sp_gap = 0.0;
                is_sp_opt = true;

                for (int dd = 0; dd < ins.d_T; dd++)
                {
                    lbbdInfo->lE.push_back(spInfo1->E[dd]);
                    lbbdInfo->lO.push_back(spInfo1->O[dd]);
                    lbbdInfo->lA.push_back(spInfo1->A[dd]);

                }
            }

        }




        remain_time -= loop_time_sp * 0.001;


        if (lbbdInfo->UB - (lbbdInfo->obj_MP + lbbdInfo->obj_SP) > 1e-2)
        {
            lbbdInfo->UB = lbbdInfo->obj_MP + lbbdInfo->obj_SP;
            OptimalSol->obj_mp = lbbdInfo->obj_MP;
            OptimalSol->obj_sp = lbbdInfo->obj_SP;
            OptimalSol->obj_opt = lbbdInfo->obj_MP + lbbdInfo->obj_SP;

            if (OptimalSol->bodyList.size() != 0)
            {
                OptimalSol->bodyList.clear();
                OptimalSol->colorList.clear();
                OptimalSol->configList.clear();
                OptimalSol->downPosList.clear();
                OptimalSol->laneList.clear();
            }
            for (int opti = 0; opti < lbbdInfo->lB.size(); opti++)
            {
                OptimalSol->bodyList.push_back(lbbdInfo->lB[opti]);
                OptimalSol->colorList.push_back(lbbdInfo->lC[opti]);
                OptimalSol->configList.push_back(lbbdInfo->lO[opti]);
                OptimalSol->downPosList.push_back(lbbdInfo->lE[opti]);
                OptimalSol->laneList.push_back(lbbdInfo->lA[opti]);
            }

        }


        lbbdInfo->LB = lbbdInfo->obj_MP + lbbdInfo->LB_delta_obj;

        // record csv2
        out_csv_f2 << " " << ',' << std::to_string(lbbdInfo->LBBD_iter) << ',' << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
            << std::to_string(loop_time_mp * 0.001) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->obj_MP) << ','
            << std::to_string(lbbdInfo->mp_gap) << ','
            << std::to_string(loop_time_sp * 0.001) << ',' << std::to_string(is_sp_opt) << ',' << std::to_string(lbbdInfo->obj_SP) << ','
            << std::to_string(lbbdInfo->sp_gap) << ','
            << std::to_string(lbbdInfo->lastSP->num_lazy1) << ',' << std::to_string(lbbdInfo->lastSP->num_lazy2) << ','
            << std::to_string(lbbdInfo->lastSP->num_lazy3) << ',' << std::to_string(lbbdInfo->lastSP->num_nogood_cut) << ','
            << std::to_string(lbbdInfo->lastSP->num_sp2_solved) << ',' << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001) << ',';
        if (lbbdInfo->lastSP->num_sp2_solved == 0)
            out_csv_f2 << std::to_string(0) << ',';
        else
            out_csv_f2 << std::to_string(lbbdInfo->lastSP->time_solve_sp2 * 0.001 / lbbdInfo->lastSP->num_sp2_solved) << ',';


        if (is_other_Heu == true)
            out_csv_f2 << std::to_string(1) << ',' << std::to_string(1) << ',';
        //out_csv_f2 << std::to_string(1) << ',' << std::to_string(1) << ','
        //<< std::to_string(0) << ',' << std::to_string(0) << ',';
        else
        {
            out_csv_f2 << std::to_string(lbbdInfo->lastSP->num_integer_sol) << ',' << std::to_string(lbbdInfo->lastSP->num_feasible_sol) << ',';

        }


        out_csv_f2 << std::endl;




        if (fabs(lbbdInfo->obj_SP - lbbdInfo->LB_delta_obj) <= 1e-2 || (lbbdInfo->obj_MP + lbbdInfo->LB_delta_obj - lbbdInfo->UB) > 1e-2)
        {
            // find the optimal solution
            lbbdInfo->SP_status = 1;
            lbbdInfo->OP_status = 1;
            lbbdInfo->lbbd_gap = 0.0;
            find_opt = true;
            delete spInfo1;

            if (is_other_Heu == true)
                lbbdInfo->Heu_sp = true;

            break;

        }
        else
        {
            k_shortest_path_node.clear();
            k_shortest_path_arc.clear();
            c_n = 0;
            V_bm1.clear();

            for (auto it = SP_feasibleSolList1.begin(); it != SP_feasibleSolList1.end(); it++)
                for (int inode = 0; inode < it->second.size(); inode++)
                    delete it->second[inode];
            SP_feasibleSolList1.clear();

            spInfo1->spInfoUpdate();

            lbbdInfo->lbbd_gap = (lbbdInfo->UB - lbbdInfo->LB) / lbbdInfo->LB;

        }

        auto tt2 = std::chrono::steady_clock::now();
        auto tt_dur = std::chrono::duration_cast<std::chrono::milliseconds>(tt2 - tt1);

        // stop_time += tt_dur.count();

        // std::cout << "Opt: " << OptimalSol->obj_mp << " " << OptimalSol->obj_opt << " " << OptimalSol->obj_sp << std::endl;




    } while (remain_time >= 1.0);




    if (find_opt == false)
    {
        // do not get the optimal solution within K-SP
        // write csv1 and txt
        out_txt_f << "Do not find the optimal solution within k iterations." << std::endl;
        out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
        out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "MP Status: " << std::to_string(lbbdInfo->MP_status) << std::endl;
        out_txt_f << "SP Status: " << std::to_string(lbbdInfo->SP_status) << std::endl;
        out_txt_f << "OP Status: " << std::to_string(lbbdInfo->OP_status) << std::endl;
        out_txt_f << "Current upper bound: " << std::to_string(lbbdInfo->UB) << std::endl;
        out_txt_f << "Current lower bound: " << std::to_string(lbbdInfo->LB) << std::endl;
        out_txt_f << "Total Time used to solve MP: " << to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
        out_txt_f << "Total Time used to solve SP: " << to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;
    }
    else
    {
        // write optimal solution
        // write csv1 and txt
        out_txt_f << "Get the global optimal solution in the KSPp!" << std::endl;
        out_txt_f << "Total time: " << std::to_string(lbbdInfo->Total_time * 0.001) << std::endl;
        out_txt_f << "Num.LBBD.Iterations: " << std::to_string(lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "Total Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP * 0.001) << std::endl;
        out_txt_f << "Total Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP * 0.001) << std::endl;
        out_txt_f << "Avg.Time used to solve MP: " << std::to_string(lbbdInfo->total_time_MP * 0.001 / lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "Avg.Time used to solve SP: " << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->LBBD_iter) << std::endl;


        out_txt_f << "Total Num.SP-Lazy1 cons: " << std::to_string(lbbdInfo->total_lazy1) << std::endl;
        out_txt_f << "Total Num.SP-Lazy2 cons: " << std::to_string(lbbdInfo->total_lazy2) << std::endl;
        out_txt_f << "Total Num.SP-Clique cons: " << std::to_string(lbbdInfo->total_lazy3) << std::endl;
        out_txt_f << "Total Num.SP-No-good cons: " << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << std::endl;

        out_txt_f << "Total Num.SP-S2 solved: " << std::to_string(lbbdInfo->total_num_SP2_solved) << std::endl;
        out_txt_f << "Total Time used to solve SP-S2: " << std::to_string(lbbdInfo->total_time_SP2 * 0.001) << std::endl;

        if (lbbdInfo->total_num_SP2_solved == 0)
            out_txt_f << "Avg.Time used to solve SP-S2: " << std::to_string(0) << std::endl;
        else
            out_txt_f << "Avg.Time used to solve SP-S2: " << std::to_string(lbbdInfo->total_time_SP2 * 0.001 / lbbdInfo->total_num_SP2_solved) << std::endl;

        out_txt_f << "Avg Num.SP-Lazy1 cons: " << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "Avg Num.SP-Lazy2 cons: " << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "Avg Num.SP-Clique cons: " << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->LBBD_iter) << std::endl;
        out_txt_f << "Avg Num.SP-No-good cons: " << std::to_string(lbbdInfo->total_num_SP2_NogoodCut / lbbdInfo->LBBD_iter) << std::endl;

        out_txt_f << "--------------------------------------------------------" << std::endl;
        out_txt_f << "The body type of each car: " << std::endl;
        for (int bb = 0; bb < OptimalSol->bodyList.size(); bb++) out_txt_f << std::to_string(OptimalSol->bodyList[bb]) << " - ";
        out_txt_f << std::endl;

        out_txt_f << "The color of each car: " << std::endl;
        for (int cc = 0; cc < OptimalSol->colorList.size(); cc++) out_txt_f << std::to_string(OptimalSol->colorList[cc]) << " - ";
        out_txt_f << std::endl;

        out_txt_f << "The configuration of each car: " << std::endl;
        for (int oo = 0; oo < OptimalSol->configList.size(); oo++) out_txt_f << std::to_string(OptimalSol->configList[oo]) << " - ";
        out_txt_f << std::endl;

        out_txt_f << "The lane to which each car is allcoated: " << std::endl;
        for (int ll = 0; ll < OptimalSol->laneList.size(); ll++) out_txt_f << std::to_string(OptimalSol->laneList[ll]) << " - ";
        out_txt_f << std::endl;

        out_txt_f << "The upstream car_id in each downstream position: " << std::endl;
        for (int ee = 0; ee < OptimalSol->downPosList.size(); ee++) out_txt_f << std::to_string(OptimalSol->downPosList[ee]) << " - ";
        out_txt_f << std::endl;

    }

    // f1: solve by Heu - Kahn_Heu - Kahn+SP2 - mpObj - spObj - total time - num.LBBD.Iter - UB - LB
       // OP Status - MP Status - SP Status - Gap - mp.gap - sp.gap - 
       // total.MP.time - total.SP.time - avg.time.MP - ave.time.SP - num.lazy1 - num.lazy2 - num.lazy3 - num.sp.nogood
       // avg.num.lazy1 - ave.num.lazy2 - avg.num.lazy3 - avg.num.sp.nogood - total.time.sp2 - total.num.sp2 - avg.time.sp2 
       // num.lazy3 added - num.reduced.lazy3.added


    out_csv_f1 << std::to_string(lbbdInfo->Heu_sp) << ','
        << std::to_string(lbbdInfo->Kahn_Heu) << ',' << std::to_string(lbbdInfo->Kahn_SP2) << ','
        << std::to_string(OptimalSol->obj_mp) << ',' << std::to_string(OptimalSol->obj_sp) << ','
        << std::to_string(lbbdInfo->Total_time * 0.001) << ',' << std::to_string(lbbdInfo->LBBD_iter) << ','
        << std::to_string(lbbdInfo->UB) << ',' << std::to_string(lbbdInfo->LB) << ','
        << std::to_string(lbbdInfo->OP_status) << ',' << std::to_string(lbbdInfo->MP_status) << ',' << std::to_string(lbbdInfo->SP_status) << ','
        << std::to_string(lbbdInfo->lbbd_gap) << ',' << std::to_string(lbbdInfo->mp_gap) << ',' << std::to_string(lbbdInfo->sp_gap) << ','
        << std::to_string(lbbdInfo->total_time_MP * 0.001) << ',' << std::to_string(lbbdInfo->total_time_SP * 0.001) << ','
        << std::to_string(lbbdInfo->total_SP_solved) << ','
        << std::to_string(lbbdInfo->total_time_MP * 0.001 / lbbdInfo->LBBD_iter) << ',';

    out_csv_f1 << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->LBBD_iter) << ',';

    //if (lbbdInfo->total_SP_solved == 0)
    //    out_csv_f1 << std::to_string(0) << ',';
    //else
    //    out_csv_f1 << std::to_string(lbbdInfo->total_time_SP * 0.001 / lbbdInfo->total_SP_solved) << ',';


    out_csv_f1 << std::to_string(lbbdInfo->total_lazy1) << ',' << std::to_string(lbbdInfo->total_lazy2) << ','
        << std::to_string(lbbdInfo->total_lazy3) << ','
        << std::to_string(lbbdInfo->total_num_SP2_NogoodCut) << ',';


    out_csv_f1 << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->LBBD_iter) << ','
        << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->LBBD_iter) << ','
        << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->LBBD_iter) << ','
        << std::to_string(lbbdInfo->total_num_SP2_NogoodCut / lbbdInfo->LBBD_iter) << ',';

    //if (lbbdInfo->total_SP_solved == 0)
    //    out_csv_f1 << std::to_string(0) << ',' << std::to_string(0) << ',' << std::to_string(0) << ',' << std::to_string(0) << ',';
    //else
    //    out_csv_f1 << std::to_string(lbbdInfo->total_lazy1 / lbbdInfo->total_SP_solved) << ','
    //    << std::to_string(lbbdInfo->total_lazy2 / lbbdInfo->total_SP_solved) << ','
    //    << std::to_string(lbbdInfo->total_lazy3 / lbbdInfo->total_SP_solved) << ','
    //    << std::to_string(lbbdInfo->total_num_SP2_NogoodCut / lbbdInfo->total_SP_solved) << ',';

    out_csv_f1 << std::to_string(lbbdInfo->total_time_SP2 * 0.001) << ',' << std::to_string(lbbdInfo->total_num_SP2_solved) << ',';

    if (lbbdInfo->total_num_SP2_solved == 0)
        out_csv_f1 << std::to_string(0) << ',';
    else
        out_csv_f1 << std::to_string(lbbdInfo->total_time_SP2 * 0.001 / lbbdInfo->total_num_SP2_solved) << ',';

    out_csv_f1 << std::endl;


    std::cout << "Finish Solving " << ins_name << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;


    out_csv_f1.close();
    out_csv_f2.close();
    out_txt_f.close();


    //g1.Clear();
    delete spInfo;

    lbbdInfo->Clear();
    delete lbbdInfo;
    delete OptimalSol;

    if (SP_feasibleSolList.size() != 0)
        for (auto it = SP_feasibleSolList.begin(); it != SP_feasibleSolList.end(); it++)
            for (int i = 0; i < it->second.size(); i++)
                delete it->second[i];
    SP_feasibleSolList.clear();


    if (k_shortest_path_arc.size() != 0) k_shortest_path_arc.clear();
    if (k_shortest_path_node.size() != 0) k_shortest_path_node.clear();



    if (SP_feasibleSolList1.size() != 0)
        for (auto it = SP_feasibleSolList1.begin(); it != SP_feasibleSolList1.end(); it++)
            if (it->second.size() != 0)
                for (int i = 0; i < it->second.size(); i++)
                    delete it->second[i];
    SP_feasibleSolList1.clear();


    return;






}