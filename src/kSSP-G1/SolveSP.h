#pragma once
#include "gurobi_c++.h"
#include "LiftCuts.h"
#include "FindCliques.h"
#include "Elements.h"
#include <chrono>
#include <string>


using namespace std;

class SP1_Node_Sol
{
public:
    double obj_val;
    std::vector<int> O;
    std::vector<int> E;
    std::vector<int> A;

    SP1_Node_Sol() : obj_val(0.0) {}

};


struct SP2_Sol
{
    int sp2_status = -1;				// the status got by solving sp-s2
    double sp2_solve_time = 0.0;		// the time for each sp-s2 solving
    std::vector<int> A;			        // the allocation between car and lane
    int num_nogood_cut = 0;			    // if sp2-status is -1, the number of no good cuts added to sp-s1
    std::vector< std::vector<std::tuple<int, int>> > sp2_nogood_set;		// the items involved in the no-good cut
};


SP2_Sol Solve_SP2(GRBEnv& env, std::map<std::tuple<int, int>, double>& S1_sol, int num_car, int num_l, int q_0, int method)
{
    // method1: a feasible problem
    // method2: a maximization problem with untouched |L|
    // method3: a minimization problem with |L|+1

    SP2_Sol sp2_solution;

    std::vector<std::tuple<int, int>> E_pos;
    for (auto s11 = S1_sol.begin(); s11 != S1_sol.end(); s11++)
        if (s11->second > 0.5)
            E_pos.push_back(s11->first);

    // get parameters z
    std::map<std::tuple<int, int>, int> z;
    for (int i = 1; i <= num_car; i++)
        for (int p = 1; p <= num_car; p++)
            z[std::make_tuple(i, p)] = 0;
    int sb = 1;
    int cn = 0;
    for (int sp = 1; sp <= num_car; sp++)
    {
        double sum_e1 = 0;
        for (int i = 1; i <= sp; i++)
            sum_e1 += S1_sol[std::make_tuple(i, sb)];

        if (sum_e1 < 1e-4)
        {
            sb = sb;
            cn = cn;

            if (cn < 1e-4)
                for (int p = 1; p < sp + 1; p++) z[std::make_tuple(sp, p)] = 0;
            else
            {
                for (int p = 1; p < cn + 1; p++) z[std::make_tuple(sp, p)] = 1;
                for (int p = cn + 1; p < sp + 1; p++) z[std::make_tuple(sp, p)] = 0;
            }
        }
        else
        {
            sb = sb + 1;
            cn = cn + 1;

            for (int p = 1; p < cn + 1; p++) z[std::make_tuple(sp, p)] = 1;
            for (int p = cn + 1; p < sp + 1; p++) z[std::make_tuple(sp, p)] = 0;
        }
    }


    // build the 2nd stage model
    switch (method)
    {
    case 1:
    {
        auto sp2_t1 = std::chrono::steady_clock::now();


        GRBModel SP_S2_model = GRBModel(env);
        SP_S2_model.set(GRB_IntParam_OutputFlag, 0);
        //SP_S2_model.set(GRB_IntParam_Threads, 8);

        // add variables
        std::vector< std::vector<GRBVar> > A(num_car);
        for (int i = 0; i < num_car; i++)
        {
            std::vector<GRBVar> row(num_l);
            for (int l = 0; l < num_l; l++)
            {
                std::string var_name = "a[" + std::to_string(i + 1) + "][" + std::to_string(l + 1) + "]";
                GRBVar a = SP_S2_model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                row[l] = a;
            }
            A[i] = row;
        }

        // add constraints
        for (int i = 0; i < num_car; i++)
        {
            GRBLinExpr sp2_con1 = 0;
            for (int l = 0; l < num_l; l++)
                sp2_con1 += A[i][l];
            SP_S2_model.addConstr(sp2_con1, GRB_EQUAL, 1.0);
        }

        for (int i2 = 2; i2 < num_car + 1; i2++)
            for (int i1 = 1; i1 < i2; i1++)
                for (int l = 1; l <= num_l; l++)
                    for (int p1 = 1; p1 <= num_car; p1++)
                    {
                        GRBLinExpr sp2_con2 = 0;
                        for (int p = p1 + 1; p < num_car + 1; p++) sp2_con2 += S1_sol[std::make_tuple(i1, p)];
                        for (int p = 1; p < p1 + 1; p++) sp2_con2 += S1_sol[std::make_tuple(i2, p)];
                        sp2_con2 += A[i2 - 1][l - 1];
                        sp2_con2 += A[i1 - 1][l - 1];
                        SP_S2_model.addConstr(sp2_con2, GRB_LESS_EQUAL, 3.1);
                    }
        for (int i = 1; i <= num_car; i++)
            for (int l = 1; l <= num_l; l++)
            {
                GRBLinExpr sp2_con3 = 0;
                for (int i_ = 1; i_ < i + 1; i_++) sp2_con3 += A[i_ - 1][l - 1];
                for (int p = 1; p <= i; p++)
                    for (int k = 1; k <= i; k++)
                        sp2_con3 -= z[std::make_tuple(i, p)] * S1_sol[std::make_tuple(k, p)] * A[k - 1][l - 1];
                SP_S2_model.addConstr(sp2_con3, GRB_LESS_EQUAL, q_0 + 0.05);
            }

        // set objective function
        GRBLinExpr sp2_obj = 1.0;
        SP_S2_model.setObjective(sp2_obj, GRB_MINIMIZE);
        SP_S2_model.optimize();

        // std::cout << "******sp2 status " << SP_S2_model.get(GRB_IntAttr_Status) << std::endl;

        if (SP_S2_model.get(GRB_IntAttr_Status) == 2)
        {
            auto sp2_t2 = std::chrono::steady_clock::now();
            auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(sp2_t2 - sp2_t1);
            sp2_solution.sp2_solve_time = duration1.count();


            sp2_solution.num_nogood_cut = 0;
            sp2_solution.sp2_status = 1;

            for (int i = 0; i < num_car; i++)
                for (int l = 0; l < num_l; l++)
                    if (A[i][l].get(GRB_DoubleAttr_X) > 0.05)
                        sp2_solution.A.push_back(l + 1);
        }
        else
        {
            auto sp2_t2 = std::chrono::steady_clock::now();
            auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(sp2_t2 - sp2_t1);
            sp2_solution.sp2_solve_time = duration1.count();

            sp2_solution.sp2_status = -1;
            sp2_solution.num_nogood_cut = 1;

            std::vector<std::tuple<int, int>> no_good_s;
            for (int i = 1; i <= num_car; i++)
                for (int p = 1; p <= num_car; p++)
                    if (S1_sol[std::make_tuple(i, p)] > 0.5) no_good_s.push_back(std::make_tuple(i, p));
            sp2_solution.sp2_nogood_set.push_back(no_good_s);

        }
        break;
    }

    case 2:
    {
        auto sp2_t1 = std::chrono::steady_clock::now();

        GRBModel SP_S2_model = GRBModel(env);
        SP_S2_model.set(GRB_IntParam_OutputFlag, 0);
        //SP_S2_model.set(GRB_IntParam_Threads, 8);

        // add variables
        std::vector< std::vector<GRBVar> > A(num_car);
        for (int i = 0; i < num_car; i++)
        {
            std::vector<GRBVar> row(num_l);
            for (int l = 0; l < num_l; l++)
            {
                std::string var_name = "a[" + std::to_string(i + 1) + "][" + std::to_string(l + 1) + "]";
                GRBVar a = SP_S2_model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                row[l] = a;
            }
            A[i] = row;
        }

        // add constraints
        for (int i = 0; i < num_car; i++)
        {
            GRBLinExpr sp2_con1 = 0;
            for (int l = 0; l < num_l; l++)
                sp2_con1 += A[i][l];
            SP_S2_model.addConstr(sp2_con1, GRB_LESS_EQUAL, 1.0);
        }

        for (int i2 = 2; i2 < num_car + 1; i2++)
            for (int i1 = 1; i1 < i2; i1++)
                for (int l = 1; l <= num_l; l++)
                    for (int p1 = 1; p1 <= num_car; p1++)
                    {
                        GRBLinExpr sp2_con2 = 0;
                        for (int p = p1 + 1; p < num_car + 1; p++) sp2_con2 += S1_sol[std::make_tuple(i1, p)];
                        for (int p = 1; p < p1 + 1; p++) sp2_con2 += S1_sol[std::make_tuple(i2, p)];
                        sp2_con2 += A[i2 - 1][l - 1];
                        sp2_con2 += A[i1 - 1][l - 1];
                        SP_S2_model.addConstr(sp2_con2, GRB_LESS_EQUAL, 3.1);
                    }
        for (int i = 1; i <= num_car; i++)
            for (int l = 1; l <= num_l; l++)
            {
                GRBLinExpr sp2_con3 = 0;
                for (int i_ = 1; i_ < i + 1; i_++) sp2_con3 += A[i_ - 1][l - 1];
                for (int p = 1; p <= i; p++)
                    for (int k = 1; k <= i; k++)
                        sp2_con3 -= z[std::make_tuple(i, p)] * S1_sol[std::make_tuple(k, p)] * A[k - 1][l - 1];
                SP_S2_model.addConstr(sp2_con3, GRB_LESS_EQUAL, q_0 + 0.05);
            }

        // set objective function
        GRBLinExpr sp2_obj = 0;
        for (int i = 0; i < num_car; i++)
            for (int l = 0; l < num_l; l++)
                sp2_obj += A[i][l];

        SP_S2_model.setObjective(sp2_obj, GRB_MAXIMIZE);

        SP_S2_model.optimize();


        auto sp2_t2 = std::chrono::steady_clock::now();
        auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(sp2_t2 - sp2_t1);
        sp2_solution.sp2_solve_time = duration1.count();



        std::map<int, int> E_pos;
        for (int i = 1; i <= num_car; i++)
            for (int p = 1; p <= num_car; p++)
                if (S1_sol[std::make_tuple(i, p)] > 0.5)
                {
                    E_pos[i] = p;
                }


        if (SP_S2_model.get(GRB_IntAttr_Status) == 2)
        {
            if (SP_S2_model.get(GRB_DoubleAttr_ObjVal) >= double(num_car) - 0.1)
            {

                sp2_solution.sp2_status = 1;
                sp2_solution.num_nogood_cut = 0;
                for (int i = 0; i < num_car; i++)
                    for (int l = 0; l < num_l; l++)
                        if (A[i][l].get(GRB_DoubleAttr_X) > 0.5)
                            sp2_solution.A.push_back(l + 1);

            }
            else
            {
                sp2_solution.sp2_status = -1;

                std::vector<std::tuple<int, int>> car_allo;
                std::vector<std::tuple<int, int>> car_unallo;


                for (int i = 0; i < num_car; i++)
                {
                    double car_v = 0;
                    for (int l = 0; l < num_l; l++) car_v += A[i][l].get(GRB_DoubleAttr_X);

                    if (car_v > 0.5)
                        car_allo.push_back(std::make_tuple(i + 1, E_pos[i + 1]));
                    else
                        car_unallo.push_back(std::make_tuple(i + 1, E_pos[i + 1]));
                }

                sp2_solution.num_nogood_cut = car_unallo.size();
                std::vector<std::tuple<int, int>> new_no_good;
                for (int no = 0; no < car_unallo.size(); no++)
                {
                    new_no_good.assign(car_allo.begin(), car_allo.end());
                    new_no_good.push_back(car_unallo[no]);
                    sp2_solution.sp2_nogood_set.push_back(new_no_good);
                    new_no_good.clear();
                }
                std::cout << std::endl;
            }
        }
        else
        {
            sp2_solution.sp2_status = -1;
            sp2_solution.num_nogood_cut = 1;
            std::vector<std::tuple<int, int>> no_good_s;
            for (int i = 1; i <= num_car; i++)
                for (int p = 1; p <= num_car; p++)
                    if (S1_sol[std::make_tuple(i, p)] > 0.5) no_good_s.push_back(std::make_tuple(i, p));
            sp2_solution.sp2_nogood_set.push_back(no_good_s);
        }
        break;
    }

    case 3:
    {
        auto sp2_t1 = std::chrono::steady_clock::now();

        GRBModel SP_S2_model = GRBModel(env);
        SP_S2_model.set(GRB_IntParam_OutputFlag, 0);
        //SP_S2_model.set(GRB_IntParam_Threads, 8);

        // add variables
        std::vector< std::vector<GRBVar> > A(num_car);
        for (int i = 0; i < num_car; i++)
        {
            std::vector<GRBVar> row(num_l + 1);
            for (int l = 0; l < num_l + 1; l++)
            {
                std::string var_name = "a[" + std::to_string(i + 1) + "][" + std::to_string(l + 1) + "]";
                GRBVar a = SP_S2_model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                row[l] = a;
            }
            A[i] = row;
        }

        // add constraints
        for (int i = 0; i < num_car; i++)
        {
            GRBLinExpr sp2_con1 = 0;
            for (int l = 0; l < num_l + 1; l++)
                sp2_con1 += A[i][l];
            SP_S2_model.addConstr(sp2_con1, GRB_EQUAL, 1.0);
        }

        for (int i2 = 2; i2 < num_car + 1; i2++)
            for (int i1 = 1; i1 < i2; i1++)
                for (int l = 1; l <= num_l; l++)
                    for (int p1 = 1; p1 <= num_car; p1++)
                    {
                        GRBLinExpr sp2_con2 = 0;
                        for (int p = p1 + 1; p < num_car + 1; p++) sp2_con2 += S1_sol[std::make_tuple(i1, p)];
                        for (int p = 1; p < p1 + 1; p++) sp2_con2 += S1_sol[std::make_tuple(i2, p)];
                        sp2_con2 += A[i2 - 1][l - 1];
                        sp2_con2 += A[i1 - 1][l - 1];
                        SP_S2_model.addConstr(sp2_con2, GRB_LESS_EQUAL, 3.1);
                    }
        for (int i = 1; i <= num_car; i++)
            for (int l = 1; l <= num_l; l++)
            {
                GRBLinExpr sp2_con3 = 0;
                for (int i_ = 1; i_ < i + 1; i_++) sp2_con3 += A[i_ - 1][l - 1];
                for (int p = 1; p <= i; p++)
                    for (int k = 1; k <= i; k++)
                        sp2_con3 -= z[std::make_tuple(i, p)] * S1_sol[std::make_tuple(k, p)] * A[k - 1][l - 1];
                SP_S2_model.addConstr(sp2_con3, GRB_LESS_EQUAL, q_0 + 0.05);
            }

        // set objective function
        GRBLinExpr sp2_obj = 0;
        for (int i = 0; i < num_car; i++)
            sp2_obj += A[i][num_l];

        SP_S2_model.setObjective(sp2_obj, GRB_MINIMIZE);

        SP_S2_model.optimize();


        auto sp2_t2 = std::chrono::steady_clock::now();
        auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(sp2_t2 - sp2_t1);
        sp2_solution.sp2_solve_time = duration1.count();


        std::map<int, int> E_pos;
        for (int i = 1; i <= num_car; i++)
            for (int p = 1; p <= num_car; p++)
                if (S1_sol[std::make_tuple(i, p)] > 0.5) E_pos[i] = p;

        if (SP_S2_model.get(GRB_IntAttr_Status) == 2)
        {
            if (SP_S2_model.get(GRB_DoubleAttr_ObjVal) < 0.5)
            {

                sp2_solution.sp2_status = 1;
                sp2_solution.num_nogood_cut = 0;
                for (int i = 0; i < num_car; i++)
                    for (int l = 0; l < num_l; l++)
                        if (A[i][l].get(GRB_DoubleAttr_X) > 0.5)
                            sp2_solution.A.push_back(l + 1);

            }
            else if (SP_S2_model.get(GRB_DoubleAttr_ObjVal) > 0.5)
            {
                sp2_solution.sp2_status = -1;

                std::vector<std::tuple<int, int>> car_allo;
                std::vector<std::tuple<int, int>> car_unallo;

                for (int i = 0; i < num_car; i++)
                {
                    if (A[i][num_l].get(GRB_DoubleAttr_X) > 0.5)
                        car_unallo.push_back(std::make_tuple(i + 1, E_pos[i + 1]));
                    else
                        car_allo.push_back(std::make_tuple(i + 1, E_pos[i + 1]));
                }

                sp2_solution.num_nogood_cut = car_unallo.size();

                std::vector<std::tuple<int, int>> new_no_good;
                for (int no = 0; no < car_unallo.size(); no++)
                {
                    new_no_good.assign(car_allo.begin(), car_allo.end());
                    new_no_good.push_back(car_unallo[no]);
                    sp2_solution.sp2_nogood_set.push_back(new_no_good);
                    new_no_good.clear();
                }
            }
        }
        else
        {
            sp2_solution.sp2_status = -1;
            sp2_solution.num_nogood_cut = 1;
            std::vector<std::tuple<int, int>> no_good_s;
            for (int i = 1; i <= num_car; i++)
                for (int p = 1; p <= num_car; p++)
                    if (S1_sol[std::make_tuple(i, p)] > 0.5) no_good_s.push_back(std::make_tuple(i, p));
            sp2_solution.sp2_nogood_set.push_back(no_good_s);
        }
        break;
    }
    default:
        break;
    }


    return sp2_solution;
}




// define callback function
class SP_Callback : public GRBCallback
{
public:
    Instance Ins;
    std::map<std::tuple<int, int, int>, GRBVar> vars;
    std::map<double, std::vector<SP1_Node_Sol*>>* feasibleNodeList;


    GRBEnv* Env;
    bool addCut1;
    bool addCut2;
    bool addCut3;
    bool liftCut2;
    bool liftCut3;

    int addCut3_method;     // 0: add all cliques; 1: add cliques until all positions occpuied, 2: add no more than k cliques
    int k_lic3;             // number of lazy3 added
    int sp2method;

    InformationSolveSP* OptSPInformation;


    SP_Callback(GRBEnv* env, Instance ins, std::map<std::tuple<int, int, int>, GRBVar> Evars,
        std::map<double, std::vector<SP1_Node_Sol*>>& feasibleSolList, bool addL1, bool addL2, bool addL3, InformationSolveSP& SP_SolvedInfo,
        bool liftc2, bool liftc3, int lic3_m, int k_lic3_num, int sp2m)
    {
        vars = Evars;
        Ins = ins;
        Env = env;
        feasibleNodeList = &feasibleSolList;
        addCut1 = addL1;
        addCut2 = addL2;
        addCut3 = addL3;
        OptSPInformation = &SP_SolvedInfo;
        liftCut2 = liftc2;
        liftCut3 = liftc3;
        addCut3_method = lic3_m;
        k_lic3 = k_lic3_num;
        sp2method = sp2m;
    }


    void addLazyCut1(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& Posk2,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt);
    void addLazyCut2(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt,
        bool& liftCut2);
    void addLazyCut3(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt,
        bool& liftCut3, int& how_to_addCut3, int& num_Cut3_added, std::vector<std::tuple<int, int>>& all_pairs);

    void ifAddNogoodCut(GRBEnv& env, Instance& ins, std::map<std::tuple<int, int>, double>& sp_sol, int& sp2_approach, std::map<double,
        std::vector<SP1_Node_Sol*>>*feasibleSolList, double& sp_obj);



    void K2(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt);

    void K3(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt);

    void K4(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt);

    // Denest: add Cut1 
    void FindAndAdd_Cut1(int& num_j, std::vector<int>& C1_set, std::vector<std::tuple<int, int>>& car_pos,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt);
    // K=2
    void K2_Find_i1(int& ik, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt);
    void newK2(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt);

    // K=3
    void K3_Find_i1(int& i2, int& p2, int& ik, int& pk, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt, bool& add_cut1);

    void K3_Find_i2(int& ik, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt, bool& add_cut1);

    void newK3(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt);

    // K=4
    void K4_Find_i1(int& i2, int& p2, int& i3, int& p3, int& ik, int& pk, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt, bool& add_cut1);

    void K4_Find_i2(int& i3, int& p3, int& ik, int& pk, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt, bool& add_cut1);

    void K4_Find_i3(int& ik, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt, bool& add_cut1);

    void newK4(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
        std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt);


protected:
    void callback()
    {
        try
        {
            if (where == GRB_CB_MIPSOL)
            {
                std::cout << "Enter callback: " << std::endl;
                // std::cout << "Num SolCount: " << getIntInfo(GRB_CB_MIPSOL_SOLCNT) << std::endl;
                OptSPInformation->num_callback_cut = 0;
                OptSPInformation->num_integer_sol += 1;
                // std::cout << "Num.integer solution: " << OptSPInformation->num_integer_sol << std::endl;


                std::vector<std::tuple<int, int>> down_seq(Ins.d_T);    // (i,p): car i in downstream position p
                std::vector<std::tuple<int, int>> all_matches;          // all (i,p), |V|*|V|
                std::map<std::tuple<int, int>, double> S1_sol;          // (i,p): 0/1

                std::vector<int> PosK2;
                //std::vector<int> PosK3;
                //std::vector<int> PosK4;

                OptSPInformation->num_integer_sol += 1;
                double obj = getDoubleInfo(GRB_CB_MIPSOL_OBJ);

                //std::cout << "Sol: " << std::endl;

                for (int i = 1; i <= Ins.d_T; i++)
                    for (int p = 1; p <= Ins.d_T; p++)
                    {

                        all_matches.push_back(std::make_tuple(i, p));
                        double ev = 0;

                        for (int j = 1; j <= Ins.num_j; j++)
                            ev += getSolution(vars[std::make_tuple(i, p, j)]);

                        S1_sol[std::make_tuple(i, p)] = ev;

                        if (ev > 0.5)
                        {
                            //std::cout << "(" << i << "," << p << ")-";
                            down_seq[i - 1] = std::make_tuple(i, p);

                            if (i > Ins.q0 * (Ins.num_l - 1) && i - p >= Ins.q0 * (Ins.num_l - 1))
                            {
                                PosK2.push_back(i);     // i is ik in newk2

                                //if (i > 3 + (Ins.num_l - 3) * Ins.q0)
                                //    PosK3.push_back(i);     // i is ik in newk3

                                //if (i > 4 + (Ins.num_l - 4) * Ins.q0)
                                //    PosK4.push_back(i);     // i is ik in newk4
                            }

                        }
                    }


                // ---------------  add lazycon1 or not
                if (addCut1 == true && Ins.num_l > 2 && (Ins.d_T > (Ins.num_l - 1) * Ins.q0))
                    addLazyCut1(Ins, down_seq, PosK2, vars, OptSPInformation);


                // ---------------  add lazycon2 or not
                if (addCut2 == true && OptSPInformation->num_callback_cut == 0)
                    addLazyCut2(Ins, down_seq, vars, OptSPInformation, liftCut2);

                // ---------------  add lazycon3 or not (maximal cliques)
                if (addCut3 == true && OptSPInformation->num_callback_cut == 0)
                    addLazyCut3(Ins, down_seq, vars, OptSPInformation, liftCut3, addCut3_method, k_lic3, all_matches);


                // ---------------  whether to add no good cuts
                if (OptSPInformation->num_callback_cut == 0)
                    ifAddNogoodCut(*Env, Ins, S1_sol, sp2method, feasibleNodeList, obj);



            } // find an integer solution
        }
        catch (GRBException e) {
            std::cout << "Error number: " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        }

    }
};



void SP_Callback::addLazyCut1(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& Posk2,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt)
{
    switch (ins.num_l)
    {
    case 3:
    {
        // K2(ins, car_pos, Posk2, var, Opt);
        newK2(ins, car_pos, Posk2, var, Opt);
        break;
    }
    case 4:
    {

        // K2(ins, car_pos, Posk2, var, Opt);
        newK2(ins, car_pos, Posk2, var, Opt);
        //if (Opt->num_callback_cut == 0)
            // K3(ins, car_pos, Posk2, var, Opt);
        break;
    }
    case 5:
    {

        // K2(ins, car_pos, Posk2, var, Opt);
        newK2(ins, car_pos, Posk2, var, Opt);

        //if (Opt->num_callback_cut == 0)
            // K3(ins, car_pos, Posk2, var, Opt);
            //newK3(ins, car_pos, Posk2, var, Opt);

        //if (Opt->num_callback_cut == 0)
            // K4(ins, car_pos, Posk2, var, Opt);
            //newK4(ins, car_pos, Posk2, var, Opt);

        break;
    }
    default:
        break;
    }


}


void SP_Callback::addLazyCut2(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::map<std::tuple<int, int, int>, GRBVar>& var,
    InformationSolveSP* Opt, bool& liftCut2)
{
    for (int k1 = 1; k1 <= ins.d_T; k1++)
    {
        int p1 = std::get<1>(car_pos[k1 - 1]);
        std::vector<std::tuple<int, int>> C1;   // k < k1, p > p1
        std::vector<std::tuple<int, int>> C2;   // k > k1, p < p1

        int C1_size = 0;
        int C2_size = 0;

        if (k1 == 1)
        {
            // no C1
            // generate C2
            for (int k = k1 + 1; k <= Ins.d_T; k++)
                for (int p = 1; p <= p1; p++)
                    C2.push_back(std::make_tuple(k, p));

            //std::cout << "C2: ";
            for (int k = k1 + 1; k <= Ins.d_T; k++)
            {
                int p = std::get<1>(car_pos[k - 1]);
                if (p1 > p)
                {
                    C2_size += 1;
                }
            }

        }
        else
        {
            // generate C1
            for (int k = 1; k < k1; k++)
                for (int p = p1; p <= Ins.d_T; p++)
                    C1.push_back(std::make_tuple(k, p));
            //std::cout << "C1: ";
            for (int k = 1; k < k1; k++)
            {
                int p = std::get<1>(car_pos[k - 1]);
                if (p > p1)
                {
                    C1_size += 1;
                }
            }

            // generate C2
            for (int k = k1 + 1; k <= Ins.d_T; k++)
                for (int p = 1; p <= p1; p++)
                    C2.push_back(std::make_tuple(k, p));

            for (int k = k1 + 1; k <= Ins.d_T; k++)
            {
                int p = std::get<1>(car_pos[k - 1]);
                if (p1 > p)
                {
                    C2_size += 1;
                }

            }
        }

        if (Ins.num_l >= 2)
        {
            if (C2_size > 0 && C1_size > (Ins.num_l - 2) * Ins.q0)
            {
                OptSPInformation->num_callback_cut += 1;
                OptSPInformation->num_lazy2 += 1;

                if (liftCut2 == false)
                {
                    // no lifting on cut 2
                    GRBLinExpr lazy2 = 0;
                    for (int j = 1; j <= Ins.num_j; j++)
                    {
                        lazy2 += vars[make_tuple(k1, p1, j)];
                    }


                    for (int c1 = 0; c1 < C1.size(); c1++)
                    {
                        for (int j = 1; j <= Ins.num_j; j++)
                        {
                            lazy2 += vars[make_tuple(std::get<0>(C1[c1]), std::get<1>(C1[c1]), j)];
                        }
                    }


                    for (int c2 = 0; c2 < C2.size(); c2++)
                    {
                        for (int j = 1; j <= Ins.num_j; j++)
                        {
                            lazy2 += vars[make_tuple(std::get<0>(C2[c2]), std::get<1>(C2[c2]), j)];
                        }
                    }

                    addLazy(lazy2, GRB_LESS_EQUAL, C1_size + C2_size + 0.01);
                    std::cout << "Add lazy2: " << OptSPInformation->num_lazy2 << std::endl;

                }
                else
                {
                    // lifting on cut 2
                    std::vector<std::tuple<int, int>> C2_Set;
                    C2_Set.push_back(std::make_tuple(k1, p1));
                    for (int c1 = 0; c1 < C1.size(); c1++)
                        C2_Set.push_back(C1[c1]);
                    for (int c2 = 0; c2 < C2.size(); c2++)
                        C2_Set.push_back(C2[c2]);

                    LiftCon2(C2_Set, Ins.d_T);

                    GRBLinExpr lazy2 = 0;
                    for (int cc = 0; cc < C2_Set.size(); cc++)
                        for (int j = 1; j <= Ins.num_j; j++)
                        {
                            lazy2 += vars[std::make_tuple(std::get<0>(C2_Set[cc]), std::get<1>(C2_Set[cc]), j)];
                        }

                    addLazy(lazy2, GRB_LESS_EQUAL, C1_size + C2_size + 0.01);
                    std::cout << "Add lazy2: " << OptSPInformation->num_lazy2 << std::endl;
                }
            }
        }   // |L| >= 2
        else
        {
            // |L| = 1
            if (C2_size > 0 || C1_size > 0)
            {
                OptSPInformation->num_callback_cut += 1;
                OptSPInformation->num_lazy2 += 1;

                if (liftCut2 == false)
                {
                    // no lifting on cut 2
                    GRBLinExpr lazy2 = 0;
                    for (int j = 1; j <= Ins.num_j; j++)
                    {
                        lazy2 += vars[make_tuple(k1, p1, j)];
                    }
                    for (int c1 = 0; c1 < C1.size(); c1++)
                    {
                        for (int j = 1; j <= Ins.num_j; j++)
                        {
                            lazy2 += vars[make_tuple(std::get<0>(C1[c1]), std::get<1>(C1[c1]), j)];
                        }
                    }
                    for (int c2 = 0; c2 < C2.size(); c2++)
                    {
                        for (int j = 1; j <= Ins.num_j; j++)
                        {
                            lazy2 += vars[make_tuple(std::get<0>(C2[c2]), std::get<1>(C2[c2]), j)];
                        }
                    }
                    addLazy(lazy2, GRB_LESS_EQUAL, C1_size + C2_size);
                    std::cout << "Add lazy2: " << OptSPInformation->num_lazy2 << std::endl;

                }
                else
                {
                    // lifting on lazy cut 2
                    std::vector<std::tuple<int, int>> C2_Set;
                    C2_Set.push_back(std::make_tuple(k1, p1));
                    for (int c1 = 0; c1 < C1.size(); c1++)
                        C2_Set.push_back(C1[c1]);
                    for (int c2 = 0; c2 < C2.size(); c2++)
                        C2_Set.push_back(C2[c2]);

                    LiftCon2(C2_Set, Ins.d_T);

                    GRBLinExpr lazy2 = 0;
                    for (int cc = 0; cc < C2_Set.size(); cc++)
                        for (int j = 1; j <= Ins.num_j; j++)
                        {
                            lazy2 += vars[std::make_tuple(std::get<0>(C2_Set[cc]), std::get<1>(C2_Set[cc]), j)];
                        }

                    addLazy(lazy2, GRB_LESS_EQUAL, C1_size + C2_size);
                    std::cout << "Add lazy2: " << OptSPInformation->num_lazy2 << std::endl;
                }
            }
        }   // |L| = 1

    } // finish adding lazy cut 2

}

void SP_Callback::addLazyCut3(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::map<std::tuple<int, int, int>, GRBVar>& var,
    InformationSolveSP* Opt, bool& liftCut3, int& how_to_addCut3, int& num_Cut3_added, std::vector<std::tuple<int, int>>& all_pairs)
{
    std::vector< std::tuple<int, int> > adj;
    for (int i1 = 1; i1 < Ins.d_T; i1++)
        for (int i2 = i1 + 1; i2 <= Ins.d_T; i2++)
        {
            int p1 = std::get<1>(car_pos[i1 - 1]);
            int p2 = std::get<1>(car_pos[i2 - 1]);
            if (p1 > p2) adj.push_back(make_tuple(i1, i2));
        }

    BK cliques(Ins.d_T, adj);
    cliques.run();

    std::vector< std::vector<std::tuple<int, int>> > invalid_cliques;
    if (cliques.maxCliques.size() != 0)
        for (auto it = cliques.maxCliques.begin(); it != cliques.maxCliques.end(); it++)
        {
            if ((*it).size() > Ins.num_l + 0.01)
            {
                std::vector<std::tuple<int, int>> a_inclique;
                for (auto it2 = (*it).begin(); it2 != (*it).end(); it2++)       //  it is a clique
                {
                    a_inclique.push_back(std::make_tuple((*it2), std::get<1>(car_pos[(*it2) - 1])));
                }
                invalid_cliques.push_back(a_inclique);
            }
        }

    if (invalid_cliques.size() > 0.1)
    {
        switch (how_to_addCut3)
        {
        case 0:
        {
            // ******************add all maximal clique constraints
            OptSPInformation->numLazy3List.push_back(invalid_cliques.size());
            OptSPInformation->reduced_numlazy3List.push_back(invalid_cliques.size());
            for (int s = 0; s < invalid_cliques.size(); s++)
            {
                OptSPInformation->num_lazy3 += 1;
                OptSPInformation->num_callback_cut += 1;
                if (liftCut3 == false)
                {
                    // no lifting
                    GRBLinExpr lazy3 = 0;
                    for (int si = 0; si < invalid_cliques[s].size(); si++)
                    {
                        for (int j = 1; j <= Ins.num_j; j++)
                            lazy3 += vars[std::make_tuple(get<0>(invalid_cliques[s][si]), get<1>(invalid_cliques[s][si]), j)];
                    }

                    addLazy(lazy3, GRB_LESS_EQUAL, Ins.num_l);
                    std::cout << "Add lazy3: " << OptSPInformation->num_lazy3 << std::endl;
                }
                else
                {
                    // lifting con3
                    std::vector<std::tuple<int, int>> liftC3;
                    GRBLinExpr lazy3 = 0;
                    liftCon3(liftC3, invalid_cliques[s], all_pairs);

                    for (int si = 0; si < liftC3.size(); si++)
                    {
                        for (int j = 1; j <= Ins.num_j; j++)
                        {
                            lazy3 += vars[std::make_tuple(std::get<0>(liftC3[si]), std::get<1>(liftC3[si]), j)];
                        }
                    }

                    addLazy(lazy3, GRB_LESS_EQUAL, Ins.num_l);
                    std::cout << "Add lazy3: " << OptSPInformation->num_lazy3 << std::endl;
                }
            }

            break;
        }
        case 1:
        {
            // *******************pick some maximal clique constraint to add(by position occupied)
            // (*OptSPInformation).numLazy3List.push_back(invalid_cliques.size());
            OptSPInformation->numLazy3List.push_back(invalid_cliques.size());
            std::map<int, int> pos_occupied;
            int reduce_num_lazy3 = 0;

            for (int s = 0; s < invalid_cliques.size(); s++)
            {
                OptSPInformation->num_lazy3 += 1;
                OptSPInformation->num_callback_cut += 1;
                reduce_num_lazy3 += 1;

                if (liftCut3 == false)
                {
                    // no lifting
                    GRBLinExpr lazy3 = 0;
                    for (int si = 0; si < invalid_cliques[s].size(); si++)
                    {
                        if (pos_occupied[std::get<1>(invalid_cliques[s][si])] == 0)
                            pos_occupied[std::get<1>(invalid_cliques[s][si])] = 1;

                        for (int j = 1; j <= Ins.num_j; j++)
                            lazy3 += vars[std::make_tuple(get<0>(invalid_cliques[s][si]), get<1>(invalid_cliques[s][si]), j)];
                    }
                    addLazy(lazy3, GRB_LESS_EQUAL, Ins.num_l);
                    std::cout << "Add lazy3: " << OptSPInformation->num_lazy3 << std::endl;
                }
                else
                {
                    // lifting con3
                    std::vector<std::tuple<int, int>> liftC3;
                    GRBLinExpr lazy3 = 0;
                    liftCon3(liftC3, invalid_cliques[s], all_pairs);

                    for (int si = 0; si < invalid_cliques[s].size(); si++)
                    {
                        if (pos_occupied[std::get<1>(invalid_cliques[s][si])] == 0)
                            pos_occupied[std::get<1>(invalid_cliques[s][si])] = 1;

                        for (int j = 1; j <= Ins.num_j; j++)
                            lazy3 += vars[std::make_tuple(get<0>(invalid_cliques[s][si]), get<1>(invalid_cliques[s][si]), j)];
                    }
                    addLazy(lazy3, GRB_LESS_EQUAL, Ins.num_l);
                    std::cout << "Add lazy3: " << OptSPInformation->num_lazy3 << std::endl;
                }

                if (pos_occupied.size() == Ins.d_T)
                    break;
            }
            OptSPInformation->reduced_numlazy3List.push_back(reduce_num_lazy3);

            break;
        }
        case 2:
        {
            OptSPInformation->numLazy3List.push_back(invalid_cliques.size());
            int reduce_num_lazy3 = 0;

            for (int s = 0; s < invalid_cliques.size(); s++)
            {
                OptSPInformation->num_lazy3 += 1;
                OptSPInformation->num_callback_cut += 1;
                reduce_num_lazy3 += 1;

                if (liftCut3 == false)
                {
                    // no lifting
                    GRBLinExpr lazy3 = 0;
                    for (int si = 0; si < invalid_cliques[s].size(); si++)
                    {
                        for (int j = 1; j <= Ins.num_j; j++)
                            lazy3 += vars[std::make_tuple(get<0>(invalid_cliques[s][si]), get<1>(invalid_cliques[s][si]), j)];
                    }
                    addLazy(lazy3, GRB_LESS_EQUAL, Ins.num_l);
                    std::cout << "Add lazy3: " << OptSPInformation->num_lazy3 << std::endl;
                }
                else
                {
                    // lifting con3
                    std::vector<std::tuple<int, int>> liftC3;
                    GRBLinExpr lazy3 = 0;
                    liftCon3(liftC3, invalid_cliques[s], all_pairs);

                    for (int si = 0; si < invalid_cliques[s].size(); si++)
                    {
                        for (int j = 1; j <= Ins.num_j; j++)
                            lazy3 += vars[std::make_tuple(get<0>(invalid_cliques[s][si]), get<1>(invalid_cliques[s][si]), j)];
                    }
                    addLazy(lazy3, GRB_LESS_EQUAL, Ins.num_l);
                    std::cout << "Add lazy3: " << OptSPInformation->num_lazy3 << std::endl;
                }

                if (reduce_num_lazy3 >= k_lic3)
                    break;
            }
            OptSPInformation->reduced_numlazy3List.push_back(reduce_num_lazy3);
            break;
        }
        default:
            break;
        }


    }
}


void SP_Callback::ifAddNogoodCut(GRBEnv& env, Instance& ins, std::map<std::tuple<int, int>, double>& SP1_sol, int& sp2_approach, std::map<double,
    std::vector<SP1_Node_Sol*>>*feasibleSolList, double& sp_obj)
{
    // --------------------solve SP_2
    SP2_Sol sp2_sol = Solve_SP2(*Env, SP1_sol, ins.d_T, ins.num_l, ins.q0, sp2_approach);

    // std::cout << "Enter sp2 check" << std::endl;
    OptSPInformation->num_sp2_solved += 1;
    OptSPInformation->time_solve_sp2 += sp2_sol.sp2_solve_time;

    if (sp2_sol.sp2_status == 1)
    {
        // get a feasible solution
        std::cout << "sp2 is feasible, find a feasible solution to the original problem. " << std::endl;
        OptSPInformation->num_feasible_sol += 1;
        SP1_Node_Sol* feasibleNode = new SP1_Node_Sol();


        for (int ll = 0; ll < sp2_sol.A.size(); ll++)
        {
            feasibleNode->A.push_back(sp2_sol.A[ll]);
        }


        for (int pp = 1; pp <= ins.d_T; pp++)
            for (int ii = 1; ii <= ins.d_T; ii++)
            {
                double ev1 = 0;
                for (int jj = 1; jj <= ins.num_j; jj++)
                    ev1 += getSolution(vars[std::make_tuple(ii, pp, jj)]);
                if (ev1 > 0.5) feasibleNode->E.push_back(ii);
            }

        for (int ii = 1; ii <= ins.d_T; ii++)
            for (int jj = 1; jj <= ins.num_j; jj++)
            {
                double ev1 = 0;
                for (int pp = 1; pp <= ins.d_T; pp++)
                    ev1 += getSolution(vars[std::make_tuple(ii, pp, jj)]);
                if (ev1 > 0.5) feasibleNode->O.push_back(jj);
            }
        feasibleNode->obj_val = sp_obj;
        (*feasibleSolList)[sp_obj].push_back(feasibleNode);

    }
    else
    {
        // add no-good cut
        (*OptSPInformation).num_nogood_cut += sp2_sol.num_nogood_cut;
        OptSPInformation->num_callback_cut += sp2_sol.num_nogood_cut;
        std::cout << "Add no good cut: " << (*OptSPInformation).num_nogood_cut << std::endl;
        for (int no = 0; no < sp2_sol.sp2_nogood_set.size(); no++)
        {
            GRBLinExpr nogood_con = 0;
            for (int nn = 0; nn < sp2_sol.sp2_nogood_set[no].size(); nn++)
            {
                for (int j = 1; j <= ins.num_j; j++)
                {
                    nogood_con += vars[std::make_tuple(std::get<0>(sp2_sol.sp2_nogood_set[no][nn]), std::get<1>(sp2_sol.sp2_nogood_set[no][nn]), j)];
                }
            }

            addLazy(nogood_con, GRB_LESS_EQUAL, sp2_sol.sp2_nogood_set[no].size() - 1);
        }
    } // add no-good cut

}




// ------------------------------ Add Lazy cut 1 ----------------------------------------------------
void SP_Callback::FindAndAdd_Cut1(int& num_j, std::vector<int>& C1_set, std::vector<std::tuple<int, int>>& car_pos,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt)
{
    // add lazy cut1
    GRBLinExpr lazy1 = 0;
    for (int ci = 0; ci < C1_set.size(); ci++)
        for (int j = 1; j <= num_j; j++)
            lazy1 += var[std::make_tuple(C1_set[ci], std::get<1>(car_pos[C1_set[ci] - 1]), j)];

    addLazy(lazy1, GRB_LESS_EQUAL, C1_set.size() - 1);
}

// ------------------------------ Find cut1 - K=2
void SP_Callback::K2(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt)
{
    for (size_t k = 0; k < PosK2.size(); k++)
    {
        int ik = PosK2[k];
        int pk = std::get<1>(car_pos[ik - 1]);

        bool add_cut1 = false;
        for (size_t i1 = 1; i1 < ik; i1++)
        {
            int p1 = std::get<1>(car_pos[i1 - 1]);
            if (p1 < pk) continue;

            std::vector<int> C1_set;

            for (size_t ii = 1; ii < ik; ii++)
            {
                if (ii == i1) continue;
                int pp = std::get<1>(car_pos[ii - 1]);

                if ((ii < i1 && pp < p1) || (ii > i1 && pp > p1)) continue;     // ii is not in conflict with i1
                if (pp < pk) continue;                                      // ii is not in conflict with ik

                C1_set.push_back(ii);
            }

            if (C1_set.size() > (ins.num_l - 2) * ins.q0)
            {
                // add lazy cut1
                GRBLinExpr lazy1 = 0;
                add_cut1 = true;

                for (int j = 1; j <= ins.num_j; j++)
                {
                    lazy1 += vars[std::make_tuple(ik, pk, j)];
                    lazy1 += vars[std::make_tuple(i1, p1, j)];

                    for (int ci = 0; ci < C1_set.size(); ci++)
                        lazy1 += vars[std::make_tuple(C1_set[ci], std::get<1>(car_pos[C1_set[ci] - 1]), j)];
                }
                addLazy(lazy1, GRB_LESS_EQUAL, C1_set.size() + 1);
                Opt->num_lazy1 += 1;
                Opt->num_callback_cut += 1;
                std::cout << "Add lazy1 - k2: " << Opt->num_lazy1 << std::endl;
            }

            if (add_cut1) break;
        }
    }
}

// ------------------------------ Denest Version: Find cut1 - K=2 -------------------------
void K2_FindConflictPairs(int& i1, int& p1, int& ik, int& pk, std::vector<std::tuple<int, int>>& car_pos,
    std::vector<int>& C1_set)
{
    for (int ii = 1; ii < ik; ii++)
    {
        if (ii == i1) continue;
        int pp = std::get<1>(car_pos[ii - 1]);

        if ((ii < i1 && pp < p1) || (ii > i1 && pp > p1)) continue;     // ii is not in conflict with i1
        if (pp < pk) continue;                                          // ii is not in conflict with ik

        C1_set.push_back(ii);
    }
}


void SP_Callback::K2_Find_i1(int& ik, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt)
{
    int pk = std::get<1>(car_pos[ik - 1]);
    bool add_cut1 = false;

    for (int i1 = 1; i1 < ik; i1++)
    {
        int p1 = std::get<1>(car_pos[i1 - 1]);
        if (p1 < pk) continue;

        std::vector<int> C1_set;
        K2_FindConflictPairs(i1, p1, ik, pk, car_pos, C1_set);

        if (C1_set.size() > (ins.num_l - 2) * ins.q0)
        {
            add_cut1 = true;
            C1_set.push_back(ik);
            C1_set.push_back(i1);
            FindAndAdd_Cut1(ins.num_j, C1_set, car_pos, var, Opt);

            Opt->num_lazy1 += 1;
            Opt->num_callback_cut += 1;
            std::cout << "Add lazy1 - k2: " << Opt->num_lazy1 << std::endl;
        }

        if (add_cut1) break;
    }
}


void SP_Callback::newK2(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt)
{
    for (int k = 0; k < PosK2.size(); k++)
    {
        int ik = PosK2[k];
        K2_Find_i1(ik, ins, car_pos, var, Opt);
    }
}


// ------------------------------ Find cut1 - K=3
void SP_Callback::K3(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt)
{
    for (size_t k = 0; k < PosK2.size(); k++)
    {
        int ik = PosK2[k];
        int pk = std::get<1>(car_pos[ik - 1]);

        bool add_cut1 = false;
        for (size_t i2 = 2; i2 < ik; i2++)
        {
            int p2 = std::get<1>(car_pos[i2 - 1]);
            if (p2 < pk) continue;

            for (size_t i1 = 1; i1 < i2; i1++)
            {
                int p1 = std::get<1>(car_pos[i1 - 1]);
                if (p1 < p2) continue;

                std::vector<int> C1_set;
                for (int ii = 1; ii < ik; ii++)
                {
                    if (ii == i1 || ii == i2) continue;
                    int pp = std::get<1>(car_pos[ii - 1]);

                    if ((ii < i1 && pp < p1) || (ii > i1 && pp > p1)) continue;
                    if ((ii < i2 && pp < p2) || (ii > i2 && pp > p2)) continue;
                    if (pp < pk) continue;

                    C1_set.push_back(ii);
                }

                if (C1_set.size() > (ins.num_l - 3) * ins.q0)
                {
                    // add lazy cut1
                    GRBLinExpr lazy1 = 0;
                    add_cut1 = true;

                    for (int j = 1; j <= ins.num_j; j++)
                    {
                        lazy1 += vars[std::make_tuple(ik, pk, j)];
                        lazy1 += vars[std::make_tuple(i1, p1, j)];
                        lazy1 += vars[std::make_tuple(i2, p2, j)];

                        for (int ci = 0; ci < C1_set.size(); ci++)
                            lazy1 += vars[std::make_tuple(C1_set[ci], std::get<1>(car_pos[C1_set[ci] - 1]), j)];
                    }
                    addLazy(lazy1, GRB_LESS_EQUAL, C1_set.size() + 2);
                    Opt->num_lazy1 += 1;
                    Opt->num_callback_cut += 1;
                    std::cout << "Add lazy1 - k3: " << Opt->num_lazy1 << std::endl;
                }


                if (add_cut1) break;
            }

            if (add_cut1) break;


        }
    }

}

// ------------------------------ Denest Version: Find cut1 - K=3 -------------------------
void K3_FindConflictPairs(int& i1, int& p1, int& i2, int& p2, int& ik, int& pk, std::vector<std::tuple<int, int>>& car_pos,
    std::vector<int>& C1_set)
{
    for (int ii = 1; ii < ik; ii++)
    {
        if (ii == i1 || ii == i2) continue;
        int pp = std::get<1>(car_pos[ii - 1]);

        if ((ii < i1 && pp < p1) || (ii > i1 && pp > p1)) continue;
        if ((ii < i2 && pp < p2) || (ii > i2 && pp > p2)) continue;
        if (pp < pk) continue;

        C1_set.push_back(ii);
    }
}


void SP_Callback::K3_Find_i1(int& i2, int& p2, int& ik, int& pk, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt, bool& add_cut1)
{
    for (int i1 = 1; i1 < i2; i1++)
    {
        int p1 = std::get<1>(car_pos[i1 - 1]);
        if (p1 < p2) continue;

        std::vector<int> C1_set;
        K3_FindConflictPairs(i1, p1, i2, p2, ik, pk, car_pos, C1_set);

        if (C1_set.size() > (ins.num_l - 3) * ins.q0)
        {
            add_cut1 = true;
            C1_set.push_back(i1);
            C1_set.push_back(i2);
            C1_set.push_back(ik);
            FindAndAdd_Cut1(ins.num_j, C1_set, car_pos, var, Opt);

            Opt->num_lazy1 += 1;
            Opt->num_callback_cut += 1;
            std::cout << "Add lazy1 - k3: " << Opt->num_lazy1 << std::endl;
        }

        if (add_cut1) return;
    }
}


void SP_Callback::K3_Find_i2(int& ik, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt, bool& add_cut1)
{
    int pk = std::get<1>(car_pos[ik - 1]);
    for (int i2 = 2; i2 < ik; i2++)
    {
        int p2 = std::get<1>(car_pos[i2 - 1]);
        if (p2 < pk) continue;

        K3_Find_i1(i2, p2, ik, pk, ins, car_pos, var, Opt, add_cut1);
        if (add_cut1) return;
    }
}


void SP_Callback::newK3(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt)
{
    for (int k = 0; k < PosK2.size(); k++)
    {
        int ik = PosK2[k];
        bool add_cut1 = false;
        K3_Find_i2(ik, ins, car_pos, var, Opt, add_cut1);
    }
}



// ------------------------------ Find cut1 - K=4
void SP_Callback::K4(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt)
{
    for (size_t k = 0; k < PosK2.size(); k++)
    {
        int ik = PosK2[k];
        int pk = std::get<1>(car_pos[ik - 1]);

        bool add_cut1 = false;
        for (size_t i3 = 3; i3 < ik; i3++)
        {
            int p3 = std::get<1>(car_pos[i3 - 1]);
            if (p3 < pk) continue;

            for (size_t i2 = 2; i2 < i3; i2++)
            {
                int p2 = std::get<1>(car_pos[i2 - 1]);
                if (p2 < p3) continue;

                for (size_t i1 = 1; i1 < i2; i1++)
                {
                    int p1 = std::get<1>(car_pos[i1 - 1]);
                    if (p1 < p2) continue;

                    std::vector<int> C1_set;
                    for (int ii = 1; ii < ik; ii++)
                    {
                        if (ii == i1 || ii == i2 || ii == i3) continue;
                        int pp = std::get<1>(car_pos[ii - 1]);

                        if ((ii < i1 && pp < p1) || (ii > i1 && pp > p1)) continue;
                        if ((ii < i2 && pp < p2) || (ii > i2 && pp > p2)) continue;
                        if ((ii < i3 && pp < p3) || (ii > i3 && pp > p3)) continue;
                        if (pp < pk) continue;

                        C1_set.push_back(ii);
                    }

                    if (C1_set.size() > (ins.num_l - 4) * ins.q0)
                    {
                        // add lazy cut1
                        GRBLinExpr lazy1 = 0;
                        add_cut1 = true;

                        for (int j = 1; j <= ins.num_j; j++)
                        {
                            lazy1 += vars[std::make_tuple(ik, pk, j)];
                            lazy1 += vars[std::make_tuple(i1, p1, j)];
                            lazy1 += vars[std::make_tuple(i2, p2, j)];
                            lazy1 += vars[std::make_tuple(i3, p3, j)];

                            for (int ci = 0; ci < C1_set.size(); ci++)
                                lazy1 += vars[std::make_tuple(C1_set[ci], std::get<1>(car_pos[C1_set[ci] - 1]), j)];
                        }
                        addLazy(lazy1, GRB_LESS_EQUAL, C1_set.size() + 3);
                        Opt->num_lazy1 += 1;
                        Opt->num_callback_cut += 1;
                        std::cout << "Add lazy1 - k4: " << Opt->num_lazy1 << std::endl;
                    }
                    if (add_cut1) break;
                }
                if (add_cut1) break;
            }
            if (add_cut1) break;
        }
    }
}



// ------------------------------ Denest Version: Find cut1 - K=4 -------------------------
void K4_FindConflictPairs(int& i1, int& p1, int& i2, int& p2, int& i3, int& p3, int& ik, int& pk, std::vector<std::tuple<int, int>>& car_pos,
    std::vector<int>& C1_set)
{
    for (int ii = 1; ii < ik; ii++)
    {
        if (ii == i1 || ii == i2 || ii == i3) continue;
        int pp = std::get<1>(car_pos[ii - 1]);

        if ((ii < i1 && pp < p1) || (ii > i1 && pp > p1)) continue;
        if ((ii < i2 && pp < p2) || (ii > i2 && pp > p2)) continue;
        if ((ii < i3 && pp < p3) || (ii > i3 && pp > p3)) continue;
        if (pp < pk) continue;

        C1_set.push_back(ii);
    }
}


void SP_Callback::K4_Find_i1(int& i2, int& p2, int& i3, int& p3, int& ik, int& pk, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt, bool& add_cut1)
{
    for (int i1 = 1; i1 < i2; i1++)
    {
        int p1 = std::get<1>(car_pos[i1 - 1]);
        if (p1 < p2) continue;

        std::vector<int> C1_set;
        K4_FindConflictPairs(i1, p1, i2, p2, i3, p3, ik, pk, car_pos, C1_set);

        if (C1_set.size() > (ins.num_l - 4) * ins.q0)
        {
            add_cut1 = true;
            C1_set.push_back(i1);
            C1_set.push_back(i2);
            C1_set.push_back(i3);
            C1_set.push_back(ik);
            FindAndAdd_Cut1(ins.num_j, C1_set, car_pos, var, Opt);

            Opt->num_lazy1 += 1;
            Opt->num_callback_cut += 1;
            std::cout << "Add lazy1 - k4: " << Opt->num_lazy1 << std::endl;
        }

        if (add_cut1) return;
    }
}

void SP_Callback::K4_Find_i2(int& i3, int& p3, int& ik, int& pk, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt, bool& add_cut1)
{
    for (int i2 = 2; i2 < i3; i2++)
    {
        int p2 = std::get<1>(car_pos[i2 - 1]);
        if (p2 < p3) continue;

        K4_Find_i1(i2, p2, i3, p3, ik, pk, ins, car_pos, var, Opt, add_cut1);
        if (add_cut1) return;
    }

}

void SP_Callback::K4_Find_i3(int& ik, Instance& ins, std::vector<std::tuple<int, int>>& car_pos,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt, bool& add_cut1)
{
    int pk = std::get<1>(car_pos[ik - 1]);
    for (int i3 = 3; i3 < ik; i3++)
    {
        int p3 = std::get<1>(car_pos[i3 - 1]);
        if (p3 < pk) continue;

        K4_Find_i2(i3, p3, ik, pk, ins, car_pos, var, Opt, add_cut1);
        if (add_cut1) return;
    }
}


void SP_Callback::newK4(Instance& ins, std::vector<std::tuple<int, int>>& car_pos, std::vector<int>& PosK2,
    std::map<std::tuple<int, int, int>, GRBVar>& var, InformationSolveSP* Opt)
{
    for (int k = 0; k < PosK2.size(); k++)
    {
        int ik = PosK2[k];
        bool add_cut1 = false;
        K4_Find_i3(ik, ins, car_pos, var, Opt, add_cut1);
    }
}


// define function to solve SP by callback
void Solve_SP(InformationSolveSP& spSol_Info, std::map<double, std::vector<SP1_Node_Sol*>>& feasibleSolList,
    Instance& ins, std::map< std::tuple<int, int>, std::vector<int> >& Vbm, GRBEnv& env, bool& addC1, bool& addC2, bool& addC3, bool& LiftC2, bool& LiftC3,
    int sp2_method, double sp_time_limit)
{
    auto sp_t1 = std::chrono::steady_clock::now();
    // create the model
    GRBModel SP_model = GRBModel(env);

    // add variables
    std::vector< std::vector<GRBVar> > U(ins.d_T);     //U[p][w]: for linearization
    for (int p = 0; p < ins.d_T; p++)
    {
        std::vector<GRBVar> row(ins.num_w);
        for (int w = 0; w < ins.num_w; w++)
        {
            std::string var_name = "u[" + std::to_string(p + 1) + "][" + std::to_string(w + 1) + "]";
            GRBVar u = SP_model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, var_name);
            row[w] = u;
        }
        U[p] = row;
    }


    std::map<std::tuple<int, int, int>, GRBVar> E;
    for (int i = 1; i <= ins.d_T; i++)
    {
        for (int p = 1; p <= ins.d_T; p++)
        {
            for (int j = 1; j <= ins.num_j; j++)
            {
                std::string var_name = "e[" + std::to_string(i) + "][" + std::to_string(p) + "][" + std::to_string(j) + "]";
                GRBVar e = SP_model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
                E[std::make_tuple(i, p, j)] = e;
            }
        }
    }



    for (int i = 1; i <= ins.d_T; i++)
    {
        GRBLinExpr sp_con1 = 0;
        for (int p = 1; p <= ins.d_T; p++)
            for (int j = 1; j <= ins.num_j; j++)
                sp_con1 += E[std::make_tuple(i, p, j)];
        SP_model.addConstr(sp_con1, GRB_EQUAL, 1);

    }



    for (int p = 1; p <= ins.d_T; p++)
    {
        GRBLinExpr sp_con2 = 0;
        for (int i = 1; i <= ins.d_T; i++)
            for (int j = 1; j <= ins.num_j; j++)
                sp_con2 += E[make_tuple(i, p, j)];
        SP_model.addConstr(sp_con2, GRB_EQUAL, 1);
    }


    for (int p = 1; p <= ins.d_T; p++)
        for (int w = 1; w <= ins.num_w; w++)
        {
            GRBLinExpr sp_con3 = 0;

            for (int j = 1; j <= ins.num_j; j++)
            {
                for (int p_ = 1; p_ <= p; p_++)
                    for (int i = 1; i <= ins.d_T; i++)
                        sp_con3 += ins.r_j_w[j][w] * E[make_tuple(i, p_, j)];
            }
            SP_model.addConstr(U[p - 1][w - 1] - sp_con3 + p * ins.sigma_w[w], GRB_GREATER_EQUAL, 0.0);
            SP_model.addConstr(U[p - 1][w - 1] + sp_con3 - p * ins.sigma_w[w], GRB_GREATER_EQUAL, 0.0);
        }


    // demand constraints
    for (int i = 0; i < ins.d_bmj.size(); i++)
    {
        GRBLinExpr sp_con4 = 0;
        for (int v = 0; v < Vbm[std::make_tuple(ins.d_bmj[i].body, ins.d_bmj[i].color)].size(); v++)
            for (int p = 0; p < ins.d_T; p++)
                sp_con4 += E[std::make_tuple(Vbm[std::make_tuple(ins.d_bmj[i].body, ins.d_bmj[i].color)][v], p + 1, ins.d_bmj[i].config)];

        SP_model.addConstr(sp_con4, GRB_EQUAL, ins.d_bmj[i].demand);
    }

    // symmetry breaking
    //for (auto it = Vbm.begin(); it != Vbm.end(); it++)
    //{
    //    for (int i1 = 0; i1 < it->second.size() - 1; i1++)
    //        for (int i2 = i1 + 1; i2 < it->second.size(); i2++)
    //        {
    //            GRBLinExpr sp_con6 = 0;
    //            int car_1 = it->second[i1];
    //            int car_2 = it->second[i2];

    //            for (int p = 1; p <= ins.d_T; p++)
    //                for (int j = 1; j <= ins.J.size(); j++)
    //                    sp_con6 += p * (E[std::make_tuple(car_1, p, j)] - E[std::make_tuple(car_2, p, j)]);
    //            SP_model.addConstr(sp_con6, GRB_LESS_EQUAL, 1.0);
    //        }
    //}


    // Preprocess
    if (ins.d_T > ins.num_l * ins.q0)
    {
        for (int i = 1 + (ins.num_l - 1) * ins.q0 + 1; i < 1 + (ins.num_l - 1) * ins.q0 + ins.q0; i++)
        {
            int p_i = i - (1 + (ins.num_l - 1) * ins.q0);
            for (int p = 1; p <= p_i; p++)
            {
                GRBLinExpr sp_con5 = 0;
                for (int j = 1; j <= ins.num_j; j++)
                    sp_con5 += E[std::make_tuple(i, p, j)];
                SP_model.addConstr(sp_con5, GRB_EQUAL, 0.0);
            }
        }


        for (int i = 1 + ins.q0 * ins.num_l; i <= ins.d_T; i++)
        {
            int ti = (i - ins.num_l * ins.q0) / ins.num_l;
            int ri = (i - ins.num_l * ins.q0) % ins.num_l;
            int pi = ins.q0 + ti * ins.num_l + ri;

            for (int p = 1; p < pi; p++)
            {
                GRBLinExpr sp_con6 = 0;
                for (int j = 1; j <= ins.num_j; j++)
                {
                    sp_con6 += E[std::make_tuple(i, p, j)];
                }
                SP_model.addConstr(sp_con6, GRB_EQUAL, 0.0);
            }
        }
    }

    // set objective function
    GRBLinExpr sp_obj = 0;
    for (int p = 0; p < ins.d_T; p++)
        for (int w = 0; w < ins.num_w; w++)
            sp_obj += U[p][w];
    SP_model.setObjective(ins.gamma * sp_obj, GRB_MINIMIZE);


    // *******************************************************
    // optimze by Callback

    //SP_model.set(GRB_IntParam_PreCrush, 1);
    //SP_model.set(GRB_IntParam_Threads, 8);
    //SP_model.set(GRB_IntParam_Presolve, 0);
    //SP_model.set(GRB_DoubleParam_Heuristics, 0.0);
    //SP_model.set(GRB_IntParam_PreDual, 0);
    //SP_model.set(GRB_DoubleParam_WorkLimit, sp_time_limit);

    SP_model.set(GRB_IntParam_LazyConstraints, 1);
    SP_model.set(GRB_DoubleParam_TimeLimit, sp_time_limit);
    // SP_model.set(GRB_StringParam_SolFiles, "solution/mysol");
    int lif3_method = 0; // 0: add all cliques; 1: add until all pos occupied; 2: add no more than k cons
    int k_lif3 = 30;

    SP_Callback cb = SP_Callback(&env, ins, E, feasibleSolList, addC1, addC2, addC3, spSol_Info, LiftC2, LiftC3, lif3_method, k_lif3, sp2_method);
    SP_model.setCallback(&cb);
    //SP_model.write("model.msp");
    SP_model.optimize();
    // SP_model.write("model.mps");
    //SP_model.set(GRB_StringParam_SolFiles, "solution/mysol");

    std::cout << "Num.Cut1: " << spSol_Info.num_lazy1 << std::endl;
    std::cout << "Num.Cut2: " << spSol_Info.num_lazy2 << std::endl;
    std::cout << "Num.Cut3: " << spSol_Info.num_lazy3 << std::endl;
    std::cout << "Num.NoGoodCut: " << spSol_Info.num_nogood_cut << std::endl;
    std::cout << "Num.Lazy cuts: " << spSol_Info.num_lazy1 + spSol_Info.num_lazy2 + spSol_Info.num_lazy3 + spSol_Info.num_nogood_cut << std::endl;
    std::cout << "Num.feasible.solutions: " << feasibleSolList.size() << std::endl;
    // std::cout << "Num.feasible.solutions: " << spSol_Info.num_integer_sol << std::endl;
    std::cout << "---Num.feasible.solution.by GRB_IntAttr_SolCount: " << SP_model.get(GRB_IntAttr_SolCount) << std::endl;


    //for (int i = 1; i <= ins.d_T; i++)
    //    for (int p = 1; p <= ins.d_T; p++)
    //        for (int j = 1; j <= ins.num_j; j++)
    //            if (E[std::make_tuple(i, p, j)].get(GRB_DoubleAttr_X) > 0.5)
    //                std::cout << "(" << i << "," << p << "," << j << ")-";
    //std::cout << std::endl;



    if (SP_model.get(GRB_IntAttr_Status) != 2)
    {
        std::cout << "Do not solve optimally." << std::endl;
        auto sp_t3 = std::chrono::steady_clock::now();
        auto duration3 = std::chrono::duration_cast<std::chrono::milliseconds>(sp_t3 - sp_t1);
        spSol_Info.SP_total_time = SP_model.get(GRB_DoubleAttr_Runtime) * 1000;

        spSol_Info.SP_gap = -1;
        spSol_Info.SP_status = -1;
        spSol_Info.SP_OptObjVal = -1;
        return;
    }


    if (feasibleSolList.size() != 0)
    {
        // SP is solved optimally
        spSol_Info.SP_total_time = SP_model.get(GRB_DoubleAttr_Runtime) * 1000;
        spSol_Info.SP_gap = 0.0;
        spSol_Info.SP_status = 1;

        spSol_Info.SP_OptObjVal = feasibleSolList.begin()->first;

        for (int i = 0; i < ins.d_T; i++)
        {
            spSol_Info.E.push_back(feasibleSolList.begin()->second[0]->E[i]);
            spSol_Info.O.push_back(feasibleSolList.begin()->second[0]->O[i]);
            spSol_Info.A.push_back(feasibleSolList.begin()->second[0]->A[i]);
        }

        std::cout << "Solve SP optimally." << std::endl;
        std::cout << "SP Obj: " << spSol_Info.SP_OptObjVal << std::endl;
        return;
    }
    else
    {
        double loop_time = 0.0;
        SP2_Sol sp2_sol;


        while (feasibleSolList.size() == 0)
        {
            if (loop_time > sp_time_limit * 1000)
            {
                if (sp2_sol.sp2_status != 1)
                {
                    spSol_Info.SP_gap = -1;
                    spSol_Info.SP_OptObjVal = -1;
                    spSol_Info.SP_status = -1;
                }
                break;
            }

            auto tt1 = std::chrono::steady_clock::now();
            // do not check the feasiblity of sp2
            std::map<std::tuple<int, int>, double> s1_solution;
            for (int i = 1; i <= ins.d_T; i++)
                for (int p = 1; p <= ins.d_T; p++)
                {
                    double ev = 0.0;
                    for (int j = 1; j <= ins.num_j; j++)
                        ev += E[std::make_tuple(i, p, j)].get(GRB_DoubleAttr_X);
                    s1_solution[std::make_tuple(i, p)] = ev;
                }

            sp2_sol = Solve_SP2(env, s1_solution, ins.d_T, ins.num_l, ins.q0, sp2_method);
            spSol_Info.num_sp2_solved += 1;
            spSol_Info.time_solve_sp2 += sp2_sol.sp2_solve_time;

            if (sp2_sol.sp2_status == 1)
            {
                SP1_Node_Sol* feasibleNode = new SP1_Node_Sol();
                for (int ll = 0; ll < sp2_sol.A.size(); ll++)
                    feasibleNode->A.push_back(sp2_sol.A[ll]);

                for (int pp = 1; pp <= ins.d_T; pp++)
                    for (int ii = 1; ii <= ins.d_T; ii++)
                        if (s1_solution[std::make_tuple(ii, pp)] > 0.5)
                            feasibleNode->E.push_back(ii);


                for (int ii = 1; ii <= ins.d_T; ii++)
                    for (int jj = 1; jj <= ins.num_j; jj++)
                    {
                        double ev1 = 0;
                        for (int pp = 1; pp <= ins.d_T; pp++)
                            ev1 += E[std::make_tuple(ii, pp, jj)].get(GRB_DoubleAttr_X);
                        if (ev1 > 0.5) feasibleNode->O.push_back(jj);
                    }

                double sp_obj = SP_model.get(GRB_DoubleAttr_ObjVal);
                feasibleNode->obj_val = SP_model.get(GRB_DoubleAttr_ObjVal);
                feasibleSolList[sp_obj].push_back(feasibleNode);
                break;
            }
            else
            {
                // add no-good cut and reoptimize SP again
                spSol_Info.num_nogood_cut += sp2_sol.num_nogood_cut;
                for (int no = 0; no < sp2_sol.sp2_nogood_set.size(); no++)
                {
                    GRBLinExpr nogood_con = 0;
                    for (int nn = 0; nn < sp2_sol.sp2_nogood_set[no].size(); nn++)
                        for (int j = 1; j <= ins.num_j; j++)
                            nogood_con += E[std::make_tuple(std::get<0>(sp2_sol.sp2_nogood_set[no][nn]), std::get<1>(sp2_sol.sp2_nogood_set[no][nn]), j)];
                    SP_model.addConstr(nogood_con, GRB_LESS_EQUAL, sp2_sol.sp2_nogood_set[no].size() - 1);
                }


                // solve SP again
                std::cout << "------------------------------" << std::endl;
                std::cout << "!!!!!!!!! SP reoptimize again" << std::endl;
                std::cout << "------------------------------" << std::endl;

                spSol_Info.num_SP_solved += 1;

                SP_model.optimize();

            }

            auto tt2 = std::chrono::steady_clock::now();
            auto tt_dur = std::chrono::duration_cast<std::chrono::milliseconds>(tt2 - tt1);
            loop_time += tt_dur.count();

        }


        auto sp_t2 = std::chrono::steady_clock::now();
        auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(sp_t2 - sp_t1);
        spSol_Info.SP_total_time = duration2.count();


        std::cout << "+++++sp2 status: " << sp2_sol.sp2_status << std::endl;
        if (SP_model.get(GRB_IntAttr_Status) == 2)
        {
            spSol_Info.SP_gap = 0.0;
            spSol_Info.SP_status = 1;


            spSol_Info.SP_OptObjVal = feasibleSolList.begin()->first;

            for (int i = 0; i < ins.d_T; i++)
            {
                spSol_Info.E.push_back(feasibleSolList.begin()->second[0]->E[i]);
                spSol_Info.O.push_back(feasibleSolList.begin()->second[0]->O[i]);
                spSol_Info.A.push_back(feasibleSolList.begin()->second[0]->A[i]);
            }
            std::cout << "Solve SP optimally." << std::endl;
            std::cout << "SP Obj: " << spSol_Info.SP_OptObjVal << std::endl;

        }
        else
        {
            std::cout << "Do not solve SP optimally." << std::endl;
            spSol_Info.SP_gap = -1;
            spSol_Info.SP_status = -1;
            spSol_Info.SP_OptObjVal = -1;
        }

    }





}