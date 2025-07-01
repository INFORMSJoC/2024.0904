#pragma once
#include<vector>
#include<tuple>
#include "FindCliques.h"
#include "ReadIns.h"
#include "BuildG2.h"


// -------------- (E,O,A) the optimal solution
class optSol
{
public:
    double obj_opt;
    double obj_mp;
    double obj_sp;

    std::vector<int> bodyList;
    std::vector<int> colorList;
    std::vector<int> configList;
    std::vector<int> laneList;
    std::vector<int> downPosList;

    optSol() : obj_opt(-1), obj_mp(-1), obj_sp(-1) {}

};



// -------------- the information of subproblem
class InformationSolveSP
{
public:
    int SP_status;
    double SP_gap;
    double SP_total_time;

    int num_integer_sol;
    int num_feasible_sol;

    int num_lazy1;
    int num_lazy2;
    int num_lazy3;
    int num_nogood_cut;     // according to the feasibility of the SP-S2
    std::vector<int> numLazy3List;
    std::vector<int> reduced_numlazy3List;

    double SP_OptObjVal;
    std::vector<int> O;
    std::vector<int> E;
    std::vector<int> A;

    int num_sp2_solved;
    double time_solve_sp2;
    int num_SP_solved;

    int num_callback_cut;

    InformationSolveSP() : SP_status(1), SP_gap(0.0), SP_total_time(0.0), num_integer_sol(0), num_feasible_sol(0),
        num_lazy1(0), num_lazy2(0), num_lazy3(0), num_nogood_cut(0), SP_OptObjVal(0.0), num_sp2_solved(0), time_solve_sp2(0.0), num_SP_solved(1),
        num_callback_cut(0)
    {}

    void spInfoUpdate();

};

void InformationSolveSP::spInfoUpdate()
{
    SP_status = 1;
    SP_gap = 0.0;
    SP_total_time = 0.0;
    num_integer_sol = 0;
    num_feasible_sol = 0;
    num_lazy1 = 0;
    num_lazy2 = 0;
    num_lazy3 = 0;
    num_nogood_cut = 0;
    SP_OptObjVal = 0.0;
    num_sp2_solved = 0;
    time_solve_sp2 = 0.0;
    num_SP_solved = 1;
    num_callback_cut = 0;

    numLazy3List.clear();
    reduced_numlazy3List.clear();
    O.clear();
    E.clear();
    A.clear();

}



// -------------- the information of the final result
class LBBD_Info
{
public:
    Graph2* graph2; // get after building G1

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

    int total_SP_solved;
    double total_time_SP;

    int total_lazy1;
    int total_lazy2;
    int total_lazy3;

    int total_num_SP2_solved;
    double total_time_SP2;
    int total_num_SP2_NogoodCut;

    double total_time_MP;

    double Total_time;

    InformationSolveSP* lastSP;

    double LB_delta_obj;    // get after obtain delta
    double LB_delta_time;
    std::vector<int> LB_Config_Seq;

    bool Heu_sp;    // if the total problem is solved optimally by the heuristic

    // final optimal solution
    std::vector<int> lB;
    std::vector<int> lC;
    std::vector<int> lO;
    std::vector<int> lE;
    std::vector<int> lA;

    int lbbd_cut_opt;
    int lbbd_cut_nogood;

    double lbbd_gap;
    double mp_gap;
    double sp_gap;

    bool Kahn_Heu;   // if the Kahn is feasible by the heuristic-E + heuristic-A
    bool Kahn_SP2;   // if the Kahn is feasible by Kahn + the heuristic-E + sp2


    LBBD_Info() : LBBD_iter(0), MP_status(-1), SP_status(-1), OP_status(-1), UB(1000000), LB(1e-3), obj_MP(-1), obj_SP(-1),
        obj_mp_theta(0.0), obj_mp_f1f2(-1),
        total_SP_solved(0), total_time_SP(0.0),
        total_lazy1(0), total_lazy2(0), total_lazy3(0),
        total_num_SP2_solved(0), total_time_SP2(0.0), total_num_SP2_NogoodCut(0), total_time_MP(0.0), Total_time(0.0),
        LB_delta_obj(0.0), LB_delta_time(0.0), Heu_sp(0),
        lbbd_cut_opt(-1), lbbd_cut_nogood(0),
        lbbd_gap(-1), mp_gap(-1), sp_gap(-1), Kahn_Heu(0), Kahn_SP2(0),
        graph2(nullptr), lastSP(nullptr) {}


    ~LBBD_Info() { Clear(); }
    void Clear()
    {
        graph2->Clear();
    }


};




// -------------- judge if the solution E gotten by solving SP-Stage1 will violate lazy cut 1
bool is_E_violateCut1(std::vector<std::tuple<int, int>>& UpDownSeq, Instance& ins)
{
    bool Cut1_feasible = true;
    // just use k = 2 to judge
    for (int ik = 1 + ins.q0 * (ins.num_l - 1); ik <= ins.d_T; ik++)
    {
        int pk = std::get<1>(UpDownSeq[ik - 1]);
        if (ik - pk >= ins.q0 * (ins.num_l - 1))
        {
            for (int i1 = ik - 2 - ins.q0 * (ins.num_l - 2); i1 >= 1; i1--)
            {
                int p1 = std::get<1>(UpDownSeq[i1 - 1]);
                if (p1 - pk > (ins.num_l - 2) * ins.q0)
                {
                    int middle = 0;
                    for (int ii = i1 + 1; ii < ik; ii++)
                    {
                        int pp = std::get<1>(UpDownSeq[ii - 1]);
                        if (pp < p1 && pp > pk)
                            middle += 1;
                    }

                    if (middle > (ins.num_l - 2) * ins.q0)
                    {
                        Cut1_feasible = false;
                        return Cut1_feasible;
                    }
                }
            }
        }

    }
    return Cut1_feasible;
}


// -------------- judge if the solution E gotten by solving SP-Stage1 will violate lazy cut 2
bool is_E_violateCut2(std::vector<std::tuple<int, int>>& UpDownSeq, Instance& Ins)
{
    bool Cut2_feasible = true;

    for (int k1 = 1; k1 <= Ins.d_T; k1++)
    {
        int p1 = std::get<1>(UpDownSeq[k1 - 1]);
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


            for (int k = k1 + 1; k <= Ins.d_T; k++)
            {
                int p = std::get<1>(UpDownSeq[k - 1]);
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
                int p = std::get<1>(UpDownSeq[k - 1]);
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
                int p = std::get<1>(UpDownSeq[k - 1]);
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
                Cut2_feasible = false;
                return Cut2_feasible;
            }
        }   // |L| >= 2
        else
        {
            // |L| = 1
            if (C2_size > 0 || C1_size > 0)
            {
                Cut2_feasible = false;
                return Cut2_feasible;
            }
        }   // |L| = 1

    } // finish adding lazy cut 2

    return Cut2_feasible;
}



// -------------- judge if the solution E gotten by solving SP-Stage1 will violate lazy cut 3
bool is_E_violateCut3(std::vector<std::tuple<int, int>>& UpDownSeq, Instance& Ins)
{
    bool Cut3_feasible = true;

    std::vector<std::tuple<int, int>> adj;
    for (int i1 = 1; i1 < Ins.d_T; i1++)
        for (int i2 = i1 + 1; i2 <= Ins.d_T; i2++)
        {
            int p1 = std::get<1>(UpDownSeq[i1 - 1]);
            int p2 = std::get<1>(UpDownSeq[i2 - 1]);
            if (p1 > p2) adj.push_back(std::make_tuple(i1, i2));
        }
    BK cliques(Ins.d_T, adj);
    cliques.run();
    for (auto it = cliques.maxCliques.begin(); it != cliques.maxCliques.end(); it++)
        if ((*it).size() > Ins.num_l)
        {
            Cut3_feasible = false;
            return Cut3_feasible;
        }

    return Cut3_feasible;
}



// -------------- judge if the solution E gotten by solving SP-Stage1 will violate the three valid inequalities
bool is_E_feasible(std::vector<std::tuple<int, int>>& UpDownSeq, Instance ins)
{
    bool E_feasible = true;

    if (is_E_violateCut1(UpDownSeq, ins) == true)
    {
        if (is_E_violateCut2(UpDownSeq, ins) == true)
        {
            if (is_E_violateCut2(UpDownSeq, ins) == false)
            {
                E_feasible = false;
                return E_feasible;
            }
        }
        else
        {
            E_feasible = false;
            return E_feasible;
        }
    }
    else
    {
        E_feasible = false;
        return E_feasible;
    }

    return E_feasible;
}



// -------------- judge if two vector are same
bool is_equal(std::vector<int>& va, std::vector<int>& vb)
{
    bool eq = false;
    if (va.size() == vb.size())
    {
        auto ap = va.begin();
        auto bp = vb.begin();
        for (int ei = 0; ei < va.size(); ei++)
        {
            if ((*ap) != (*bp)) return eq;
            ap++;
            bp++;
        }
        return eq = true;
    }
    return eq;
}

























