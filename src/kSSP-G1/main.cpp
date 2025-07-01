// #include "solution_newg1.h"
// #include "solution2.h"  // new code to generate the new G1 2024/7/13
#include "LBBDg1.h"

int main(int arg, char** argv)
{
    std::string ins_name = argv[1];
    std::string ins_f = "../BS/" + ins_name + ".txt";
    const char* ins_f_name = ins_f.c_str();
    //const char* ins_f_name = "C://CppProjects//DeFG1//x64//BS//10_1.txt";


    // read instance
    Instance ins;
    readIns(ins_f_name, ins);


    // CAR_COUNT = ins.cars.size();
    // solution_newG1 solution_new_g1;
    // Solve(ins.cars, ins.c0, ins.alpha, ins.beta, solution_new_g1);
    // std::cout << "Solution of the new g1 size: " << std::to_string(solution_new_g1.solution.size()) << std::endl;


    std::string txt_file_name = "../Re/";
    // std::string txt_file_name = "../RTest1/";
    txt_file_name += ins_name;
    txt_file_name += ".txt";


    std::string g1_f = "../GBS/" + ins_name + ".txt";
    const char* g1_path = g1_f.c_str();

    std::string csv_file_name1 = "../NLBBD1-Com.csv";
    std::string csv_file_name2 = "../NLBBD1-Iter.csv";
    std::ofstream out_csv_f1;
    std::ofstream out_csv_f2;
    out_csv_f1.open(csv_file_name1, ios::app);
    out_csv_f2.open(csv_file_name2, ios::app);


    out_csv_f1 << "Ins" << ',' << "LB.2nd.Stage" << ',' << "Time.LB" << ','
        << "Num.Nodes" << ',' << "Num.Arcs" << ',' << "Time.Build.G1" << ',' << "Time.Find.SP" << ',' << "ShortestPath.Cost" << ','
        << "Solve optimally by Heu" << ',' << "Kahn_Heu is feasible" << ',' << "Kahn_SP2 is feasible" << ','
        << "Obj.MP" << ',' << "Obj.SP" << ',' << "Total.Time" << ',' << "Num.LBBD.Iters" << ',' << "UB" << ',' << "currentLB" << ','
        << "OP.Status" << ',' << "MP.Status" << ',' << "SP.Status" << ',' << "Gap" << ',' << "MP.gap" << ',' << "SP.gap" << ','
        << "Total.Time.MP" << ','
        << "Total.Time.SP" << ',' << "Num.SP-1.Solved" << ','
        << "Avg.Time.MP" << ','
        << "Avg.Time.SP" << ','
        << "Num.Lazy1" << ',' << "Num.Lazy2" << ',' << "Num.Lazy3" << ',' << "Num.SP-Nogood.Cuts" << ','
        << "Avg.Num.Lazy1" << ',' << "Avg.Num.Lazy2" << ',' << "Avg.Num.Lazy3" << ',' << "Avg.Num.SP-No.Cuts" << ','
        << "Time.Solve.SP2-S2" << ',' << "Num.SP-S2.Solved" << ',' << "Avg.Time.Solve.SP-S2"
        << std::endl;


    // f1: solve by Heu - Kahn_Heu - Kahn+SP2 - mpObj - spObj - total time - num.LBBD.Iter - UB - LB
           // OP Status - MP Status - SP Status - Gap - mp.gap - sp.gap - 
           // total.MP.time - total.SP.time - avg.time.MP - ave.time.SP - num.lazy1 - num.lazy2 - num.lazy3 - num.sp.nogood
           // avg.num.lazy1 - ave.num.lazy2 - avg.num.lazy3 - avg.num.sp.nogood - total.time.sp2 - total.num.sp2 - avg.time.sp2 
           // num.lazy3 added - num.reduced.lazy3.added


    out_csv_f2 << "Ins" << ',' << "Iter" << ',' << "UB" << ',' << "currentLB" << ',' << "Time.MP" << ','
        << "MP.Status" << ',' << "Obj.MP" << ',' << "MP.Gap" << ',' << "Time.SP" << ','
        << "SP.Status" << ',' << "Obj.SP" << ',' << "SP.Gap" << ','
        << "Num.Lazy1" << ',' << "Num.Lazy2" << ','
        << "Num.Lazy3" << ',' << "Num.SP-Nogood" << ',' << "Num.SP-S2.solved" << ',' << "Time.Solve.SP-S2" << ','
        << "Avg.Time.SolveSP2" << ','
        << "Num.IntegerSol.Found" << ',' << "Num.FeasibleSol.Found" << ',' << "Num.lazy3.added" << ','
        << "Num.reduced.lazy3.added"
        << std::endl;

    // f2: Iter th - UB - LB - MP.time - mp.status - mp.obj - mp.gap - SP.time - sp.status - sp.obj - sp.gap - num.lazy1 - num.lazy2 
    // - num.lazy3 - num.sp.nogood
    // num.sp2.solved - time.sp2.solved - avg.time.sp2 - num.integer.found - num.feasible.solution.found - 
    // num.lazy3.added - num.reduced.lazy3.added



    out_csv_f1.close();
    out_csv_f2.close();

    //LBBDg1(Instance ins, std::string const& ins_name, std::string const& txtfile,
    //    std::string const& csv1, std::string const& csv2, std::string const& KSP_method,
    //    bool& addCut1, bool& addCut2, bool& addCut3, bool& LiftCut2, bool& LiftCut3, int& sp2_method, bool& Kahn_Heu, bool other_Heu, double timeLimit)

    bool addCut1 = 1;
    bool addCut2 = 1;
    bool addCut3 = 1;
    bool LiftCut2 = 0;
    bool LiftCut3 = 0;
    int sp2_method = 1;     // 1: feasible problem; 2: maximization problem; 3: minimization problem
    bool useHeuristic = 0;
    double timeLimit = 3600.0;


    LBBDg1(ins, ins_name, txt_file_name, csv_file_name1, csv_file_name2, "PSB",
        addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, useHeuristic, useHeuristic, timeLimit);

    remove(g1_path);


}