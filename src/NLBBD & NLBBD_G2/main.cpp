#include<iostream>
#include "LBBDG2.h"
#include "ReadIns.h"
#include "NLBBD.h"

using namespace std;

int main(int arg, char** argv)
{
    std::string ins_name = argv[1];
    std::string ins_f = "../BS/" + ins_name + ".txt";
    const char* ins_f_name = ins_f.c_str();
    //const char* ins_f_name = "C://CppProjects//DeFG1//x64//BS//10_1.txt";


    // read instance
    Instance ins;
    readIns(ins_f_name, ins);



    std::string txt_file_name = "../Re1/";      // results for NLBBD
    // std::string txt_file_name = "../Re/";    // results for LBBDg2
    txt_file_name += ins_name;
    txt_file_name += ".txt";


    std::string csv_file_name1 = "../NLBBD-Com.csv";  // results for NLBBD
    std::string csv_file_name2 = "../NLBBD-Iter.csv"; // results for NLBBD
    // std::string csv_file_name1 = "../Test-Com.csv";  // results for LBBDg2
    // std::string csv_file_name2 = "../Test-Iter.csv"; // results for LBBDg2
    std::ofstream out_csv_f1;
    std::ofstream out_csv_f2;
    out_csv_f1.open(csv_file_name1, ios::app);
    out_csv_f2.open(csv_file_name2, ios::app);

    // file for LBBDg2
    /*
    * out_csv_f1 << "Ins" << ',' << "LB.2nd.Stage" << ',' << "Time.LB" << ','
        << "Num.Nodes" << ',' << "Num.Arcs" << ',' << "Time.Build.G2" << ',' << "Get opt sol by SP-Heu" << ','
        << "Obj.MP" << ',' << "Obj.SP" << ',' << "Obj.F1F2" << ',' << "Obj.theta" << ','
        << "Total.Time" << ',' << "Num.LBBD.Iters" << ','
        << "UB" << ',' << "LB" << ','
        << "OP.Status" << ',' << "MP.Status" << ',' << "SP.Status" << ','
        << "Gap" << ',' << "MP.gap" << ',' << "SP.gap" << ',' << "Num.LBBD.OptCut" << ',' << "Num.LBBD.NogoodCut" << ','
        << "Total.Time.MP" << ','
        << "Total.Time.SP" << ',' << "Num.SP.Solved" << ','
        << "Avg.Time.MP" << ','
        << "Avg.Time.SP" << ','
        << "Num.Lazy1" << ',' << "Num.Lazy2" << ',' << "Num.Lazy3" << ',' << "Num.SP-Nogood.Cuts" << ','
        << "Avg.Num.Lazy1" << ',' << "Avg.Num.Lazy2" << ',' << "Avg.Num.Lazy3" << ',' << "Avg.Num.SP-No.Cuts" << ','
        << "Time.Solve.SP2-S2" << ',' << "Num.SP-S2.Solved" << ',' << "Avg.Time.Solve.SP-S2"
        << std::endl;
    */

    // file for NLBBD
    out_csv_f1 << "Ins" << ',' << "LB.2nd.Stage" << ',' << "Time.LB" << ','
        << "Get opt sol by SP-Heu" << ','
        << "Obj.MP" << ',' << "Obj.SP" << ',' << "Obj.F1F2" << ',' << "Obj.theta" << ','
        << "Total.Time" << ',' << "Num.LBBD.Iters" << ','
        << "UB" << ',' << "LB" << ','
        << "OP.Status" << ',' << "MP.Status" << ',' << "SP.Status" << ','
        << "Gap" << ',' << "MP.gap" << ',' << "SP.gap" << ',' << "Num.LBBD.OptCut" << ',' << "Num.LBBD.NogoodCut" << ','
        << "Total.Time.MP" << ','
        << "Total.Time.SP" << ',' << "Num.SP.Solved" << ','
        << "Avg.Time.MP" << ','
        << "Avg.Time.SP" << ','
        << "Num.Lazy1" << ',' << "Num.Lazy2" << ',' << "Num.Lazy3" << ',' << "Num.SP-Nogood.Cuts" << ','
        << "Avg.Num.Lazy1" << ',' << "Avg.Num.Lazy2" << ',' << "Avg.Num.Lazy3" << ',' << "Avg.Num.SP-No.Cuts" << ','
        << "Time.Solve.SP2-S2" << ',' << "Num.SP-S2.Solved" << ',' << "Avg.Time.Solve.SP-S2"
        << std::endl;


    // f1: solve by Heu - mpObj - spObj - total time - num.LBBD.Iter - UB - LB
    // OP Status - MP Status - SP Status - Gap - mp.gap - sp.gap - num.opt.cut - num.nogood.cut -
    // total.MP.time - total.SP.time - avg.time.MP - ave.time.SP - num.lazy1 - num.lazy2 - num.lazy3 - num.sp.nogood
    // avg.num.lazy1 - ave.num.lazy2 - avg.num.lazy3 - avg.num.sp.nogood - total.time.sp2 - total.num.sp2 - avg.time.sp2 



    out_csv_f2 << "Ins" << ',' << "Iter" << ',' << "UB" << ',' << "LB" << ',' << "Time.MP" << ','
        << "MP.Status" << ',' << "Obj.MP" << ',' << "MP.theta" << ',' << "MP.F1F2" << ',' << "MP.gap" << ','
        << "Time.SP" << ','
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


    bool addCut1 = true;
    bool addCut2 = true;
    bool addCut3 = true;
    bool LiftCut2 = false;
    bool LiftCut3 = false;
    int sp2_method = 1;
    bool use_Heu = 1;
    double time_limit = 3600.0;


    out_csv_f1.close();
    out_csv_f2.close();

    NLBBD(ins, ins_name, txt_file_name, csv_file_name1, csv_file_name2, addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, use_Heu, time_limit);
    //LBBDg2(ins, ins_name, txt_file_name, csv_file_name1, csv_file_name2, addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, use_Heu, time_limit);



    //// use shell script to solve
    //std::string ins_name = argv[1];
    //std::string ins_f = "../BS/" + ins_name + ".txt";
    //const char* ins_f_name = ins_f.c_str();


    //// read instance
    //Instance ins;
    //int insal = std::stoi(argv[2]);
    //int insbe = std::stoi(argv[3]);
    //int insga = std::stoi(argv[4]);
    //int insq = std::stoi(argv[5]);
    //int insl = std::stoi(argv[6]);
    //readIns(ins_f_name, ins, insq, insl, insal, insbe, insga);


    //std::string a1 = "1";
    //std::string a2 = "1";
    //std::string a3 = "1";
    //std::string li2 = "1";
    //std::string li3 = "1";
    //std::string uheu = "1";


    //bool addCut1 = (a1 == argv[7]);
    //bool addCut2 = (a2 == argv[8]);
    //bool addCut3 = (a3 == argv[9]);
    //bool LiftCut2 = (li2 == argv[10]);
    //bool LiftCut3 = (li3 == argv[11]);
    //int sp2_method = std::stoi(argv[12]);


    //bool useHeuristic = (uheu == argv[13]);
    //double timeLimit = std::stof(argv[14]);



    //std::string ad1 = argv[7];
    //std::string ad2 = argv[8];
    //std::string ad3 = argv[9];
    //std::string lift2 = argv[10];
    //std::string lift3 = argv[11];
    //std::string sp2m = argv[12];
    //std::string uheus = argv[13];


    //std::string txt_file_name = "../Re-G2/" + ad1 + ad2 + ad3 + lift2 + lift3 + sp2m + uheus + "/" + ins_name + ".txt";
    //std::string csv_file_name1 = "./" + ad1 + ad2 + ad3 + lift2 + lift3 + sp2m + uheus + "-Com.csv";
    //std::string csv_file_name2 = "./" + ad1 + ad2 + ad3 + lift2 + lift3 + sp2m + uheus + "-Iter.csv";
    //std::ofstream out_csv_f1;
    //std::ofstream out_csv_f2;
    //out_csv_f1.open(csv_file_name1, ios::app);
    //out_csv_f2.open(csv_file_name2, ios::app);


    //out_csv_f1 << "Ins" << ',' << "LB.2nd.Stage" << ',' << "Time.LB" << ','
    //    << "Num.Nodes" << ',' << "Num.Arcs" << ',' << "Time.Build.G2" << ',' << "Get opt sol by SP-Heu" << ','
    //    << "Obj.MP" << ',' << "Obj.SP" << ',' << "Obj.F1F2" << ',' << "Obj.theta" << ','
    //    << "Total.Time" << ',' << "Num.LBBD.Iters" << ','
    //    << "UB" << ',' << "LB" << ','
    //    << "OP.Status" << ',' << "MP.Status" << ',' << "SP.Status" << ','
    //    << "Gap" << ',' << "MP.gap" << ',' << "SP.gap" << ',' << "Num.LBBD.OptCut" << ',' << "Num.LBBD.NogoodCut" << ','
    //    << "Total.Time.MP" << ','
    //    << "Total.Time.SP" << ',' << "Num.SP.Solved" << ','
    //    << "Avg.Time.MP" << ','
    //    << "Avg.Time.SP" << ','
    //    << "Num.Lazy1" << ',' << "Num.Lazy2" << ',' << "Num.Lazy3" << ',' << "Num.SP-Nogood.Cuts" << ','
    //    << "Avg.Num.Lazy1" << ',' << "Avg.Num.Lazy2" << ',' << "Avg.Num.Lazy3" << ',' << "Avg.Num.SP-No.Cuts" << ','
    //    << "Time.Solve.SP2-S2" << ',' << "Num.SP-S2.Solved" << ',' << "Avg.Time.Solve.SP-S2"
    //    << std::endl;




    //out_csv_f2 << "Ins" << ',' << "Iter" << ',' << "UB" << ',' << "LB" << ',' << "Time.MP" << ','
    //    << "MP.Status" << ',' << "Obj.MP" << ',' << "MP.theta" << ',' << "MP.F1F2" << ',' << "MP.gap" << ','
    //    << "Time.SP" << ','
    //    << "SP.Status" << ',' << "Obj.SP" << ',' << "SP.Gap" << ','
    //    << "Num.Lazy1" << ',' << "Num.Lazy2" << ','
    //    << "Num.Lazy3" << ',' << "Num.SP-Nogood" << ',' << "Num.SP-S2.solved" << ',' << "Time.Solve.SP-S2" << ','
    //    << "Avg.Time.SolveSP2" << ','
    //    << "Num.IntegerSol.Found" << ',' << "Num.FeasibleSol.Found" << ',' << "Num.lazy3.added" << ','
    //    << "Num.reduced.lazy3.added"
    //    << std::endl;


    //out_csv_f1.close();
    //out_csv_f2.close();


    // LBBDg2(ins, ins_name, txt_file_name, csv_file_name1, csv_file_name2, addCut1, addCut2, addCut3, LiftCut2, LiftCut3, sp2_method, useHeuristic, timeLimit);





}