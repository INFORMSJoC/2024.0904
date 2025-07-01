#include "ReadIns.h"
#include<iostream>
#include "BasicLBBD.h"


using namespace std;

int main(int arg, char** argv)
{
    std::string ins_name = argv[1];
    std::string ins_f = "../BS/" + ins_name + ".txt";
    const char* ins_f_name = ins_f.c_str();



    // read instance
    Instance ins;
    readIns(ins_f_name, ins);


    // ------------------------------- Baisc MP + Basic SP
    std::string txt_file_name = "../Re-BasicLBBD/";
    txt_file_name += ins_name;
    txt_file_name += ".txt";


    std::string csv_file_name1 = "../Basic-C.csv";
    std::string csv_file_name2 = "../Basic-I.csv";
    std::ofstream out_csv_f1;
    std::ofstream out_csv_f2;
    out_csv_f1.open(csv_file_name1, ios::app);
    out_csv_f2.open(csv_file_name2, ios::app);



    out_csv_f1 << "Ins" << ',' << "LB.2nd.Stage" << ',' << "Time.LB" << ',' << "OP.status" << ',' << "MP.Status" << ',' << "SP.Status" << ','
        << "Num.Iteration" << ',' << "Gap" << ',' << "Total.Time" << ',' << "Total.Time.MP" << ',' << "Total.Time.SP" << ','
        << "Avg.Time.MP" << ',' << "Avg.Time.SP" << ','
        << "UB" << ',' << "LB" << ',' << "Obj.MP" << ',' << "Obj-F1F2" << ',' << "Obj-theta" << ','
        << "Obj-SP" << ',' << "Num.OptCut" << ',' << "Num.NogoodCut" << std::endl;


    out_csv_f2 << "Ins" << ',' << "Iter" << ',' << "UB" << ',' << "LB" << ',' << "MP.Time" << ',' << "MP.Status" << ',' << "Obj.MP" << ','
        << "Obj.theta" << ',' << "Obj.f1f2" << ',' << "MP.Gap" << ',' << "SP.Time" << ',' << "SP.Status" << ',' << "Obj.SP" << ',' << "SP.Gap" << std::endl;


    BasicLBBD(ins, ins_name, txt_file_name, csv_file_name1, csv_file_name2, 80.0);





}