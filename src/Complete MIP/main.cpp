#include "FullModel.h"
#include <iostream>

using namespace std;

int main(int arg, char** argv)
{
    std::string ins_name = argv[1];
    std::string ins_f = "../BS/" + ins_name + ".txt";
    const char* ins_f_name = ins_f.c_str();

    // read instance
    Instance ins;
    readIns(ins_f_name, ins);



    std::string txt_file_name = "../Re-Complete/";
    txt_file_name += ins_name;
    txt_file_name += ".txt";


    std::string csv_file_name1 = "./CompleteResults.csv";
    std::ofstream out_csv_f1;
    out_csv_f1.open(csv_file_name1, ios::app);

    //double timeLim = std::stof(argv[2]);

    out_csv_f1 << "Ins" << ',' << "GRB Status" << ',' << "Solve optimally" << ',' << "Gap" << ',' << "Total.Time" << ',' << "Opt.obj" << ','
        << "Obj.F1" << ',' << "Obj.F2" << ',' << "Obj.F3" << std::endl;


    SolutionInfo* solution = new SolutionInfo();
    solution->num_car = ins.d_T;
    CompleteModel(ins, ins_name, csv_file_name1, txt_file_name, 3600.0, solution);

}