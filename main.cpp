#include<iostream>
#include "ProblemReader.h"
#include "ResultReader.h"
#include "Conic.h"
#include "CuttingPlane.h"
#include <string>
#include <iomanip>
#include <chrono>
using namespace std;

int main(int argc, char* argv[]) {
    string instance_name = argv[1];
    string no_pay = argv[2];
    string model = argv[3];
    string instance_file = "AO_data//" + instance_name + ".dat";
    double time_limit = 3600;
    double noPay = stod(no_pay);
    Data data;
    data.read_data(instance_file, noPay);

    //run command line: dataFile noPay modelName
    //===============================================================================
    //==============================SOLVE MODELS=====================================
    //===============================================================================
    //Basic model
    if (model == "CP") {
        string out_file = "result_cp//" + instance_name + "_" + no_pay + ".txt";
        CuttingPlane cpoa(data, time_limit, out_file);
        cpoa.solve(data);
        //cpoa.solve_multicut(data, 5);
    }
    if (model == "Conic") {
        string out_file = "result_conic//" + instance_name + "_" + no_pay + ".txt";
        Conic conicmcoa(data, time_limit, out_file);
        conicmcoa.solve(data);
}
