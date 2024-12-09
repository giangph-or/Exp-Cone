#include <iostream>
#include <math.h>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <gurobi_c++.h>
#include "ProblemReader.h"
using namespace std;

class Conic {
public:
    Data data;
    string output;
    string log_file;
    string out_res_csv;
    double time_limit;
    double obj_val;
    double gap;
    double time_for_solve;
    char var_name[1000];
    Conic();
    Conic(Data data, double time_limit, string outfile);
    double calculate_master_obj(Data data, vector<int> x);
    double calculate_optimal_bound_denominator(Data data, int i);
    int calculate_bound_y(Data data);
    double optimal_bound_y_in(Data data, int i, int j, double alpha);
    double optimal_bound_y_notin(Data data, int i, int j, double alpha);
    void solve(Data data);
};
