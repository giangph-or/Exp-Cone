#include<iostream>
#include "ProblemReader.h"
#include "ResultReader.h"
#include "Conic.h"
#include "BranchAndCut.h"
#include "MILP.h"
#include "CuttingPlane.h"
#include "CPMarketExpansion.h"
#include "BCMarketExpansion.h"
#include "StochasticCuttingPlane.h"
#include "AEMarketExpansion.h"
#include "QCAEMarketExpansion.h"
#include "ApproximateExp.h"
#include "Quadratic.h"
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
    if (model == "BC") {
        string out_file = "result_bc//" + instance_name + "_" + no_pay + ".txt";
        BC bcoa(data, time_limit, out_file);
        bcoa.solve(data);
        //bcoa.solve_approximation(data);
        //bcoa.solve_multicut(data, 5);
    }
    if (model == "MILP") {
        string out_file = "result_milp//" + instance_name + "_" + no_pay + ".txt";
        MILP milpmcoa(data, time_limit, out_file);
        milpmcoa.solve(data);
    }
    if (model == "Conic") {
        string out_file = "result_conic//" + instance_name + "_" + no_pay + ".txt";
        Conic conicmcoa(data, time_limit, out_file);
        conicmcoa.solve(data);
    }
    if (model == "AE") {
        string out_file = "result_ae//" + instance_name + "_" + no_pay + ".txt";
        ApproximateExp aeoa(data, time_limit, out_file);
        aeoa.solve(data);
    }
    if (model == "Q") {
        string out_file = "result_q//" + instance_name + "_" + no_pay + ".txt";
        Quadratic qoa(data, time_limit, out_file);
        qoa.solve(data);
    }

    //Market Expansion model
    if (model == "CPME") {
        string out_file = "result_cpme//" + instance_name + "_" + no_pay + ".txt";
        CPMarketExpansion cpmeoa(data, time_limit, out_file);
        vector<double> beta = cpmeoa.create_beta(data);
        vector<double> gamma = cpmeoa.create_gamma(data);
        vector<double> kappa = cpmeoa.calculate_kappa(data, beta);
        cpmeoa.solve_market_expansion(data, kappa, beta, gamma);
    }
    if (model == "BCME") {
        string out_file = "result_bcme//" + instance_name + "_" + no_pay + ".txt";
        BCMarketExpansion bcmeoa(data, time_limit, out_file);
        vector<double> beta = bcmeoa.create_beta(data);
        vector<double> gamma = bcmeoa.create_gamma(data);
        vector<double> kappa = bcmeoa.calculate_kappa(data, beta);
        bcmeoa.solve_market_expansion(data, kappa, beta, gamma);
    }
    if (model == "AEME") {
        //data.print_data();
        string out_file = "result_aeme//" + instance_name + "_" + no_pay + ".txt";
        AEMarketExpansion aemeoa(data, time_limit, out_file);
        vector<double> beta = aemeoa.create_beta(data);
        vector<double> gamma = aemeoa.create_gamma(data);
        vector<double> kappa = aemeoa.calculate_kappa(data, beta);
        aemeoa.solve_market_expansion(data, kappa, beta, gamma);
    }
    if (model == "QCAEME") {
        //data.print_data();
        string out_file = "result_qcaeme//" + instance_name + "_" + no_pay + ".txt";
        AEMarketExpansion qcaemeoa(data, time_limit, out_file);
        vector<double> beta = qcaemeoa.create_beta(data);
        vector<double> gamma = qcaemeoa.create_gamma(data);
        vector<double> kappa = qcaemeoa.calculate_kappa(data, beta);
        qcaemeoa.solve_market_expansion(data, kappa, beta, gamma);
    }
    if (model == "QCCPME") {
        //data.print_data();
        string out_file = "result_qccpme//" + instance_name + "_" + no_pay + ".txt";
        AEMarketExpansion qccpmeoa(data, time_limit, out_file);
        vector<double> beta = qccpmeoa.create_beta(data);
        vector<double> gamma = qccpmeoa.create_gamma(data);
        vector<double> kappa = qccpmeoa.calculate_kappa(data, beta);
        qccpmeoa.solve_market_expansion(data, kappa, beta, gamma);
    }


    //run command line: dataFile 0 Read
    //===============================================================================
    //==============================READ RESULTS=====================================
    //===============================================================================
    //Basic models and Market Expansion models on new dataset
    if (model == "Read") {
        Report report;
        vector<string> nopay;
        vector<string> nGroup;
        vector<string> capacity;
        if (instance_name == "50_5_overlap" || instance_name == "50_10_overlap" 
            || instance_name == "50_15_overlap" || instance_name == "50_20_overlap") {
            nopay = {"1", "2"};
            nGroup = {"2", "5", "10"};
            capacity = {"5", "10", "20"};
        }
        if (instance_name == "100_10_overlap" || instance_name == "100_20_overlap" 
            || instance_name == "100_15_overlap" || instance_name == "100_5_overlap") {
            nopay = {"1", "2"};
            nGroup = {"5", "10", "20"};
            capacity = {"5", "10", "20"};
        }
        if (instance_name == "200_20_overlap" || instance_name == "200_10_overlap" 
            || instance_name == "200_15_overlap" || instance_name == "200_5_overlap") {
             nopay = {"5", "10"};
             nGroup = {"5", "10", "20"};
             capacity = {"5", "10", "20"};
        }
        if (instance_name == "500_50_overlap") {
            nopay = {"10", "20"};
            nGroup = {"5", "10", "20"};
            capacity = {"5", "10", "20"};
        }
        if (instance_name == "500_100_overlap") {
            nopay = { "10", "20"};
            nGroup = { "5", "10", "20"};
            capacity = { "5", "10", "20"};
        }
        if (instance_name == "200_100_overlap") {
            nopay = { "5", "10" };
            nGroup = { "5", "10", "20" };
            capacity = { "5", "10", "20" };
        }
        report.create_report(instance_name, nGroup, capacity, nopay);
    }

    if (model == "ReadME") {
        Report report;
        vector<string> nopay;
        vector<string> nGroup;
        vector<string> capacity;
        if (instance_name == "50_5_overlap" || instance_name == "50_10_overlap"
            || instance_name == "50_15_overlap" || instance_name == "50_20_overlap") {
            nopay = { "2", "4"};
            nGroup = { "2", "5", "10" };
            capacity = { "5", "10", "20" };
        }
        if (instance_name == "100_10_overlap" || instance_name == "100_20_overlap"
            || instance_name == "100_15_overlap" || instance_name == "100_5_overlap") {
            nopay = { "2", "4" };
            nGroup = { "5", "10", "20" };
            capacity = { "5", "10", "20" };
        }
        if (instance_name == "200_20_overlap" || instance_name == "200_10_overlap"
            || instance_name == "200_15_overlap" || instance_name == "200_5_overlap") {
            nopay = { "5", "10" };
            nGroup = { "5", "10", "20" };
            capacity = { "5", "10", "20" };
        }
        if (instance_name == "500_50_overlap") {
            nopay = { "10", "20" };
            nGroup = { "5", "10", "20" };
            capacity = { "5", "10", "20" };
        }
        if (instance_name == "500_100_overlap") {
            nopay = { "10", "20" };
            nGroup = { "5", "10", "20" };
            capacity = { "5", "10", "20" };
        }
        if (instance_name == "200_100_overlap") {
            nopay = { "5", "10" };
            nGroup = { "5", "10", "20" };
            capacity = { "5", "10", "20" };
        }
        report.create_report(instance_name, nGroup, capacity, nopay);
    }

    //CP and Conic models on Sen et al (2018) dataset
    if (model == "SenCP") {
        string budget = argv[4];
        instance_file = "Sen_data//" + instance_name + ".dat";
        Data data_Sen;
        data_Sen.read_data_Sen(instance_file, noPay, stod(budget));
        string out_file = "result_cp//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        CuttingPlane cpoa(data_Sen, time_limit, out_file);
        cpoa.solve(data_Sen);
    }

    if (model == "SenConic") {
        string budget = argv[4];
        //string subBudget = argv[5];
        instance_file = "Sen_data//" + instance_name + ".dat";
        Data data_Sen;
        data_Sen.read_data_Sen(instance_file, noPay, stod(budget));
        //data_Sen.read_data_Sen_general(instance_file, noPay, stod(budget), stod(subBudget));
        //data_Sen.print_data();
        //string out_file = "result_conic//" + instance_name + "_" + no_pay + "_" + budget + "_" + subBudget + ".txt";
        string out_file = "result_conic//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        Conic conicmcoa(data_Sen, time_limit, out_file);
        conicmcoa.solve(data_Sen);

        //string budget = argv[4];
        //instance_file = "Sen_data//" + instance_name + ".dat";
        //Data data_Sen;
        //data_Sen.read_data_Sen(instance_file, noPay, stod(budget));

        //vector<double> alpha(data_Sen.number_customers, -1);
        //for (int i = 0; i < data_Sen.number_customers; ++i)
        //    for (int j = 0; j < data_Sen.number_products; ++j)
        //        if (data_Sen.revenue[i][j] > alpha[i])
        //            alpha[i] = data_Sen.revenue[i][j];

        //double total_alpha = 0;
        //for (int i = 0; i < data_Sen.number_customers; ++i)
        //    total_alpha += data_Sen.fraction[i] * alpha[i];

        //instance_file = "lpmodel//" + instance_name + "_" + no_pay + "_cap" + budget + ".lp";
        //GRBEnv env = GRBEnv(true);
        //env.start();
        //GRBModel model = GRBModel(env, instance_file);
        //model.set(GRB_DoubleParam_TimeLimit, 600);
        //model.set(GRB_IntParam_Threads, 1);
        //auto start = chrono::steady_clock::now();        
        //model.optimize();
        //auto end = chrono::steady_clock::now();
        //std::chrono::duration<double> solvingTime = end - start;
        //string out_file = "result_conic//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        //ofstream report_results(out_file, ofstream::out);
        //report_results.precision(10);
        //report_results << model.get(GRB_DoubleAttr_ObjVal) << " " << total_alpha - model.get(GRB_DoubleAttr_ObjVal) << " " << solvingTime.count();
        //report_results.close();
    }

    if (model == "SenMILP") {
        //string budget = argv[4];
        //instance_file = "Sen_data//" + instance_name + ".dat";
        //Data data_Sen;
        //data_Sen.read_data_Sen(instance_file, noPay, stod(budget));
        ////data_Sen.print_data();
        //string out_file = "result_conic//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        //Conic conicmcoa(data_Sen, time_limit, out_file);
        //conicmcoa.solve(data_Sen);

        string budget = argv[4];
        instance_file = "Sen_data//" + instance_name + ".dat";
        Data data_Sen;
        data_Sen.read_data_Sen(instance_file, noPay, stod(budget));

        vector<double> alpha(data_Sen.number_customers, -1);
        for (int i = 0; i < data_Sen.number_customers; ++i)
            for (int j = 0; j < data_Sen.number_products; ++j)
                if (data_Sen.revenue[i][j] > alpha[i])
                    alpha[i] = data_Sen.revenue[i][j];

        double total_alpha = 0;
        for (int i = 0; i < data_Sen.number_customers; ++i)
            total_alpha += alpha[i];

        instance_file = "milpmodel//" + instance_name + "_" + no_pay + "_cap" + budget + ".lp";
        GRBEnv env = GRBEnv(true);
        env.start();
        GRBModel model = GRBModel(env, instance_file);
        model.set(GRB_DoubleParam_TimeLimit, 600);
        auto start = chrono::steady_clock::now();
        model.optimize();
        auto end = chrono::steady_clock::now();
        std::chrono::duration<double> solvingTime = end - start;
        string out_file = "result_milp//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        ofstream report_results(out_file, ofstream::out);
        report_results.precision(10);
        report_results << model.get(GRB_DoubleAttr_ObjVal) << " " << total_alpha - model.get(GRB_DoubleAttr_ObjVal) << " " << solvingTime.count();
        report_results.close();
    }

    if (model == "SenQ") {
        string budget = argv[4];
        instance_file = "Sen_data//" + instance_name + ".dat";
        Data data_Sen;
        data_Sen.read_data_Sen(instance_file, noPay, stod(budget));
        //data_Sen.print_data();
        string out_file = "result_q//" + instance_name + "_" + no_pay + "_" + budget + ".txt";
        Quadratic qoa(data_Sen, time_limit, out_file);
        qoa.solve(data_Sen);
    }

    if (model == "SenRead") {
        Report report;
        vector<string> nopay;
        vector<string> capacity;
        if (instance_name == "100_100_hard") {
            nopay = { "1", "2" };
            capacity = { "10", "20", "50", "100" };
            //nopay = { "1", "2" };
            //capacity = { "5_2", "10_4", "25_10", "50_20" };
            //capacity = { "5_2", "10_4", "25_10"};
        }
        if (instance_name == "200_20") {
            nopay = { "5", "10" };
            capacity = { "10", "20", "50", "100", "200" };
        }
        if (instance_name == "500_50") {
            nopay = { "10", "20" };
            capacity = { "20", "50", "100", "200", "500" };
        }
        if (instance_name == "generalized_capacity") {
            nopay = { "10", "20" };
            //capacity = { "5_2", "10_4", "25_10", "50_20" };
            capacity = { "5_2", "10_4", "25_10", "50_20", "100_40" };
        }
        if (instance_name == "500_50_generalized") {
            nopay = { "10", "20" };
            //capacity = { "5_2", "10_4", "25_10", "50_20", "125_50" };
            capacity = { "10_4", "25_10", "50_20", "125_50", "250_100"};
        }
        if (instance_name == "1000_100_generalized") {
            nopay = { "10", "20" };
            capacity = { "25", "50", "100", "250", "500" };
            //capacity = { "25_10", "50_20", "125_50", "250_100", "500_200"};
            //capacity = { "10_4", "25_10", "50_20", "125_50", "250_100"};
            //capacity = { "5_2", "10_4", "25_10", "50_20", "125_50" };
        }
        if (instance_name == "100_1000_generalized") {
            nopay = { "1", "2"};
            capacity = { "10", "20", "50"};
        }
        if (instance_name == "100_5000_generalized") {
            nopay = { "1", "2" };
            capacity = { "10", "20", "50" };
        }
        if (instance_name == "200_2000_generalized") {
            nopay = { "5", "10"};
            capacity = { "10", "20", "50"};
        }
        if (instance_name == "200_4000_generalized") {
            nopay = { "5", "10" };
            capacity = { "10", "20", "50" };
        }
        report.create_report_Sen(instance_name, nopay, capacity);
    }
}