#include "Conic.h"
#include <chrono>
#include <algorithm>

Conic::Conic() {

}

Conic::Conic(Data data, double time_limit, string outfile) {
    this->data = data;
    this->time_limit = time_limit;
    this->out_res_csv = outfile;
}

CBSen::CBSen() {

}

CBSen::CBSen(GRBVar* varx, GRBVar* vary, int p, int c, vector<double> n, vector<vector<double>> u, vector<vector<double>> r) {
    x = varx;
    y = vary;
    products = p;
    customers = c;
    noPay = n;
    util = u;
    ren = r;
}

double Conic::calculate_master_obj(Data data, vector<int> x) {
    //double obj = 0;
    //for (int i = 0; i < data.number_customers; ++i) {
    //    double ts = 0, ms = data.no_purchase[i];
    //    for (int j = 0; j < data.number_products; ++j) {
    //        ts += data.revenue[i][j] * x[j] * data.utilities[i][j];
    //        ms += x[j] * data.utilities[i][j];
    //    }
    //    obj += data.fraction[i] * ts / ms;
    //}
    //return obj;

    double obj = 0;
    for (int i = 0; i < data.number_customers; ++i) {
        double ts = 0, ms = data.no_purchase[i];
        for (int j = 0; j < data.number_products; ++j) {
            ts += data.revenue[i][j] * x[j] * data.utilities[i][j];
            ms += x[j] * data.utilities[i][j];
        }
        obj += ts / ms;
    }
    return obj;
}

double Conic::calculate_optimal_bound_denominator(Data data, int i) {
    GRBEnv env = GRBEnv(true);
    env.start();

    GRBModel model = GRBModel(env);

    //Decison variables: x_j = 1 if product j is chosen, 0 otherwise
    GRBVar* x = 0;
    x = model.addVars(data.number_products, GRB_BINARY);

    //Set capacity constraints
    for (int s = 0; s < data.number_sets; ++s) {
        GRBLinExpr sum;
        for (int j = 0; j < data.number_products; ++j)
            if (data.in_set[j][s] == 1)
                sum += data.cost[j] * x[j];
        model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
    }

    //GRBLinExpr cost;
    //for (int j = 0; j < data.number_products; ++j)
    //    cost += data.fraction2[j] * x[j];
    //model.addConstr(cost <= data.capacity_each_set, "ct_set_cap");

    GRBLinExpr obj;
    for (int j = 0; j < data.number_products; ++j)
        obj += data.utilities[i][j] * x[j];
    model.setObjective(obj, GRB_MAXIMIZE);

    model.set(GRB_DoubleParam_TimeLimit, 600);
    model.set(GRB_IntParam_OutputFlag, 0);

    model.optimize();

    return model.get(GRB_DoubleAttr_ObjVal) + data.no_purchase[i];
}

int Conic::calculate_bound_y(Data data) {
    vector<int> chosen(data.number_products);
    vector<pair<double, int>> u(data.number_products);
    for (int j = 0; j < data.number_products; ++j)
        u[j] = make_pair(data.cost[j], j);

    sort(u.begin(), u.end());
    int number_chosen = 0;
    for (int s = 0; s < data.number_sets; ++s) {
        double current_capacity = 0;
        int count = 0;
        while (count < data.number_products) {
            if (data.in_set[u[count].second][s] == 1 && current_capacity + data.cost[u[count].second] <= data.capacity_each_set) {
                current_capacity += data.cost[u[count].second];
                number_chosen++;
            }
            count++;
        }
    }
    for (int j = 0; j < data.number_products; ++j) {
        int check_in_set = 0;
        for (int s = 0; s < data.number_sets; ++s)
            if (data.in_set[j][s] == 1) check_in_set = 1;
        if (check_in_set == 0) number_chosen++;
    }
    if (number_chosen >= data.number_products) number_chosen = data.number_products;
    return number_chosen;
}

vector<vector<double>> Conic::calculate_bound_y_in(Data data) {
    vector<vector<double>> lb_in(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        lb_in[i].resize(data.number_products);

    for (int i = 0; i < data.number_customers; ++i) {
        for (int k = 0; k < data.number_products; ++k) {
            vector<pair<double, int>> u(data.number_products);
            for (int j = 0; j < data.number_products; ++j)
                u[j] = make_pair(data.utilities[i][j], j);

            sort(u.begin(), u.end(), greater<pair<double, int>>());
            for (int j = 0; j < data.number_products; ++j)
                if (u[j].second == k) {
                    u.erase(u.begin() + j);
                    break;
                }

            lb_in[i][k] = data.utilities[i][k];

            int count = 0;
            //while (count < calculate_bound_y(data)) {
            while (count < data.capacity_each_set-1) {
                lb_in[i][k] += u[count].first;
                count++;
            }

            ////General
            //double cost = 0;
            //for (int a = 0; a < u.size(); ++a)
            //    if (cost + data.fraction2[u[a].second] <= data.capacity_each_set) {
            //        lb_in[i][k] += u[a].first;
            //        cost += data.fraction2[u[a].second];
            //    }
        }
    }
    return lb_in;
}

vector<vector<double>> Conic::calculate_bound_y_notin(Data data) {
    vector<vector<double>> lb_notin(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        lb_notin[i].resize(data.number_products);

    for (int i = 0; i < data.number_customers; ++i) {
        for (int k = 0; k < data.number_products; ++k) {
            vector<pair<double, int>> u(data.number_products);
            for (int j = 0; j < data.number_products; ++j)
                u[j] = make_pair(data.utilities[i][j], j);

            sort(u.begin(), u.end(), greater<pair<double, int>>());
            for (int j = 0; j < data.number_products; ++j)
                if (u[j].second == k) {
                    u.erase(u.begin() + j);
                    break;
                }

            lb_notin[i][k] = 0;

            int count = 0;
            //while (count < calculate_bound_y(data)) {
            while (count < data.capacity_each_set) {
                lb_notin[i][k] += u[count].first;
                count++;
            }

            ////General
            //double cost = 0;
            //for (int a = 0; a < u.size(); ++a)
            //    if (cost + data.fraction2[u[a].second] <= data.capacity_each_set) {
            //        lb_notin[i][k] += u[a].first;
            //        cost += data.fraction2[u[a].second];
            //    }
        }
    }
    return lb_notin;
}

//for overlap cases
double Conic::optimal_bound_y_in(Data data, int i, int j, double alpha) {
    GRBEnv env = GRBEnv(true);
    env.start();

    GRBModel model = GRBModel(env);

    //Decison variables: x_j = 1 if product j is chosen, 0 otherwise
    GRBVar* x = 0;
    x = model.addVars(data.number_products, GRB_BINARY);

    model.addConstr(x[j] == 1);

    //Set capacity constraints
    for (int s = 0; s < data.number_sets; ++s) {
        GRBLinExpr sum;
        for (int k = 0; k < data.number_products; ++k)
            if (data.in_set[k][s] == 1)
                sum += data.cost[k] * x[k];
        model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
    }

    GRBLinExpr obj = alpha * data.no_purchase[i];
    for (int k = 0; k < data.number_products; ++k)
        obj += (alpha - data.revenue[i][k]) * data.utilities[i][k] * x[k];
    model.setObjective(obj, GRB_MAXIMIZE);

    model.set(GRB_DoubleParam_TimeLimit, 600);
    model.set(GRB_IntParam_OutputFlag, 0);

    model.optimize();

    return model.get(GRB_DoubleAttr_ObjVal);
}

//for overlap cases
double Conic::optimal_bound_y_notin(Data data, int i, int j, double alpha) {
    GRBEnv env = GRBEnv(true);
    env.start();

    GRBModel model = GRBModel(env);

    //Decison variables: x_j = 1 if product j is chosen, 0 otherwise
    GRBVar* x = 0;
    x = model.addVars(data.number_products, GRB_BINARY);

    model.addConstr(x[j] == 0);

    //Set capacity constraints
    for (int s = 0; s < data.number_sets; ++s) {
        GRBLinExpr sum;
        for (int k = 0; k < data.number_products; ++k)
            if (data.in_set[k][s] == 1)
                sum += data.cost[k] * x[k];
        model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
    }

    GRBLinExpr obj = alpha * data.no_purchase[i];
    for (int k = 0; k < data.number_products; ++k)
        obj += (alpha - data.revenue[i][k]) * data.utilities[i][k] * x[k];
    model.setObjective(obj, GRB_MAXIMIZE);

    model.set(GRB_DoubleParam_TimeLimit, 600);
    model.set(GRB_IntParam_OutputFlag, 0);

    model.optimize();

    return model.get(GRB_DoubleAttr_ObjVal);
}

vector<vector<double>> Conic::subset_bound_y_in(Data data) {
    vector<vector<double>> lb_in(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        lb_in[i].resize(data.number_products);

    for (int i = 0; i < data.number_customers; ++i) {
        for (int k = 0; k < data.number_products; ++k) {
            vector<pair<double, int>> u(data.number_products);
            for (int j = 0; j < data.number_products; ++j)
                u[j] = make_pair(data.utilities[i][j], j);

            sort(u.begin(), u.end(), greater<pair<double, int>>());
            for (int j = 0; j < data.number_products; ++j)
                if (u[j].second == k) {
                    u.erase(u.begin() + j);
                    break;
                }

            lb_in[i][k] = data.utilities[i][k];
            vector<int> in(10, data.sub_capacity_each_set);

            for (int a = 0; a < u.size(); ++a)
                for (int s = 0; s < 10; ++s)
                    if (in[s] > 0 && data.in_set[u[a].second][s] == 1) {
                        lb_in[i][k] += u[a].first;
                        in[s]--;
                    }
        }
    }
    return lb_in;
}

vector<vector<double>> Conic::subset_bound_y_notin(Data data) {
    vector<vector<double>> lb_notin(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        lb_notin[i].resize(data.number_products);

    for (int i = 0; i < data.number_customers; ++i) {
        for (int k = 0; k < data.number_products; ++k) {
            vector<pair<double, int>> u(data.number_products);
            for (int j = 0; j < data.number_products; ++j)
                u[j] = make_pair(data.utilities[i][j], j);

            sort(u.begin(), u.end(), greater<pair<double, int>>());
            for (int j = 0; j < data.number_products; ++j)
                if (u[j].second == k) {
                    u.erase(u.begin() + j);
                    break;
                }

            lb_notin[i][k] = 0;
            vector<int> in(10, data.sub_capacity_each_set);
            for (int a = 0; a < u.size(); ++a)
                for (int s = 0; s < 10; ++s)
                    if (in[s] > 0 && data.in_set[u[a].second][s] == 1) {
                        lb_notin[i][k] += u[a].first;
                        in[s]--;
                    }
        }
    }
    return lb_notin;
}

void Conic::solve(Data data) {
    auto start = chrono::steady_clock::now(); //get start time
    vector<double> alpha(data.number_customers, -1);
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            if (data.revenue[i][j] > alpha[i])
                alpha[i] = data.revenue[i][j];

    vector<double> upper_bound_denominator(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i)
        upper_bound_denominator[i] = calculate_optimal_bound_denominator(data, i);

    //vector<vector<double>> bound_in = calculate_bound_y_in(data);
    //vector<vector<double>> bound_notin = calculate_bound_y_notin(data);
    
    //for overlap cases
    vector<vector<double>> bound_in(data.number_customers);
    vector<vector<double>> bound_notin(data.number_customers);
    for (int i = 0; i < data.number_customers; ++i) {
        bound_in[i].resize(data.number_products, 0);
        bound_notin[i].resize(data.number_products, 0);
    }
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j) {
            bound_in[i][j] = optimal_bound_y_in(data, i, j, alpha[i]);
            bound_notin[i][j] = optimal_bound_y_notin(data, i, j, alpha[i]);
        }

    //General
    //vector<vector<double>> subset_bound_in = subset_bound_y_in(data);
    //vector<vector<double>> subset_bound_notin = subset_bound_y_notin(data);

    GRBEnv env = GRBEnv(true);
    env.start();

    GRBModel model = GRBModel(env);

    //cout << "Decison variables : x_j\n" << endl;
    GRBVar* x;
    x = new GRBVar[data.number_products];
    for (int j = 0; j < data.number_products; ++j)
        x[j] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(j));

    vector<double> upper_bound_y(data.number_customers, 0);
    vector<double> lower_bound_y(data.number_customers, 0);
    for (int i = 0; i < data.number_customers; ++i) {
        lower_bound_y[i] = 1 / upper_bound_denominator[i];
        upper_bound_y[i] = 1 / data.no_purchase[i];
    }

    //cout << "Slack variables : y_i\n" << endl;
    GRBVar* y;
    y = new GRBVar[data.number_customers];
    for (int i = 0; i < data.number_customers; ++i)
        y[i] = model.addVar(lower_bound_y[i], upper_bound_y[i], 0, GRB_CONTINUOUS, "y_" + to_string(i));
        //y[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "y_" + to_string(i));

    //cout << "Slack variables: z_{ij}\n" << endl;
    GRBVar** z;
    z = new GRBVar * [data.number_customers];
    for (int i = 0; i < data.number_customers; ++i)
        z[i] = new GRBVar[data.number_products];
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            //z[i][j] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "z_" + to_string(i) + "_" + to_string(j));
            z[i][j] = model.addVar(0, upper_bound_y[i], 0, GRB_CONTINUOUS, "z_" + to_string(i) + "_" + to_string(j));

    //cout << "Slack variables: w_i\n" << endl;
    GRBVar* w = 0;
    w = new GRBVar[data.number_customers];
    for (int i = 0; i < data.number_customers; ++i)
        w[i] = model.addVar(data.no_purchase[i], GRB_INFINITY, 0, GRB_CONTINUOUS, "w_" + to_string(i));

    //cout << "Constraints related to w and x\n" << endl;
    for (int i = 0; i < data.number_customers; ++i) {
        GRBLinExpr sum_x;
        for (int j = 0; j < data.number_products; ++j)
            sum_x += x[j] * data.utilities[i][j];
        model.addConstr(sum_x + data.no_purchase[i] == w[i]);
    }

    //cout << "Constraints related to x, w and z\n" << endl;
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j)
            model.addQConstr(z[i][j] * w[i] >= x[j] * x[j]);

    //cout << "Constraints related to y and w\n" << endl;
    for (int i = 0; i < data.number_customers; ++i)
        model.addQConstr(y[i] * w[i] >= 1);

    //cout << "Constraints related to y and z\n" << endl;
    for (int i = 0; i < data.number_customers; ++i) {
        GRBLinExpr sum_z;
        for (int j = 0; j < data.number_products; ++j)
            sum_z += data.utilities[i][j] * z[i][j];

        model.addConstr(data.no_purchase[i] * y[i] + sum_z >= 1);
    }

    //cout << "McCornick constraints" << endl;
    for (int i = 0; i < data.number_customers; ++i)
        for (int j = 0; j < data.number_products; ++j) {
            model.addConstr(z[i][j] <= x[j] * (1 / (data.no_purchase[i] + data.utilities[i][j])));

            //if (bound_in[i][j] >= subset_bound_in[i][j])
                model.addConstr(z[i][j] >= x[j] * (1 / (data.no_purchase[i] + bound_in[i][j])));
            //else
                //model.addConstr(z[i][j] >= x[j] * (1 / (data.no_purchase[i] + subset_bound_in[i][j])));

            //if (bound_notin[i][j] >= subset_bound_notin[i][j])
                model.addConstr(z[i][j] <= y[i] - (1 - x[j]) * (1 / (data.no_purchase[i] + bound_notin[i][j])));
            //else
                //model.addConstr(z[i][j] <= y[i] - (1 - x[j]) * (1 / (data.no_purchase[i] + subset_bound_notin[i][j])));

            model.addConstr(z[i][j] >= y[i] - (1 - x[j]) * (1 / data.no_purchase[i]));
        }

    //Set capacity constraints
    for (int s = 0; s < data.number_sets; ++s) {
        GRBLinExpr sum;
        for (int j = 0; j < data.number_products; ++j)
            if (data.in_set[j][s] == 1)
                sum += data.cost[j] * x[j];
        model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
    }

    ////General
    //for (int s = 0; s < 10; ++s) {
    //    GRBLinExpr sum;
    //    for (int j = 0; j < data.number_products; ++j)
    //        if (data.in_set[j][s] == 1)
    //            sum += data.cost[j] * x[j];
    //    model.addConstr(sum <= data.sub_capacity_each_set, "ct_set_cap" + to_string(s));
    //}

    //GRBLinExpr cost;
    //for (int j = 0; j < data.number_products; ++j)
    //    cost += data.fraction2[j] * x[j];
    //model.addConstr(cost <= data.capacity_each_set, "ct_set_cap");

    //cout << "Objective\n" << endl;
    //GRBLinExpr obj;
    //for (int i = 0; i < data.number_customers; ++i) {
    //    obj += data.fraction[i] * alpha[i] * data.no_purchase[i] * y[i];
    //    for (int j = 0; j < data.number_products; ++j)
    //        obj += data.fraction[i] * (alpha[i] - data.revenue[i][j]) * z[i][j] * data.utilities[i][j];
    //}

    //for overlap cases
    GRBLinExpr obj;
    for (int i = 0; i < data.number_customers; ++i) {
        obj += alpha[i] * data.no_purchase[i] * y[i];
        for (int j = 0; j < data.number_products; ++j)
            obj += (alpha[i] - data.revenue[i][j]) * z[i][j] * data.utilities[i][j];
    }

    model.setObjective(obj, GRB_MINIMIZE);

    auto time_now = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = time_now - start;

    model.set(GRB_DoubleParam_TimeLimit, time_limit - elapsed_seconds.count());
    model.set(GRB_IntParam_MIQCPMethod, 1);
    model.set(GRB_IntParam_Threads, 8);
    //model.set(GRB_IntParam_LazyConstraints, 1);
    //model.set(GRB_IntParam_PreCrush, 1);
    //model.set(GRB_DoubleParam_MIPGap, 1e-8);
    model.write("conic.lp");
    //model.set(GRB_IntParam_OutputFlag, 0);

    //CBSen cb = CBSen(x, y, data.number_products, data.number_customers, data.no_purchase, data.utilities, data.revenue);
    //model.setCallback(&cb);

    model.optimize();

    vector<int> x_sol(data.number_products);

    double master_obj = 0;

    if (model.get(GRB_IntAttr_SolCount) != 0) {
        cout << "\nResult product list: " << endl;
        for (int j = 0; j < data.number_products; ++j)
            if (x[j].get(GRB_DoubleAttr_X) > 0.5) {
                x_sol[j] = 1;
                cout << j << " ";
            }
            else x_sol[j] = 0;
        cout << endl;

        cout << "Conic obj = " << std::setprecision(5) << fixed << model.get(GRB_DoubleAttr_ObjVal) << endl;
        master_obj = calculate_master_obj(data, x_sol);
        cout << "Master obj = " << std::setprecision(5) << fixed << master_obj << endl;

        //check time
        auto time_now = std::chrono::steady_clock::now(); //get now time
        std::chrono::duration<double> total_time = time_now - start;
        cout << "Time now: " << total_time.count() << endl;
        cout << "--- --- --- --- --- --- ---" << endl;
    }
    else {
        cout << "No solution found..." << endl;
        auto end = chrono::steady_clock::now();
        elapsed_seconds = end - start;
        time_for_solve = elapsed_seconds.count();
    }

    auto end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    time_for_solve = elapsed_seconds.count();

    ofstream report_results(out_res_csv, ofstream::out);
    report_results.precision(10);
    report_results << model.get(GRB_DoubleAttr_ObjVal) << " " << master_obj << " " << time_for_solve << endl;
    for (int j = 0; j < data.number_products; ++j)
        if (x_sol[j] == 1)
            report_results << j << " ";
    report_results.close();
}

void CBSen::callback() {
    try {
        if (where == GRB_CB_MIPSOL) {
            double* initial_x = new double[products];
            double* initial_y = new double[customers];
            initial_x = getSolution(x, products);
            initial_y = getSolution(y, customers);

            vector<double> initial_denominator(customers, 0);
            for (int i = 0; i < customers; ++i) {
                initial_denominator[i] += noPay[i];
                for (int j = 0; j < products; ++j)
                    initial_denominator[i] += initial_x[j] * util[i][j];
            }

            //Calculate subgradient y
            vector<double> partitial_y(customers, 0);
            for (int i = 0; i < customers; ++i)
                partitial_y[i] = 1 / initial_denominator[i];

            vector<vector<double>> subgradient_y(customers);
            for (int i = 0; i < customers; ++i)
                subgradient_y[i].resize(products, 0);

            for (int i = 0; i < customers; ++i)
                for (int j = 0; j < products; ++j)
                    subgradient_y[i][j] -= util[i][j] / (initial_denominator[i] * initial_denominator[i]);

            //cout << "Outer-cuts\n" << endl;
            for (int i = 0; i < customers; ++i) {
                GRBLinExpr grad;
                for (int j = 0; j < products; ++j)
                    grad += subgradient_y[i][j] * (x[j] - initial_x[j]);

                addLazy(y[i] >= partitial_y[i] + grad);
                //cout << "Outer-cuts " << i << "\n" << endl;
            }

            ////cout << "Calculate total utility\n" << endl;
            //vector<double> sum_uti_customer(customers, 0);
            //for (int i = 0; i < customers; ++i) {
            //	sum_uti_customer[i] += noPay[i];
            //	for (int j = 0; j < products; ++j)
            //		sum_uti_customer[i] += util[i][j];
            //}

            ////cout << "Submodular-cuts\n" << endl;
            //for (int i = 0; i < customers; ++i) {
            //	if (initial_y[i] < partitial_y[i]) {
            //		GRBLinExpr submodular_cut_a_z, submodular_cut_b_z;
            //		for (int j = 0; j < products; ++j)
            //			if (initial_x[j] == 1) {
            //				submodular_cut_a_z += (1 - x[j]) * util[i][j] /
            //					(sum_uti_customer[i] * (sum_uti_customer[i] - util[i][j]));
            //				submodular_cut_b_z += (1 - x[j]) * util[i][j] /
            //					(initial_denominator[i] * (initial_denominator[i] - util[i][j]));
            //			}
            //			else {
            //				submodular_cut_a_z -= x[j] * util[i][j] /
            //					(initial_denominator[i] * (initial_denominator[i] + util[i][j]));
            //				submodular_cut_b_z -= x[j] * util[i][j] /
            //					(noPay[i] * (noPay[i] + util[i][j]));
            //			}

            //		submodular_cut_a_z += partitial_y[i];
            //		addLazy(y[i] >= submodular_cut_a_z);
            //		submodular_cut_b_z += partitial_y[i];
            //		addLazy(y[i] >= submodular_cut_b_z);
            //		//cout << "Submodular-cuts " << i << "\n" << endl;
            //	}
            //}
        }
    }
    catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Error during callback" << endl;
    }
}
