#include "CuttingPlane.h"
#include <chrono>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <cassert>

CuttingPlane::CuttingPlane() {

}

CuttingPlane::CuttingPlane(Data data, double time_limit, string outfile) {
	this->data = data;
	this->time_limit = time_limit;
	this->out_res_csv = outfile;
}

vector<double> CuttingPlane::calculate_y(Data data, vector<int> x, vector<double> alpha) {
	vector<double> y(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		double tmp_y = alpha[i] * data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j)
			tmp_y += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
		y[i] = log(tmp_y);
	}

	return y;
}

vector<double> CuttingPlane::calculate_z(Data data, vector<int> x) {
	vector<double> z(data.number_customers);
	for (int i = 0; i < data.number_customers; ++i) {
		double tmp_z = data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j)
			tmp_z += x[j] * data.utilities[i][j];
		z[i] = -log(tmp_z);
	}

	return z;
}

double  CuttingPlane::calculate_original_obj(Data data, vector<int> x, vector<double> alpha) {
	double obj = 0;
	for (int i = 0; i < data.number_customers; ++i) {
		double ts = alpha[i] * data.no_purchase[i], ms = data.no_purchase[i];
		for (int j = 0; j < data.number_products; ++j) {
			ts += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
			ms += x[j] * data.utilities[i][j];
		}
		obj += ts / ms;
	}

	return obj;
}

double CuttingPlane::calculate_master_obj(Data data, vector<int> x) {
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

double CuttingPlane::calculate_original_obj_tmp(Data data, vector<int> x, vector<double> alpha, int candidate) {
	double obj = 0;
	for (int i = 0; i < data.number_customers; ++i) {
		double ts = alpha[i] * data.no_purchase[i] + (alpha[i] - data.revenue[i][candidate]) * data.utilities[i][candidate],
			ms = data.no_purchase[i] + data.utilities[i][candidate];
		for (int j = 0; j < data.number_products; ++j) {
			ts += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
			ms += x[j] * data.utilities[i][j];
		}

		obj += ts / ms;
	}

	return obj;
}

double CuttingPlane::calculate_optimal_bound_y(Data data, int i, double alpha) {
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

	GRBLinExpr obj;
	for (int j = 0; j < data.number_products; ++j)
		obj += data.utilities[i][j] * x[j] * (alpha - data.revenue[i][j]);
	model.setObjective(obj, GRB_MAXIMIZE);

	model.set(GRB_IntParam_OutputFlag, 0);

	model.optimize();

	return model.get(GRB_DoubleAttr_ObjVal) + alpha * data.no_purchase[i];
}

double CuttingPlane::calculate_optimal_bound_z(Data data, int i) {
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

	GRBLinExpr obj;
	for (int j = 0; j < data.number_products; ++j)
		obj += data.utilities[i][j] * x[j];
	model.setObjective(obj, GRB_MAXIMIZE);

	model.set(GRB_IntParam_OutputFlag, 0);

	model.optimize();

	return model.get(GRB_DoubleAttr_ObjVal) + data.no_purchase[i];
}

double CuttingPlane::calculate_bound_y(Data data, int i, double alpha) {
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
	//cout << "Chosen y = " << number_chosen << endl;
	double bound = 0;
	vector<pair<double, int>> c(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		c[j] = make_pair(data.utilities[i][j] * (alpha - data.revenue[i][j]), j);

	sort(c.begin(), c.end(), greater<pair<double, int>>());
	if (number_chosen >= data.number_products) number_chosen = data.number_products;
	for (int j = 0; j < number_chosen; ++j)
		bound += c[j].first;
	bound += alpha * data.no_purchase[i];
	//cout << "Bound y = [" << log(alpha * data.no_purchase[i]) << "," << log(bound) << "]" << endl;
	//cout << "Bound u = [" << alpha * data.no_purchase[i] << "," << bound << "]" << endl;
	return bound;
}

double CuttingPlane::calculate_bound_z(Data data, int i) {
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
	//cout << "Chosen z = " << number_chosen << endl;
	double bound = 0;
	vector<pair<double, int>> c(data.number_products);
	for (int j = 0; j < data.number_products; ++j)
		c[j] = make_pair(data.utilities[i][j], j);

	sort(c.begin(), c.end(), greater<pair<double, int>>());
	if (number_chosen >= data.number_products) number_chosen = data.number_products;
	for (int j = 0; j < number_chosen; ++j)
		bound += c[j].first;
	bound += data.no_purchase[i];
	//cout << "Bound z = [" << -log(bound) << "," << -log(data.no_purchase[i]) << "]" << endl;
	return bound;
}

vector<int> CuttingPlane::greedy(Data data, vector<double> alpha) {
	auto start = chrono::steady_clock::now();
	vector<int> chosen(data.number_products, 0);
	double obj = calculate_original_obj(data, chosen, alpha);
	//cout << "Obj: " << obj << endl;
	vector<double> set_capacity(data.number_sets, 0);

	while (true) {
		double obj_tmp = 0;
		int inserted_product = -1;
		double min = 999999;
		for (int j = 0; j < data.number_products; ++j)
			if (chosen[j] == 0) {
				int check_in_set = 0;
				int check_valid = 0;
				for (int s = 0; s < data.number_sets; ++s) {
					if (data.in_set[j][s] == 1) {
						check_in_set++;
						if (set_capacity[s] + data.cost[j] <= data.capacity_each_set)
							check_valid++;
						else break;
					}
				}
				if (check_in_set != 0 && check_in_set == check_valid) {
					obj_tmp = calculate_original_obj_tmp(data, chosen, alpha, j);
					if (min > obj_tmp) {
						min = obj_tmp;
						inserted_product = j;
					}
				}
				if (check_in_set == 0) {
					obj_tmp = calculate_original_obj_tmp(data, chosen, alpha, j);
					if (min > obj_tmp) {
						min = obj_tmp;
						inserted_product = j;
					}
				}
			}

		if (inserted_product != -1) {
			chosen[inserted_product] = 1;
			obj = calculate_original_obj(data, chosen, alpha);
			//cout << "Best-Insertion: " << inserted_product << " " << obj << endl;
			for (int s = 0; s < data.number_sets; ++s)
				if (data.in_set[inserted_product][s] == 1)
					set_capacity[s] += data.cost[inserted_product];
		}
		else break;
	}

	auto end = chrono::steady_clock::now();
	std::chrono::duration<double> greedy_time = end - start;

	cout << "Greedy solution: " << endl;
	for (int j = 0; j < data.number_products; j++)
		if (chosen[j] == 1)
			cout << j << " ";
	cout << endl;
	for (int s = 0; s < data.number_sets; ++s)
		cout << set_capacity[s] << " ";

	cout << endl;
	cout << "Sub obj = " << obj << endl;
	cout << "Master obj = " << calculate_master_obj(data, chosen) << endl;
	cout << "Time: " << greedy_time.count() << endl;

	return chosen;
}

void CuttingPlane::solve(Data data) {
	auto start = chrono::steady_clock::now();

	vector<double> alpha(data.number_customers, -1);
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			if (data.revenue[i][j] > alpha[i])
				alpha[i] = data.revenue[i][j];

	//Calculate initial_x, initial_y, initial_z
	//vector<int> initial_x(data.number_products, 0);
	vector<int> initial_x = greedy(data, alpha);
	vector<double> initial_y = calculate_y(data, initial_x, alpha);
	vector<double> initial_z = calculate_z(data, initial_x);
	//vector<double> initial_theta(data.number_customers, 0);
	//for (int i = 0; i < data.number_customers; ++i)
	//	initial_theta[i] = exp(initial_y[i] + initial_z[i]);

	GRBEnv env = GRBEnv(true);
	env.start();

	GRBModel model = GRBModel(env);

	//cout << "Decison variables : x_j\n" << endl;
	GRBVar* x;
	x = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		x[j] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(j));

	//cout << "Slack variables : y_i\n" << endl;
	vector<double> lower_bound_y;
	vector<double> upper_bound_y;
	GRBVar* y;
	GRBVar* u;
	y = new GRBVar[data.number_customers];
	u = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i) {
		double bound = calculate_optimal_bound_y(data, i, alpha[i]);
		//double bound = calculate_bound_y(data, i, alpha[i]);
		y[i] = model.addVar(log(alpha[i] * data.no_purchase[i]), log(bound), 0, GRB_CONTINUOUS, "y_" + to_string(i));
		u[i] = model.addVar(alpha[i] * data.no_purchase[i], bound, 0, GRB_CONTINUOUS, "u_" + to_string(i));
		lower_bound_y.push_back(log(alpha[i] * data.no_purchase[i]));
		upper_bound_y.push_back(log(bound));
	}

	//cout << "Slack variables : z_i\n" << endl;
	vector<double> lower_bound_z;
	vector<double> upper_bound_z;
	GRBVar* z;
	z = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i) {
		double bound = calculate_optimal_bound_z(data, i);
		//double bound = calculate_bound_z(data, i);
		z[i] = model.addVar(-log(bound), -log(data.no_purchase[i]), 0, GRB_CONTINUOUS, "z_" + to_string(i));
		lower_bound_z.push_back(-log(bound));
		upper_bound_z.push_back(-log(data.no_purchase[i]));
	}

	GRBVar* theta;
	theta = new GRBVar[data.number_customers];
	for (int c = 0; c < data.number_customers; ++c)
		theta[c] = model.addVar(exp(lower_bound_y[c] + lower_bound_z[c]), exp(upper_bound_y[c] + upper_bound_z[c]), 0, GRB_CONTINUOUS, "theta_" + to_string(c));

	for (int i = 0; i < data.number_customers; ++i) {
		model.addGenConstrExp(y[i], u[i]);
		GRBLinExpr sum_x;
		for (int j = 0; j < data.number_products; ++j)
			sum_x += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
		sum_x += alpha[i] * data.no_purchase[i];
		model.addConstr(u[i] >= sum_x);
	}

	//Set capacity constraints
	for (int s = 0; s < data.number_sets; ++s) {
		GRBLinExpr sum;
		for (int j = 0; j < data.number_products; ++j)
			if (data.in_set[j][s] == 1)
				sum += data.cost[j] * x[j];
		model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
	}

	//Objective
	GRBLinExpr obj;
	for (int c = 0; c < data.number_customers; ++c)
		obj += theta[c];
	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	int num_iterative = 0;
	double stop_param = 1e-4;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;
	double best_sub_obj = 0;
	vector<int> best_x = initial_x;
	double best_obj = calculate_master_obj(data, best_x);

	//cout << "Add cut iteratively" << endl;
	while (sub_obj > obj_val_cplex + stop_param) {
		//compute gradient e^{y+z} at initial_x, initial_y, initial_z and set up constraints related to theta
		for (int c = 0; c < data.number_customers; ++c)
			model.addConstr(theta[c] >= exp(initial_y[c] + initial_z[c]) * (1 + y[c] - initial_y[c] + z[c] - initial_z[c]), "ct_sub_gradient_e_" + to_string(c));

		//compute gradient e^{-z} at initial_x, initial_y, initial_z and set up constraints related to e^{-z}
		for (int i = 0; i < data.number_customers; ++i) {
			GRBLinExpr sum = 0;
			for (int j = 0; j < data.number_products; ++j)
				sum += x[j] * data.utilities[i][j];
			sum += data.no_purchase[i];
			model.addConstr(exp(-initial_z[i]) * (1 - z[i] + initial_z[i]) <= sum, "ct_sub_gradient_z_" + to_string(i));
		}

		////cout << "Check theta" << endl;
		//vector<int> check_theta;
		//for (int i = 0; i < data.number_customers; ++i)
		//	if (initial_theta[i] < exp(initial_y[i] + initial_z[i]))
		//		check_theta.push_back(i);

		////cout << "Check e^z" << endl;
		//vector<int> check_z;
		//for (int i = 0; i < data.number_customers; ++i) {
		//	double sum_x = 0;
		//	for (int j = 0; j < data.number_products; ++j)
		//		sum_x += initial_x[j] * data.utilities[i][j];
		//	sum_x += data.no_purchase[i];
		//	if (exp(-initial_z[i]) > sum_x)
		//		check_z.push_back(i);
		//}

		//if (check_theta.size() >= 1) {
		//	//cout << "compute gradient e^ {y + z} at initial_x, initial_y, initial_z and set up constraints related to theta" << endl;
		//	for (int i = 0; i < check_theta.size(); ++i) {
		//		GRBLinExpr sum;
		//		sum += exp(initial_y[check_theta[i]] + initial_z[check_theta[i]]) * (1 + y[check_theta[i]] - initial_y[check_theta[i]] + z[check_theta[i]] - initial_z[check_theta[i]]);
		//		model.addConstr(theta[check_theta[i]] >= sum);
		//		//cout << "Add lazy cut theta\n";
		//	}
		//}

		//if (check_z.size() >= 1) {
		//	//cout << "compute gradient e^ {-z} at initial_x, initial_y, initial_z and set up constraints related to e^ {-z}" << endl;
		//	for (int i = 0; i < check_z.size(); ++i) {
		//		GRBLinExpr sum = 0;
		//		for (int j = 0; j < data.number_products; ++j)
		//			sum += x[j] * data.utilities[check_z[i]][j];
		//		sum += data.no_purchase[check_z[i]];
		//		model.addConstr(exp(-initial_z[check_z[i]]) * (1 - z[check_z[i]] + initial_z[check_z[i]]) <= sum);
		//		//cout << "Add lazy cut z\n";
		//	}
		//}

		//solve
		num_iterative++;
		cout << "Remaining time: " << run_time << endl;

		model.write("cutting_plane.lp");
		model.set(GRB_DoubleParam_TimeLimit, run_time);
		model.set(GRB_DoubleParam_MIPGap, 1e-3);
		model.set(GRB_IntParam_FuncPieces, -1);
		model.set(GRB_DoubleParam_FuncPieceError, 1e-1);
		//model.set(GRB_DoubleParam_FuncPieceError, 1e-2);
		//model.set(GRB_DoubleParam_FuncPieceError, 1e-3);
		model.set(GRB_IntParam_OutputFlag, 0);
		model.set(GRB_IntParam_Threads, 8);

		model.optimize();

		if (model.get(GRB_IntAttr_SolCount) > 0) {
			cout << "\nIteration " << num_iterative << endl;
			//update obj, variables
			obj_val_cplex = model.get(GRB_DoubleAttr_ObjVal);
			cout << "\nResult product list: " << endl;
			for (int j = 0; j < data.number_products; ++j)
				if (x[j].get(GRB_DoubleAttr_X) > 0.5) {
					initial_x[j] = 1;
					cout << j << " ";
				}
				else initial_x[j] = 0;
			cout << endl;

			for (int i = 0; i < data.number_customers; ++i) {
				initial_y[i] = y[i].get(GRB_DoubleAttr_X);
				//cout << initial_y[i] << " ";
				initial_z[i] = z[i].get(GRB_DoubleAttr_X);
				//cout << initial_z[i] << endl;
				//initial_theta[i] = theta[i].get(GRB_DoubleAttr_X);
			}

			//check the in equation related to theta_i and e^{y_i + z_i} for next iteration
			sub_obj = 0;

			for (int i = 0; i < data.number_customers; ++i)
				sub_obj += exp(initial_y[i] + initial_z[i]);

			cout << "Sub obj = " << std::setprecision(7) << fixed << sub_obj << endl;
			cout << "Gurobi obj = " << std::setprecision(7) << fixed << obj_val_cplex << endl;
			master_obj_val = calculate_master_obj(data, initial_x);
			cout << "Master obj = " << std::setprecision(7) << fixed << master_obj_val << endl;

			if (master_obj_val >= best_obj) {
				best_obj = master_obj_val;
				best_x = initial_x;
				best_sub_obj = obj_val_cplex;
			}

			//check time
			auto time_now = std::chrono::steady_clock::now(); //get now time
			std::chrono::duration<double> after_cut = time_now - start;
			cout << "Time now: " << after_cut.count() << endl;
			cout << "--- --- --- --- --- --- ---" << endl;

			if (after_cut.count() > time_limit) break;
			run_time = time_limit - after_cut.count();
		}
		else {
			cout << "Iteration " << num_iterative << ". No solution found..." << endl;
			auto end = chrono::steady_clock::now();
			chrono::duration<double> elapsed_seconds = end - start;
			time_for_solve = elapsed_seconds.count();
			break;
		}
	}
	auto end = chrono::steady_clock::now();
	chrono::duration<double> total_time = end - start;
	time_for_solve = total_time.count();

	cout << "\nObjective value: " << setprecision(5) << best_obj << endl;
	cout << "Solution: ";
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			cout << j << " ";
	cout << "\nTotal time: " << time_for_solve << " seconds" << endl;

	ofstream report_results(out_res_csv, ofstream::out);
	report_results.precision(10);
	report_results << best_sub_obj << " " << best_obj << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}

void CuttingPlane::solve_multicut(Data data, int number_cuts) {
	auto start = chrono::steady_clock::now();

	vector<double> alpha(data.number_customers, -1);
	for (int i = 0; i < data.number_customers; ++i)
		for (int j = 0; j < data.number_products; ++j)
			if (data.revenue[i][j] > alpha[i])
				alpha[i] = data.revenue[i][j];

	//Calculate initial_x, initial_y, initial_z
	//vector<int> initial_x(data.number_products, 0);
	vector<int> initial_x = greedy(data, alpha);
	vector<double> initial_y = calculate_y(data, initial_x, alpha);
	vector<double> initial_z = calculate_z(data, initial_x);
	vector<double> initial_theta(number_cuts,0);
	for (int k = 0; k < number_cuts; ++k)
		for (int c = 0; c < data.number_customers; ++c)
			if (c % number_cuts == k)
				initial_theta[k] += exp(initial_y[c] + initial_z[c]);

	GRBEnv env = GRBEnv(true);
	env.start();

	GRBModel model = GRBModel(env);

	//cout << "Decison variables : x_j\n" << endl;
	GRBVar* x;
	x = new GRBVar[data.number_products];
	for (int j = 0; j < data.number_products; ++j)
		x[j] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + to_string(j));

	//cout << "Slack variables : y_i\n" << endl;
	vector<double> lower_bound_y;
	vector<double> upper_bound_y;
	GRBVar* y;
	GRBVar* u;
	y = new GRBVar[data.number_customers];
	u = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i) {
		double bound = calculate_optimal_bound_y(data, i, alpha[i]);
		//double bound = calculate_bound_y(data, i, alpha[i]);
		y[i] = model.addVar(log(alpha[i] * data.no_purchase[i]), log(bound), 0, GRB_CONTINUOUS, "y_" + to_string(i));
		u[i] = model.addVar(alpha[i] * data.no_purchase[i], bound, 0, GRB_CONTINUOUS, "u_" + to_string(i));
		lower_bound_y.push_back(log(alpha[i] * data.no_purchase[i]));
		upper_bound_y.push_back(log(bound));
	}

	//cout << "Slack variables : z_i\n" << endl;
	vector<double> lower_bound_z;
	vector<double> upper_bound_z;
	GRBVar* z;
	z = new GRBVar[data.number_customers];
	for (int i = 0; i < data.number_customers; ++i) {
		double bound = calculate_optimal_bound_z(data, i);
		//double bound = calculate_bound_z(data, i);
		z[i] = model.addVar(-log(bound), -log(data.no_purchase[i]), 0, GRB_CONTINUOUS, "z_" + to_string(i));
		lower_bound_z.push_back(-log(bound));
		upper_bound_z.push_back(-log(data.no_purchase[i]));
	}

	vector<double> lower_bound_cut(number_cuts);
	vector<double> upper_bound_cut(number_cuts);
	for (int c = 0; c < data.number_customers; ++c)
		for (int k = 0; k < number_cuts; ++k)
			if (c % number_cuts == k) {
				lower_bound_cut[k] += exp(lower_bound_y[c] + lower_bound_z[c]);
				upper_bound_cut[k] += exp(upper_bound_y[c] + upper_bound_z[c]);
			}

	GRBVar* theta_cut;
	theta_cut = new GRBVar[number_cuts];
	for (int k = 0; k < number_cuts; ++k)
		theta_cut[k] = model.addVar(lower_bound_cut[k],upper_bound_cut[k], 0, GRB_CONTINUOUS, "theta_cut_" + to_string(k));

	for (int i = 0; i < data.number_customers; ++i) {
		model.addGenConstrExp(y[i], u[i]);
		GRBLinExpr sum_x;
		for (int j = 0; j < data.number_products; ++j)
			sum_x += (alpha[i] - data.revenue[i][j]) * x[j] * data.utilities[i][j];
		sum_x += alpha[i] * data.no_purchase[i];
		model.addConstr(u[i] >= sum_x);
	}

	//Set capacity constraints
	for (int s = 0; s < data.number_sets; ++s) {
		GRBLinExpr sum;
		for (int j = 0; j < data.number_products; ++j)
			if (data.in_set[j][s] == 1)
				sum += data.cost[j] * x[j];
		model.addConstr(sum <= data.capacity_each_set, "ct_set_cap" + to_string(s));
	}

	//Objective
	GRBLinExpr obj;
	for (int k = 0; k < number_cuts; ++k)
		obj += theta_cut[k];
	model.setObjective(obj, GRB_MINIMIZE);

	auto time_before_cut = chrono::steady_clock::now();
	chrono::duration<double> before_cut = time_before_cut - start;

	double run_time = time_limit - before_cut.count();

	int num_iterative = 0;
	double stop_param = 1e-4;
	double sub_obj = 1.0;
	double obj_val_cplex = 0.0;
	double best_sub_obj = 0;
	vector<int> best_x = initial_x;
	double best_obj = calculate_master_obj(data, best_x);

	//cout << "Add cut iteratively" << endl;
	while (sub_obj > obj_val_cplex + stop_param) {
		////cout << "Check theta" << endl;
		//vector<int> check_theta;
		//for (int i = 0; i < number_cuts; ++i) {
		//	double sum = 0;
		//	for (int c = 0; c < data.number_customers; ++c)
		//		if (c % number_cuts == i)
		//			sum += exp(initial_y[c] + initial_z[c]);
		//	if (initial_theta[i] < sum)
		//		check_theta.push_back(i);
		//}

		////cout << "Check e^z" << endl;
		//vector<int> check_z;
		//for (int i = 0; i < data.number_customers; ++i) {
		//	double sum_x = 0;
		//	for (int j = 0; j < data.number_products; ++j)
		//		sum_x += initial_x[j] * data.utilities[i][j];
		//	sum_x += data.no_purchase[i];
		//	if (exp(-initial_z[i]) > sum_x)
		//		check_z.push_back(i);
		//}

		//if (check_theta.size() >= 1) {
		//	//cout << "compute gradient e^ {y + z} at initial_x, initial_y, initial_z and set up constraints related to theta" << endl;
		//	for (int i = 0; i < check_theta.size(); ++i) {
		//		GRBLinExpr sum;
		//		for (int c = 0; c < data.number_customers; ++c)
		//			if (c % number_cuts == check_theta[i])
		//				sum += exp(initial_y[c] + initial_z[c]) * (1 + y[c] - initial_y[c] + z[c] - initial_z[c]);
		//		model.addConstr(theta_cut[check_theta[i]] >= sum);
		//		//cout << "Add lazy cut theta\n";
		//	}
		//}

		//if (check_z.size() >= 1) {
		//	//cout << "compute gradient e^ {-z} at initial_x, initial_y, initial_z and set up constraints related to e^ {-z}" << endl;
		//	for (int i = 0; i < check_z.size(); ++i) {
		//		GRBLinExpr sum = 0;
		//		for (int j = 0; j < data.number_products; ++j)
		//			sum += x[j] * data.utilities[check_z[i]][j];
		//		sum += data.no_purchase[check_z[i]];
		//		model.addConstr(exp(-initial_z[check_z[i]]) * (1 - z[check_z[i]] + initial_z[check_z[i]]) <= sum);
		//		//cout << "Add lazy cut z\n";
		//	}
		//}

		//compute gradient e^{y+z} at initial_x, initial_y, initial_z and set up constraints related to theta
		for (int k = 0; k < number_cuts; ++k) {
			GRBLinExpr sum = 0;
			for (int c = 0; c < data.number_customers; ++c)
				if (c % number_cuts == k)
					sum += exp(initial_y[c] + initial_z[c]) * (1 + y[c] - initial_y[c] + z[c] - initial_z[c]);
			model.addConstr(theta_cut[k] >= sum, "ct_sub_gradient_theta_cut_" + to_string(k));
		}

		//compute gradient e^{-z} at initial_x, initial_y, initial_z and set up constraints related to e^{-z}
		for (int i = 0; i < data.number_customers; ++i) {
			GRBLinExpr sum = 0;
			for (int j = 0; j < data.number_products; ++j)
				sum += x[j] * data.utilities[i][j];
			sum += data.no_purchase[i];
			model.addConstr(exp(-initial_z[i]) * (1 - z[i] + initial_z[i]) <= sum, "ct_sub_gradient_z_" + to_string(i));
		}

		//solve
		num_iterative++;
		cout << "Remaining time: " << run_time << endl;

		model.write("cutting_plane.lp");
		model.set(GRB_DoubleParam_TimeLimit, run_time);
		model.set(GRB_DoubleParam_MIPGap, 1e-3);
		model.set(GRB_IntParam_FuncPieces, -1);
		model.set(GRB_DoubleParam_FuncPieceError, 1e-1);
		//model.set(GRB_DoubleParam_FuncPieceError, 1e-2);
		//model.set(GRB_DoubleParam_FuncPieceError, 1e-3);
		//model.set(GRB_IntParam_OutputFlag, 0);

		model.optimize();

		if (model.get(GRB_IntAttr_SolCount) > 0) {
			cout << "\nIteration " << num_iterative << endl;
			//update obj, variables
			obj_val_cplex = model.get(GRB_DoubleAttr_ObjVal);
			cout << "\nResult product list: " << endl;
			for (int j = 0; j < data.number_products; ++j)
				if (x[j].get(GRB_DoubleAttr_X) > 0.5) {
					initial_x[j] = 1;
					cout << j << " ";
				}
				else initial_x[j] = 0;
			cout << endl;

			for (int i = 0; i < data.number_customers; ++i) {
				initial_y[i] = y[i].get(GRB_DoubleAttr_X);
				//cout << initial_y[i] << " ";
				initial_z[i] = z[i].get(GRB_DoubleAttr_X);
				//cout << initial_z[i] << endl;
			}
			for (int i = 0; i < number_cuts; ++i)
				initial_theta[i] = theta_cut[i].get(GRB_DoubleAttr_X);

			//check the in equation related to theta_i and e^{y_i + z_i} for next iteration
			sub_obj = 0;
			for (int i = 0; i < data.number_customers; ++i)
				sub_obj += exp(initial_y[i] + initial_z[i]);

			//sub_obj = calculate_original_obj(data, initial_x, alpha);

			cout << "Sub obj = " << std::setprecision(7) << fixed << sub_obj << endl;
			cout << "Gurobi obj = " << std::setprecision(7) << fixed << obj_val_cplex << endl;
			master_obj_val = calculate_master_obj(data, initial_x);
			cout << "Master obj = " << std::setprecision(7) << fixed << master_obj_val << endl;

			if (master_obj_val >= best_obj) {
				best_obj = master_obj_val;
				best_x = initial_x;
				best_sub_obj = obj_val_cplex;
			}

			//check time
			auto time_now = std::chrono::steady_clock::now(); //get now time
			std::chrono::duration<double> after_cut = time_now - start;
			cout << "Time now: " << after_cut.count() << endl;
			cout << "--- --- --- --- --- --- ---" << endl;

			if (after_cut.count() > time_limit) break;
			run_time = time_limit - after_cut.count();
		}
		else {
			cout << "Iteration " << num_iterative << ". No solution found..." << endl;
			auto end = chrono::steady_clock::now();
			chrono::duration<double> elapsed_seconds = end - start;
			time_for_solve = elapsed_seconds.count();
			break;
		}
	}
	auto end = chrono::steady_clock::now();
	chrono::duration<double> total_time = end - start;
	time_for_solve = total_time.count();

	cout << "\nObjective value: " << setprecision(5) << best_obj << endl;
	cout << "Solution: ";
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			cout << j << " ";
	cout << "\nTotal time: " << time_for_solve << " seconds" << endl;

	ofstream report_results(out_res_csv, ofstream::out);
	report_results.precision(10);
	report_results << best_sub_obj << " " << best_obj << " " << time_for_solve << endl;
	for (int j = 0; j < data.number_products; ++j)
		if (best_x[j] == 1)
			report_results << j << " ";
	report_results.close();
}
