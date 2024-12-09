from mosek.fusion import *
import mosek.fusion.pythonic
import numpy as np
import sys
import threading
import time

class Data:
    def __init__(self):
        self.number_products = 0
        self.number_customers = 0
        self.number_sets = 0
        self.capacity_each_set = 0
        self.utilities = []
        self.revenue = []
        self.cost = []
        self.in_set = []
        self.no_purchase = []

    def read_data(self, path, noPay):
        file = open(path, 'r').readlines()
        self.number_products = int(file[0])
        #print(self.number_products)
        self.number_customers = int(file[1])
        #print(self.number_customers)
        self.number_sets = int(file[2])
        #print(self.number_sets)
        self.capacity_each_set = int(file[3])
        #print(self.capacity_each_set)
        line = 4
        for _ in range(self.number_customers):
            self.utilities.append(list(map(float, file[line].split())))
            line += 1
        #print(self.utilities)
        for _ in range(self.number_customers):
            self.revenue.append(list(map(float, file[line].split())))
            line += 1
        #print(self.revenue)
        self.cost = list(map(float, file[line].split()))
        line += 1
        #print(self.cost)
        for _ in range(self.number_products):
            self.in_set.append(list(map(int, file[line].split())))
            line += 1
        #print(self.in_set)
        for _ in range(self.number_customers):
            self.no_purchase.append(noPay)

class Solver:
    def __init__(self):
        self.output =''
        self.log_file =''
        self.out_res_csv=''
        self.data = Data()
        self.time_limit = 0
        self.master_obj_val = 0
        self.gap = 0
        self.time_for_solve = 0

    def calculate_master_obj(self, data, x):
        obj = 0
        for i in range(data.number_customers):
            ts = 0
            ms = data.no_purchase[i]
            for j in range(data.number_products):
                ts += data.revenue[i][j] * x[j] * data.utilities[i][j]
                ms += x[j] * data.utilities[i][j]
            obj += ts / ms
        return obj

    def calculate_optimal_bound_y(self, data, i, alpha):
        with Model('bound_y') as M:
            x = M.variable('x', data.number_products, Domain.binary())
            #Set capacity constraints
            for s in range(data.number_sets):
                sum_x = 0
                for j in range(data.number_products):
                    if data.in_set[j][s] == 1:
                        sum_x += data.cost[j] * x[j]
                M.constraint(sum_x <= data.capacity_each_set)

            #Objective
            obj = 0
            for j in range(data.number_products):
                obj += data.utilities[i][j] * (alpha - data.revenue[i][j]) * x[j]
            M.objective('Objective', ObjectiveSense.Maximize, obj)
            M.solve()
            return M.primalObjValue() + alpha * data.no_purchase[i]

    def calculate_optimal_bound_z(self, data, i):
        with Model('bound_z') as M:
            x = M.variable('x', data.number_products, Domain.binary())
            #Set capacity constraints
            for s in range(data.number_sets):
                sum_x = 0
                for j in range(data.number_products):
                    if data.in_set[j][s] == 1:
                        sum_x += data.cost[j] * x[j]
                M.constraint(sum_x <= data.capacity_each_set)

            #Objective
            obj = 0
            for j in range(data.number_products):
                obj += data.utilities[i][j] * x[j]
            M.objective('Objective', ObjectiveSense.Maximize, obj)
            M.solve()
            return M.primalObjValue() + data.no_purchase[i]

    def next_approximate_point(self, b, epsilon):
        x = b + 0.001
        while np.exp(b) + (np.exp(x) - np.exp(b)) / (x - b) * (np.log((np.exp(x) - np.exp(b)) / (x - b)) - b) - (np.exp(x) - np.exp(b)) / (x - b) <= epsilon:
            x += 0.001
        return x
    
    def optimal_sub_intervals(self, data, alpha, epsilon):
        c = []
        for _ in range(data.number_customers):
            c.append(list())
        for i in range(data.number_customers):
            c[i].append(np.log(alpha[i] * data.no_purchase[i]))
            upper = np.log(self.calculate_optimal_bound_y(data, i, alpha[i]))
            while True:
                x = self.next_approximate_point(c[i][len(c[i])-1], epsilon)
                if x >= upper: 
                    break
                else:
                    c[i].append(x)
            c[i].append(upper)
        return c
    
    def solve(self, data, name):
        alpha = []
        for _ in range(data.number_customers):
            alpha.append(-1)
        for i in range(data.number_customers):
            for j in range(data.number_products):
                if data.revenue[i][j] > alpha[i]:
                    alpha[i] = data.revenue[i][j]

        c = self.optimal_sub_intervals(data, alpha, 0.001)

        number_sub_intervals = []
        for _ in range(data.number_customers):
            number_sub_intervals.append(0)
        for i in range(data.number_customers):
            number_sub_intervals[i] = len(c[i]) - 1
        print(number_sub_intervals)
        max_number_sub_intervals = 0
        for i in range(data.number_customers):
            if number_sub_intervals[i] > max_number_sub_intervals:
                max_number_sub_intervals = number_sub_intervals[i]

        # bound_y = []
        # upper_bound_y = []
        # lower_bound_y = []
        # bound_z = []
        # upper_bound_z = []
        # lower_bound_z = []
        # for i in range(data.number_customers):
        #     bound_y.append(self.calculate_optimal_bound_y(data, i, alpha[i]))
        #     bound_z.append(self.calculate_optimal_bound_z(data, i))

        # for i in range(data.number_customers):
        #     lower_bound_y.append(np.log(alpha[i] * data.no_purchase[i]))
        #     upper_bound_y.append(np.log(bound_y[i]))
        #     lower_bound_z.append(-np.log(bound_z[i]))
        #     upper_bound_z.append(-np.log(data.no_purchase[i]))

        with Model('ao') as M:
            x = M.variable('x', data.number_products, Domain.binary())
            y = M.variable('y', data.number_customers, Domain.unbounded())
            z = M.variable('z', data.number_customers, Domain.unbounded())
            theta = M.variable('theta', data.number_customers, Domain.unbounded())
            r = M.variable('r', [data.number_customers,max_number_sub_intervals], Domain.binary())
            s = M.variable('s', [data.number_customers,max_number_sub_intervals], Domain.inRange(0.0,1.0))

            # for i in range(data.number_customers):
            #     M.constraint(y[i] <= upper_bound_y[i])
            #     M.constraint(y[i] >= lower_bound_y[i])
            #     M.constraint(z[i] <= upper_bound_z[i])
            #     M.constraint(z[i] >= lower_bound_z[i])
            #     M.constraint(theta[i] <= np.exp(upper_bound_y[i] + upper_bound_z[i]))
            #     M.constraint(theta[i] >= np.exp(lower_bound_y[i] + lower_bound_z[i]))

            #Constraints related to e^{y_i}
	        #exp(c[i][0]) = alpha * no_purchase[i] => remove form both sides
            for i in range(data.number_customers):
                sum_r = 0
                sum_x = 0
                for k in range(number_sub_intervals[i]):
                    sum_r += (np.exp(c[i][k+1]) - np.exp(c[i][k])) * r[i,k]
                for j in range(data.number_products):
                    sum_x += (alpha[i] - data.revenue[i][j]) * data.utilities[i][j] * x[j]
                M.constraint(sum_r >= sum_x)

            #Constraints related to y_i
            for i in range(data.number_customers):
                sum_r = 0
                for k in range (number_sub_intervals[i]):
                    sum_r += (c[i][k+1] - c[i][k]) * r[i,k]
                M.constraint(y[i] == c[i][0] + sum_r)

            #Constraints related to s_{ik} and r_{ik}
            for i in range(data.number_customers):
                for k in range(number_sub_intervals[i]):
                    M.constraint(r[i,k] >= s[i,k])
                    if k < number_sub_intervals[i] - 1:
                        M.constraint(s[i,k] >= s[i,k+1])
                        M.constraint(s[i,k] >= r[i,k+1])

            #Set capacity constraints
            for s in range(data.number_sets):
                sum_x = 0
                for j in range(data.number_products):
                    if data.in_set[j][s] == 1:
                        sum_x += data.cost[j] * x[j]
                M.constraint(sum_x <= data.capacity_each_set)

            #theta >= e^(y+z) and e^(-z) <= sum x
            for i in range(data.number_customers):
                M.constraint(Expr.vstack(theta[i],1,y[i]+z[i]), Domain.inPExpCone())
                sum_x = data.no_purchase[i]
                for j in range(data.number_products):
                    sum_x += data.utilities[i][j] * x[j]
                M.constraint(Expr.vstack(sum_x,1,-z[i]), Domain.inPExpCone())

            #Objective
            obj = 0
            for i in range(data.number_customers):
                obj += theta[i]
            M.objective('Objective', ObjectiveSense.Minimize, obj)

            M.setLogHandler(sys.stdout)
            #M.setSolverParam("numThreads", 8)
            M.setSolverParam('mioTolRelGap', 1e-4)

            timeout = 3600
            T = threading.Thread(target=M.solve)
            T0 = time.time()
            try:
                T.start()
                while True:
                    if not T.is_alive():
                        print("Solver terminated before anything happened!")
                        break
                    elif time.time() - T0 > timeout:
                        print("Solver terminated due to timeout!")
                        M.breakSolver()
                        break
            except KeyboardInterrupt:
                print("Signalling the solver that it can give up now!")
                M.breakSolver()
            finally:
                try:
                    T.join()
                except:
                    pass

            #M.solve()
            M.writeTask("C:/Users/user/mosek/ao.ptf")
            M.acceptedSolutionStatus(AccSolutionStatus.Feasible)

            sol_x = []
            sol_x = M.getVariable('x').level()
            for i in range(len(sol_x)):
                if sol_x[i] >= 0.5:
                    sol_x[i] = 1
                else:
                    sol_x[i] = 0
            
            output = open("C:/Users/user/mosek/Result/" + name + "_" + str(data.no_purchase[0]) + ".txt", "w+")
            output.write("%.7f %.7f %.5f" % (M.primalObjValue(),self.calculate_master_obj(data, sol_x),M.getSolverDoubleInfo("optimizerTime")))
            output.write("\n")
            for j in range(data.number_products):
                if sol_x[j] == 1:
                    output.write(str(j) + " ")
            output.close()

if __name__ == "__main__":
    capacity = ["5","10","20"]
    group = ["5","10","20"]
    noPay = [1,2]
    customers = ["5"]#,"10","15","20"]
    for n in range(len(customers)):
        for g in range(len(group)):
            for c in range(len(capacity)):
                for i in range(len(noPay)):
                    name = "100_" + customers[n] + "_overlap_" + group[g] + "_" + capacity[c]
                    path = "C:/Users/user/mosek/AO_data/" + name + ".dat"
                    data = Data()
                    data.read_data(path, noPay[i])
                    solver = Solver()
                    solver.solve(data, name)