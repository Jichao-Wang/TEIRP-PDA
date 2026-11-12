import numpy as np
from gurobipy import *
import docplex.mp.model as cpx

from case33 import case33
import matplotlib.pyplot as plt

import pandas as pd
import copy
import time
import os

global over_response_ratio_low  # 允许需求响应超调功率相对于ACE信号调节值的比例
global over_response_ratio_up  # 允许需求响应超调功率相对于ACE信号调节值的比例
global colors


def E(agg):
    # E(\eta_{rn})
    # return agg.alpha / (agg.alpha + agg.beta)
    return agg.alpha / (agg.alpha + agg.beta) * 0.2 - 0.1 + 1


def Var(agg):
    # D(\eta_{rn})
    return agg.alpha * agg.beta / ((agg.alpha + agg.beta) * (agg.alpha + agg.beta) * (agg.alpha + agg.beta + 1)) * 0.04


class AGG:
    def __init__(self, cap, alpha, beta, price_param):
        self.capacity = cap  # 聚合商接入容量
        self.alpha = alpha  # 聚合商不确定性alpha
        self.beta = beta  # 聚合商不确定性beta
        self.mu = self.alpha / (self.alpha + self.beta)  # 均值
        self.price_param = price_param
        self.called = 0

    def set_beta_keep_mu(self, b):
        self.beta = b
        self.alpha = self.mu / (1 - self.mu) * b

    def set_alpha_keep_mu(self, a):
        print("Warning: Please finish the function set_alpha_keep_mu")
        self.alpha = a
        # self. = self.mu / (1 - self.mu) * b

    def set_beta(self, b):
        self.beta = b
        self.update_mu()

    def set_alpha(self, a):
        self.alpha = a

    def update_mu(self):
        self.mu = self.alpha / (self.alpha + self.beta)


class ISO:
    def __init__(self, power, agg_num):
        self.p_set = power  # ISO的总调节功率，单位MW
        self.agg_num = agg_num  # AGG数量

    def set_total_power(self, p):
        self.p_set = p


def define_mpc(mpc, aggs, nid2aid={3: [0], 21: [1], 17: [2], 28: [3]}):
    kappa = mpc["kappa"]  # 负荷节点功率因数

    big_feeders = [
        [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], [3, 17],
        [25, 26, 27, 29, 30, 31, 32], [28],
        [18, 19, 20], [21],
        [22, 23, 24], []
    ]

    # 构建一个矩阵，第i行记录表示节点i以及功率流流过节点i的所有节点的id（节点编号从0开始）
    mpc["adj"] = [
        [[0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 22, 23, 24, 25, 26, 27, 29, 30, 31,
          32],
         [3, 21, 17, 28], -1],
        [[1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32],
         [3, 21, 17, 28], 0],
        [[2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32], [3, 17, 28], 1],
        [[4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 25, 26, 27, 29, 30, 31, 32], [3, 17, 28], 2],
        [[4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 25, 26, 27, 29, 30, 31, 32], [17, 28], 3],
        [[5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 25, 26, 27, 29, 30, 31, 32], [17, 28], 4],
        [[6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], [17], 5],
        [[7, 8, 9, 10, 11, 12, 13, 14, 15, 16], [17], 6],
        [[8, 9, 10, 11, 12, 13, 14, 15, 16], [17], 7],
        [[9, 10, 11, 12, 13, 14, 15, 16], [17], 8],
        [[10, 11, 12, 13, 14, 15, 16], [17], 9],
        [[11, 12, 13, 14, 15, 16], [17], 10],
        [[12, 13, 14, 15, 16], [17], 11],
        [[13, 14, 15, 16], [17], 12],
        [[14, 15, 16], [17], 13],
        [[15, 16], [17], 14],
        [[16], [17], 15],
        [[], [17], 16],
        [[18, 19, 20], [21], 17],
        [[19, 20], [21], 18],
        [[20], [21], 19],
        [[], [21], 20],
        [[22, 23, 24], [], 21],
        [[23, 24], [], 22],
        [[24], [], 23],
        [[25, 26, 27, 29, 30, 31, 32], [28], 24],
        [[26, 27, 29, 30, 31, 32], [28], 25],
        [[27, 29, 30, 31, 32], [28], 26],
        [[29, 30, 31, 32], [28], 27],
        [[29, 30, 31, 32], [], 28],
        [[30, 31, 32], [], 29],
        [[31, 32], [], 30],
        [[32], [], 31]
    ]

    # information of feeders (only 4 feeders in IEEE-33 case)
    mpc["pre_nodes"] = [[1, 2, 22, 23, 24],  # 上路
                        [1, 2, 3, 4, 5, 25, 26, 27, 28, 29, 30, 31, 32],  # 中上路
                        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],  # 中路
                        [1, 18, 19, 20, 21]  # 下路
                        ]
    mpc["pre_lines"] = [[0, 1, 21, 22, 23],  # 上路
                        [0, 1, 2, 3, 4, 24, 25, 26, 27, 28, 29, 30, 31],  # 中上路
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],  # 中路
                        [0, 17, 18, 19, 20]  # 下路
                        ]

    mpc["nid2aid"] = nid2aid  # 节点id:[agg id]

    # 流过每个节点（以及本身）的功率压降。每条记录包含2个参数，[常数P对应的压降(P*(R+kappa*X))，各个聚合商对应的压降参数(未乘c*\eta)]
    mpc["diff_v_param"] = []
    for i in range(len(mpc["adj"])):
        if i == int(mpc["root_node"]):
            mpc["diff_v_param"].append([0] * (1 + len(aggs)))
            continue
        branch_id = mpc["adj"][i][2]
        line_rx = mpc["branch"][branch_id][2] + kappa * mpc["branch"][branch_id][3]
        tmp = [0] * (1 + len(aggs))
        tmp1 = 0
        for k in range(len(mpc["adj"][i][0])):
            tmp1 += mpc["bus"][mpc["adj"][i][0][k]][2]

        for k in range(len(mpc["adj"][i][1])):
            node_aggs = mpc["nid2aid"].get(mpc["adj"][i][1][k], [])
            if len(node_aggs) != 0:
                for j in node_aggs:
                    tmp[j + 1] = line_rx
        tmp[0] = tmp1 * line_rx
        mpc["diff_v_param"].append(copy.deepcopy(tmp))

    # print("diff_v_param")
    # for i in range(len(mpc["diff_v_param"])):
    #     print(i, mpc["diff_v_param"][i])
    # print(len(mpc["diff_v_param"]), len(mpc["diff_v_param"][0]))  # 33*5

    # 每条支路终点处的电压压降
    mpc["diff_V_terminate_nodes"] = []
    for i in range(len(mpc["pre_lines"])):  # 对于每个feeder终点
        tmp = [0] * (1 + len(aggs))  # 存储格式为[常数, c_i*\eta_i前对应的系数, ... , c_i*\eta_i前对应的系数]
        for j in range(len(mpc["pre_nodes"][i])):
            for k in range(1 + len(aggs)):
                tmp[k] += mpc["diff_v_param"][mpc["pre_nodes"][i][j]][k]

        mpc["diff_V_terminate_nodes"].append(copy.deepcopy(tmp))

    # print("diff_V")
    # for i in range(len(mpc["diff_V_terminate_nodes"])):
    #     print(i, mpc["diff_V_terminate_nodes"][i])


def assignment_gurobi(iso, aggs, cost_param, mpc, method="gurobi"):
    if method == "gurobi":
        return assignment_gurobi_only(iso, aggs, cost_param, mpc)
    elif method == "cplex":
        return assignment_cplex_only(iso, aggs, cost_param, mpc)


def assignment_cplex_only(iso, aggs, cost_param, mpc):
    global V_root_sqr, V_max_sqr, V_min_sqr
    t0 = time.time()

    # Create a model
    m = cpx.Model()
    # Create variables
    etas = m.continuous_var_list([i for i in range(len(aggs))], name='X')

    # Set objective
    m.minimize(cost_param[0] +
               cost_param[1] * (iso.p_set - sum([aggs[i].capacity * E(aggs[i]) * etas[i] for i in range(len(aggs))])) +
               cost_param[2] * (
                       pow(iso.p_set - sum([aggs[i].capacity * E(aggs[i]) * etas[i] for i in range(len(aggs))]),
                           2) +
                       sum(aggs[i].capacity * aggs[i].capacity * etas[i] * etas[i] * Var(aggs[i]) for i in
                           range(len(aggs)))
               ) + sum(
        aggs[i].price_param * aggs[i].capacity * etas[i] * E(aggs[i]) for i in range(len(aggs))))

    # Add constraints:
    m.add_constraint(iso.p_set >= sum(etas[i] * aggs[i].capacity for i in range(len(aggs))), "c_e")
    for i in range(len(etas)):
        m.add_constraint(etas[i] >= 0, "c_m" + str(i + 1))
        m.add_constraint(etas[i] <= 1, "c_l" + str(i + 1))
    for i in range(4):
        m.add_constraint(
            (V_root_sqr - V_min_sqr) >= (mpc["diff_V_terminate_nodes"][i][0] + sum(
                mpc["diff_V_terminate_nodes"][i][j + 1] * etas[j] * aggs[j].capacity * E(aggs[j]) for j in
                range(len(aggs)))), "c_v" + str(i + 1))

    # solve
    t1 = time.time()
    solution = m.solve()
    t2 = time.time()

    if solution:
        eta, cost = [], solution.objective_value
        for v in [solution.get_value(v) for v in etas]:
            if abs(v) < 1e-6:
                eta.append(0)
            else:
                eta.append(v)
        # print("eta:", eta)
        # print('cost:', m.objVal)
        t3 = time.time()
        return cost, eta, t3 - t0
    else:
        print("Not optimal solution from Cplex.")
        t3 = time.time()
        return -1, [], t3 - t0


def assignment_gurobi_only(iso, aggs, cost_param, mpc):
    # 求ISO对每个AGG的最优功率分配策略
    global V_root_sqr, V_max_sqr, V_min_sqr

    t0 = time.time()

    # Create a model
    m = Model("mip1")
    m.setParam('OutputFlag', 0)

    # Create variables
    etas = []
    for i in range(len(aggs)):
        etas.append(m.addVar(vtype=GRB.CONTINUOUS, name="eta" + str(i + 1)))

    # Set objective, like E(X^2) = D(X) + E(X)^2, where X = p_set - \eta_{ri}
    m.setObjective(cost_param[0] +
                   cost_param[1] * (iso.p_set - sum(
        [aggs[i].capacity * E(aggs[i]) * etas[i] for i in
         range(len(aggs))])) +
                   cost_param[2] * (pow(iso.p_set - sum(
        [aggs[i].capacity * E(aggs[i]) * etas[i] for i in range(len(aggs))]),
                                        2) +
                                    sum(aggs[i].capacity * aggs[i].capacity * etas[i] * etas[i] * Var(aggs[i]) for i in
                                        range(len(aggs)))
                                    ) +
                   sum(aggs[i].price_param * aggs[i].capacity * etas[i] * E(aggs[i]) for i in range(len(aggs)))
                   , GRB.MINIMIZE)

    # Add inequality constraint
    m.addConstr(iso.p_set >= sum(etas[i] * aggs[i].capacity for i in range(len(aggs))), "c_e")

    # Add inequality constraints:
    for i in range(len(etas)):
        m.addConstr(etas[i] >= 0, "c_m" + str(i + 1))
        m.addConstr(etas[i] <= 1, "c_l" + str(i + 1))

    # voltage constraints
    for i in range(len(mpc["pre_lines"])):
        m.addConstr((V_root_sqr - V_min_sqr) >= (mpc["diff_V_terminate_nodes"][i][0] + sum(
            mpc["diff_V_terminate_nodes"][i][j + 1] * etas[j] * aggs[j].capacity * E(aggs[j]) for j in
            range(len(aggs)))), name="c_v" + str(i + 1))

    # solve
    t1 = time.time()
    m.optimize()
    t2 = time.time()

    if m.status == 2:
        eta, cost = [], m.objVal
        for v in m.getVars():
            if abs(v.x) < 1e-6:
                eta.append(0)
            else:
                eta.append(v.x)
        # print("eta:", eta)
        # print('cost:', m.objVal)
        t3 = time.time()
        # print(cost, end=",")

        # 将结果写入 json 文件
        m.write("optimal_solution_Gurobi_Only_" + str(len(etas)) + "_aggs_IEEE" + str(
            len(mpc["bus"])) + "_Syetem" + ".json")

        return cost, eta, t3 - t0
    else:
        print("Not optimal solution from Gurobi.")
        t3 = time.time()
        return -1, [], t3 - t0


def draw_assignment(aggs, total_power, iso, cost_param, mpc, sample_time=600, method="gurobi"):
    p_set, etas, costs = [], [], []
    for i in range(len(aggs)):
        etas.append([])

    # dead_band = deadband(aggs, cost_param)
    # print("dead-band is at p_set =", dead_band, "MW")  # 死区：最低成本策略是不需求响应的区域

    for i in range(sample_time + 1):
        iso.set_total_power(i / sample_time * total_power)
        cost, eta, _ = assignment_gurobi(iso, aggs, cost_param, mpc, method)
        p_set.append(i / sample_time * total_power)
        for j in range(len(aggs)):
            etas[j].append(eta[j])
        costs.append(cost)

    # 设置图形大小和线条颜色
    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot(111)
    lines = []
    for i in range(len(aggs)):
        line, = ax1.plot(p_set, etas[i], label=f"$LA_{{{i + 1}}}$")
        lines.append(line)
    # ax1.set_title("Demand response assignment for " + str(len(aggs)) + " aggregator.")
    ax1.set_xlabel(r"$P_{AGC}$ (kW)")
    ax1.set_ylabel(r"$\eta_{si}$")
    ax1.set_xlim([0, iso.p_set + 1])
    xticks_pos1 = [i * 100 for i in range(int(iso.p_set / 100 + 1))]
    xticks_labels1 = [f"{int(tick)}" for tick in xticks_pos1]
    plt.xticks(xticks_pos1, xticks_labels1, rotation=0)

    ax1.set_ylim([0, 1])
    yticks_pos1 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    yticks_labels1 = [f"{int(tick * 100)}\%" for tick in yticks_pos1]
    ax1.set_yticks(yticks_pos1)
    ax1.set_yticklabels(yticks_labels1)

    ax2 = ax1.twinx()
    line, = ax2.plot(p_set, costs, color='#17becf', linestyle='--', label="cost")
    lines.append(line)
    ax2.set_ylabel(r"Expectation of DSO cost (\$)   ")
    delta_y = 250
    ax2.set_ylim([0, 5 * delta_y])
    yticks_pos2 = [0, 1 * delta_y, 2 * delta_y, 3 * delta_y, 4 * delta_y, 5 * delta_y]
    yticks_labels2 = [f'{int(tick)}' for tick in yticks_pos2]
    ax2.set_yticks(yticks_pos2)
    ax2.set_yticklabels(yticks_labels2)

    # ax1.legend(handles=lines, loc='upper left', bbox_to_anchor=(1.1, 1), borderaxespad=0.)
    ax1.legend(handles=lines, loc='upper left')
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.subplots_adjust(right=0.8)

    plt.savefig("./assignment.pdf", format="pdf")
    plt.show()

    # print("Costs:", costs)
    # print("Dispatching etas:", etas)
    name = []
    for i in range(len(aggs)):
        name.append("agg" + str(i + 1))
    test = pd.DataFrame(data=etas).transpose()
    test.columns = name
    test.to_csv('./gurobi_etas.csv', encoding='gbk')
    test = pd.DataFrame(data=costs)
    test.columns = ["cost"]
    test.to_csv('./gurobi_costs.csv', encoding='gbk')


def E_cost_function(aggs, iso, cost_param, eta):
    # 已知所有参数和分配比例情况下，求E(Cost)
    delta_p = iso.p_set - sum(aggs[i].capacity * eta[i] * E(aggs[i]) for i in range(len(aggs)))
    delta_p2 = sum(aggs[i].capacity * aggs[i].capacity * eta[i] * eta[i] * Var(aggs[i]) for i in
                   range(len(aggs))) + delta_p * delta_p
    price = sum(aggs[i].price_param * aggs[i].capacity * eta[i] * E(aggs[i]) for i in range(len(aggs)))
    cost = cost_param[0] + cost_param[1] * delta_p + cost_param[2] * delta_p2 + price
    return cost


def value_estimate_function(iso, aggs, cost_param, eta, i):
    # 边际价值
    # 衡量在不同总功率分配下的第i个聚合商的价值函数，即  E(Cost)/(eta[i]) 的偏导数
    tmp_b = aggs[i].capacity * eta[i] * Var(aggs[i]) - E(aggs[i]) * iso.p_set
    for j in range(len(aggs)):
        tmp_b += aggs[j].capacity * eta[j] * E(aggs[j]) * E(aggs[i])
    cost = E(aggs[i]) * (aggs[i].price_param - cost_param[1]) + cost_param[2] * 2 * tmp_b

    # val为聚合商帮助ISO避免的损失，所以E(loss)函数对聚合商调用量的偏导是负数，意味着聚合商用得越多，ISO损失越小.聚合商的价值是正的。
    val = 0 - cost
    return val


def draw_value_estimate_function(aggs, total_power, iso, cost_param, mpc, sample_time=200,
                                 method_name="gurobi"):
    p_set, prices = [], []
    tmp_iso = copy.deepcopy(iso)
    eta = None
    for i in range(len(aggs)):
        prices.append([])

    if method_name == "TEIRP":
        positions = get_inflection_positions(tmp_iso, aggs, cost_param, mpc)
        for i in range(sample_time + 1):
            tmp_Pset = i / sample_time * total_power
            p_set.append(tmp_Pset)
            tmp_iso.set_total_power(tmp_Pset)
            tmp_Pset, tmp_eta_init, tmp_using_aggs = i / sample_time * total_power, [0 for k in range(len(aggs))], []
            stage_id = -1
            for j in range(1, len(positions)):
                if tmp_Pset <= positions[j][0]:
                    stage_id = j - 1
                    break
            # print(Pset_list[i], "MW stage_id", stage_id, tmp_Pset, positions[stage_id][0])
            etas = []
            if stage_id == len(positions) - 1:
                tmp_iso.p_set = tmp_Pset
                etas = positions[-1][1]
            else:
                tmp_iso.p_set = tmp_Pset
                tmp_eta_init, tmp_using_aggs = positions[stage_id][1], positions[stage_id][2]
                _, eta = dispatch_with_inital_state(tmp_iso, aggs, cost_param, mpc, tmp_eta_init, tmp_using_aggs)

            for j in range(len(aggs)):
                prices[j].append(-1 * value_estimate_function(tmp_iso, aggs, cost_param, eta, j))
            #     print(value_estimate_function(tmp_iso, aggs, cost_param, eta, j), end=' ')

    else:
        for i in range(sample_time + 1):
            tmp_Pset = i / sample_time * total_power
            p_set.append(tmp_Pset)
            tmp_iso.set_total_power(tmp_Pset)
            # tmp_iso.set_total_power(313*2 + i / sample_time * total_power)
            _, eta, _ = assignment_gurobi(tmp_iso, aggs, cost_param, mpc, method=method_name)
            # p_set.append(313*2 + i / sample_time * total_power)
            # print(i, tmp_iso.p_set, eta, end="; ")
            for j in range(len(aggs)):
                prices[j].append(-1 * value_estimate_function(tmp_iso, aggs, cost_param, eta, j))
            #     print(value_estimate_function(tmp_iso, aggs, cost_param, eta, j), end=' ')
            # print("")

    # 设置图形大小和线条颜色
    fig, lines = plt.figure(figsize=(6, 4)), []
    ax1 = fig.add_subplot(111)
    for i in range(len(aggs)):
        line, = ax1.plot(p_set, prices[i], label=f"$LA_{{{i + 1}}}$")
        lines.append(line)
    # ax1.set_title("Incremental rate of aggregators.")
    ax1.set_xlabel(r"$P_{AGC}$ (kW)")
    ax1.set_ylabel(r"Incremental Rate")

    y_bottom, y_top = -10, 1
    yticks_pos1 = [i for i in range(y_bottom, y_top, 1)]
    yticks_pos1.append(y_top)
    yticks_labels1 = [f"{int(tick)}" for tick in yticks_pos1]
    ax1.set_yticks(yticks_pos1)
    ax1.set_yticklabels(yticks_labels1)
    ax1.set_ylim([y_bottom, y_top])

    ax1.set_xlim([0, iso.p_set + 1])
    xticks_pos1 = [i * 100 for i in range(int(iso.p_set / 100 + 1))]
    xticks_labels1 = [f"{int(tick)}" for tick in xticks_pos1]
    plt.xticks(xticks_pos1, xticks_labels1, rotation=0)

    ax1.legend(handles=lines, loc='lower left')
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.savefig('./incremental_rate_' + method_name + '.pdf', format="pdf")
    plt.show()

    name = []
    for i in range(len(aggs)):
        name.append("LA" + str(i + 1))
    test = pd.DataFrame(data=prices).transpose()
    test.columns = name
    test.to_csv('./incremental_rate_' + method_name + '.csv', encoding='gbk')


def deadband(aggs, cost_param, return_flag=0):
    # E(cost)对每个聚合商c_n\eta_{sn}E(\eta_rn)的偏导数的最小P_{set}
    tmp_ans = []
    for i in range(len(aggs)):
        tmp_ans.append((aggs[i].price_param - cost_param[1]) / 2 / cost_param[2])
    # print("dead_zone of each aggs:", tmp_ans)
    if return_flag == 0:
        return tmp_ans
    elif return_flag == 2:
        tmp_ids, min_deadband = [], min(tmp_ans)
        for i in range(len(tmp_ans)):
            if tmp_ans[i] == min_deadband:
                tmp_ids.append(i)
        return min_deadband, tmp_ids
    else:
        return min(tmp_ans)


# ------ For one time dispatch Start------
# E(Cost) / c_i\eta_{si} 的导数相等处，是optimal dispatching solution
def get_equal_value_dispatch(iso, aggs, a, mpc, eta_init, using_aggs):
    global V_root_sqr, V_max_sqr, V_min_sqr
    not_using_aggs = []
    for i in range(len(eta_init)):
        if i not in using_aggs:
            not_using_aggs.append(i)

    formular_params_list = []  # [[a, b]]   tmp_a * \eta_{si} = tmp_b * V + tmp_const
    i_0 = using_aggs[0]
    for i in range(len(using_aggs)):
        tmp_i = using_aggs[i]  # tmp_i is agg id.
        if i == 0:
            # tmp_a * \eta_{s0} = tmp_b * r + tmp_const, 此处是指using_aggs中的第1个agg_id
            tmp_a = 2.0 * a[2] * aggs[tmp_i].capacity * Var(aggs[tmp_i]) \
                    + 2.0 * a[2] * aggs[tmp_i].capacity * E(aggs[tmp_i]) * E(aggs[tmp_i]) \
                    + 2.0 * a[2] * E(aggs[tmp_i]) \
                    * sum((aggs[k].capacity * E(aggs[k])) * ((aggs[tmp_i].capacity * Var(aggs[tmp_i]) * E(aggs[k])) / (
                    E(aggs[tmp_i]) * aggs[k].capacity * Var(aggs[k]))) for k in using_aggs[1:])

            tmp_const = E(aggs[tmp_i]) * (a[1] - aggs[tmp_i].price_param) \
                        + 2.0 * a[2] * E(aggs[tmp_i]) * iso.p_set \
                        - 2.0 * a[2] * E(aggs[tmp_i]) * \
                        (sum(aggs[k].capacity * E(aggs[k]) *
                             E(aggs[k]) * (aggs[tmp_i].price_param - aggs[k].price_param) /
                             (2 * a[2] * aggs[k].capacity * Var(aggs[k])) for k in
                             using_aggs[1:])
                         + sum(aggs[k].capacity * eta_init[k] * E(aggs[k]) for k in not_using_aggs))

            tmp_b = 1.0 - 2.0 * a[2] * E(aggs[tmp_i]) \
                    * sum((aggs[k].capacity * E(aggs[k])) * (
                    (E(aggs[tmp_i]) - E(aggs[k])) / (
                    2.0 * a[2] * aggs[k].capacity * Var(aggs[k]) * E(aggs[tmp_i]))) for k in using_aggs[1:])

            formular_params_list.append([tmp_b / tmp_a, tmp_const / tmp_a])  # tmp_a* \eta_{s0} = tmp_b * r + tmp_const
        else:
            # tmp_a(=1) * \eta_{si} = tmp_b * r + tmp_const
            eta_coefficient = (aggs[i_0].capacity * Var(aggs[i_0]) * E(aggs[tmp_i])) \
                              / (E(aggs[i_0]) * aggs[tmp_i].capacity * Var(aggs[tmp_i]))

            tmp_b = (E(aggs[i_0]) - E(aggs[tmp_i])) \
                    / (2.0 * a[2] * aggs[tmp_i].capacity * Var(aggs[tmp_i]) * E(aggs[i_0])) \
                    + eta_coefficient * formular_params_list[0][0]

            tmp_const = (E(aggs[tmp_i]) * (aggs[i_0].price_param - aggs[tmp_i].price_param)) \
                        / (2.0 * a[2] * aggs[tmp_i].capacity * Var(aggs[tmp_i])) \
                        + eta_coefficient * formular_params_list[0][1]

            formular_params_list.append([tmp_b, tmp_const])  # tmp_a(=1) * \eta_{si} = tmp_b * r + tmp_const
            # could draw function of \eta_{si}~r by https://www.desmos.com/calculator?lang=zh-CN

    tmp_a = sum(formular_params_list[i][0] * aggs[using_aggs[i]].capacity for i in range(len(using_aggs)))
    tmp_b = sum(formular_params_list[i][1] * aggs[using_aggs[i]].capacity for i in range(len(using_aggs)))
    tmp_b += sum(eta_init[k] * aggs[k].capacity for k in not_using_aggs)
    formular_params_list.append([tmp_a / iso.p_set, tmp_b / iso.p_set])
    # print("s.t. formulars params list (0 <= a * \eta_{si} + b <= 1):\n", formular_params_list)

    # tmp_const = (V_root_sqr - V_min_sqr) * (10 * 10 * 1e6) / 1e3  # 参考assignment_gurobi函数中的电压约束的常数部分。
    tmp_const = V_root_sqr - V_min_sqr  # 参考assignment_gurobi函数中的电压约束的常数部分。
    for i in range(len(mpc["diff_V_terminate_nodes"])):  # 4为feeder数量
        tmp_a, tmp_b, tmp_c = mpc["diff_V_terminate_nodes"][i][0], 0, 0
        for j in not_using_aggs:
            tmp_a += mpc["diff_V_terminate_nodes"][i][j + 1] * aggs[j].capacity * eta_init[j] * E(aggs[j])
        for j in range(len(using_aggs)):
            tmp_agg_id = using_aggs[j]
            tmp_1 = mpc["diff_V_terminate_nodes"][i][tmp_agg_id + 1] * aggs[tmp_agg_id].capacity * E(aggs[tmp_agg_id])
            tmp_b += tmp_1 * formular_params_list[j][0]
            tmp_c += tmp_1 * formular_params_list[j][1]
            # (V_root_sqr - V_min_sqr) * (10 * 10 * 1e6) >= (1e3) * (mpc["diff_V_terminate_nodes"][i][0] + sum(
            #     mpc["diff_V_terminate_nodes"][i][j + 1] * etas[j] * aggs[j].capacity * E(aggs[i])
        formular_params_list.append([tmp_b / tmp_const, (tmp_a + tmp_c) / tmp_const])

    # print("formular_params_list")
    # for x in formular_params_list:
    #     print(x)

    # 根据\eta_{si}约束，使用线性规划求解r
    manual_LP, r = True, 0
    if manual_LP:
        # 手动求解LP
        tmp_v = 9999999999
        for i in range(len(formular_params_list)):
            if i < len(using_aggs):
                tmp_up_bound = (1 - formular_params_list[i][1]) / formular_params_list[i][0]
                tmp_low_bound = (0 - formular_params_list[i][1]) / formular_params_list[i][0]
                if formular_params_list[i][0] < 0:
                    tmp_exchange = tmp_up_bound
                    tmp_up_bound = tmp_low_bound
                    tmp_low_bound = tmp_exchange
                # print("tmp_r_range :", tmp_low_bound, "<= r <=", tmp_up_bound)
            else:
                if formular_params_list[i][0] < 0:
                    continue
                tmp_up_bound = (1 - formular_params_list[i][1]) / (formular_params_list[i][0])
                # print("tmp_r_range :", "r <=", tmp_up_bound)
            if tmp_v >= tmp_up_bound:
                tmp_v = tmp_up_bound
        # print("tmp_r = ", tmp_v)
        r = min(tmp_v, 0)  # 此处的变量r为边际价值,即-r.限制边际价值(-r)的下限,进而避免\eta超过1.

    # 如果多聚合商等价的话，计算协同使用的比例
    Equivalent_distribution_ratios = []
    # print("r=", r)
    return r, Equivalent_distribution_ratios


def dispatch_with_inital_state(iso, aggs, cost_param, mpc, eta_init, using_aggs):
    N, etas = len(aggs), []

    # 挑选需要改变的几个agg求解
    # 求边际价格k
    pi, Equivalent_distribution_ratios = 0, []
    if len(using_aggs) >= 1:
        pi, Equivalent_distribution_ratios = get_equal_value_dispatch(iso, aggs, cost_param, mpc, eta_init, using_aggs)
    # print("Incremental rate for assignment:", pi)
    # 求边际价格为k时的分配比例
    A1 = np.zeros((len(using_aggs), len(using_aggs)))
    B1 = np.zeros((len(using_aggs), 1))
    for i in range(len(using_aggs)):
        tmp_id = using_aggs[i]

        B1[i][0] = pi * aggs[tmp_id].capacity \
                   + aggs[tmp_id].capacity * E(aggs[tmp_id]) * cost_param[1] \
                   + cost_param[2] * 2 * aggs[tmp_id].capacity * E(aggs[tmp_id]) * iso.p_set \
                   - aggs[tmp_id].price_param * aggs[tmp_id].capacity * E(aggs[tmp_id])
        for j in range(N):
            if j not in using_aggs:
                B1[i][0] -= cost_param[2] * 2 * aggs[tmp_id].capacity * E(aggs[tmp_id]) * aggs[j].capacity * eta_init[
                    j] * E(aggs[j])

        A1[i][i] = 2 * cost_param[2] * aggs[tmp_id].capacity * aggs[tmp_id].capacity * Var(aggs[tmp_id])
        for k in range(len(using_aggs)):
            j = using_aggs[k]
            A1[i][k] += cost_param[2] * 2.0 * aggs[j].capacity * E(aggs[j]) * aggs[tmp_id].capacity * E(aggs[tmp_id])

    x = np.linalg.pinv(A1) @ B1
    # print(A1, "\n", B1)
    # print("x:", x)
    # print(A1.dot(x) - B1)

    # save x -> etas
    tmp_i = 0
    for i in range(N):
        if i in using_aggs:
            # etas.append(min(x[tmp_i][0], changed_Pset / aggs[i].capacity))
            etas.append(x[tmp_i][0])
            tmp_i += 1
        else:
            etas.append(eta_init[i])
    cost = E_cost_function(aggs, iso, cost_param, etas)
    # print("etas:", etas)
    # print("Cost:", E_cost_function(aggs, iso, cost_param, etas))
    return cost, etas


# ------ For one time dispatch End------


# ------ For get_next_inflection Start------
def get_equal_value_inflection(aggs, a, mpc, eta_init, using_aggs, agg_j):
    # print("eta_init", eta_init)
    # print("using_aggs", using_aggs)
    # print("agg_j", agg_j)
    global V_root_sqr, V_max_sqr, V_min_sqr

    # 判断agg_j替代using_aggs中的某一个时的Pset
    not_using_aggs = []
    for i in range(len(eta_init)):
        if i not in using_aggs:
            not_using_aggs.append(i)

    formular_params_list = []  # [[a, b, c]]   \eta_{si} = a * r + b * Pset + c
    i_0 = using_aggs[0]
    using_aggs.append(agg_j)
    for i in range(len(using_aggs)):
        tmp_i = using_aggs[i]  # tmp_i is agg id.
        # 0 <= \eta_{s0} <= 1
        if i == 0:
            # tmp_a * \eta_{s0} = tmp_b * r + tmp_const, 此处是指using_aggs中的第1个agg_id, r为price
            tmp_a = 2.0 * a[2] * aggs[tmp_i].capacity * Var(aggs[tmp_i]) \
                    + 2.0 * a[2] * aggs[tmp_i].capacity * E(aggs[tmp_i]) * E(aggs[tmp_i]) \
                    + 2.0 * a[2] * E(aggs[tmp_i]) \
                    * sum((aggs[k].capacity * E(aggs[k])) * ((aggs[tmp_i].capacity * Var(aggs[tmp_i]) * E(aggs[k])) / (
                    E(aggs[tmp_i]) * aggs[k].capacity * Var(aggs[k]))) for k in using_aggs[1:])

            tmp_const = E(aggs[tmp_i]) * (a[1] - aggs[tmp_i].price_param) \
                        - 2.0 * a[2] * E(aggs[tmp_i]) * \
                        (sum(aggs[k].capacity * E(aggs[k]) *
                             E(aggs[k]) * (aggs[tmp_i].price_param - aggs[k].price_param) /
                             (2 * a[2] * aggs[k].capacity * Var(aggs[k])) for k in
                             using_aggs[1:])
                         + sum(aggs[k].capacity * eta_init[k] * E(aggs[k]) for k in not_using_aggs))

            tmp_b = 1.0 - 2.0 * a[2] * E(aggs[tmp_i]) \
                    * sum((aggs[k].capacity * E(aggs[k])) * (
                    (E(aggs[tmp_i]) - E(aggs[k])) / (
                    2.0 * a[2] * aggs[k].capacity * Var(aggs[k]) * E(aggs[tmp_i]))) for k in using_aggs[1:])

            tmp_c = 2.0 * a[2] * E(aggs[tmp_i])

            # tmp_a * \eta_{si} = tmp_b * r + tmp_c * iso.p_set + tmp_const
            formular_params_list.append([tmp_b / tmp_a, tmp_c / tmp_a, tmp_const / tmp_a])
        else:
            # 0 <= \eta_{si} <= 1
            eta_coefficient = (aggs[i_0].capacity * Var(aggs[i_0]) * E(aggs[tmp_i])) \
                              / (E(aggs[i_0]) * aggs[tmp_i].capacity * Var(aggs[tmp_i]))

            tmp_b = (E(aggs[i_0]) - E(aggs[tmp_i])) \
                    / (2.0 * a[2] * aggs[tmp_i].capacity * Var(aggs[tmp_i]) * E(aggs[i_0])) \
                    + eta_coefficient * formular_params_list[0][0]

            tmp_c = eta_coefficient * formular_params_list[0][1]

            tmp_const = (E(aggs[tmp_i]) * (aggs[i_0].price_param - aggs[tmp_i].price_param)) / (
                    2.0 * a[2] * aggs[tmp_i].capacity * Var(aggs[tmp_i])) \
                        + eta_coefficient * formular_params_list[0][2]

            # tmp_a(=1) * \eta_{si} = tmp_b * r + tmp_c * iso.p_set + tmp_const
            formular_params_list.append([tmp_b, tmp_c, tmp_const])
            # could draw function of \eta_{si}~r by https://www.desmos.com/calculator?lang=zh-CN

    # V = v1 * Pset + v2 根据agg_j的偏导=V，找到V和Pset的关系.进而代入不等式组，求最大Pset
    tmp_formular_params_list_agg_j = copy.deepcopy(formular_params_list[-1])
    v1 = -tmp_formular_params_list_agg_j[1] / (tmp_formular_params_list_agg_j[0])
    v2 = (eta_init[agg_j] - tmp_formular_params_list_agg_j[2]) / tmp_formular_params_list_agg_j[0]
    del formular_params_list[-1]
    del using_aggs[-1]
    ans_eta_V_Pset_const = copy.deepcopy(formular_params_list)

    new_formular_params_list = []
    for i in range(len(formular_params_list)):
        tmp11 = formular_params_list[i][1] + formular_params_list[i][0] * v1
        # if abs(tmp11 - 0) < 1e-6: tmp11 = 1e-6
        new_formular_params_list.append([tmp11,
                                         formular_params_list[i][2] + formular_params_list[i][0] * v2])

    # print("new_formular_params_list:", new_formular_params_list)
    # sum(c_i*\eta_{si}) <= P_set  ==> tmp_a * P_set + tmp_b + 1 <= 1
    tmp_a = sum(new_formular_params_list[i][0] * aggs[using_aggs[i]].capacity for i in range(len(using_aggs))) - 1
    tmp_b = sum(new_formular_params_list[i][1] * aggs[using_aggs[i]].capacity for i in range(len(using_aggs)))
    tmp_b += sum(eta_init[k] * aggs[k].capacity for k in not_using_aggs) + 1
    new_formular_params_list.append([tmp_a, tmp_b])
    # 电压约束
    # sum(k_i*c_i*\eta_{si}*E_i) + k_0 <= const. 其中k_i为mpc中计算的常数参数,const为电压下限
    # ==> tmp_a * P_set + tmp_b <= 1
    # tmp_const = (V_root_sqr - V_min_sqr) * (10 * 10 * 1e6) / 1e3  # 参考assignment_gurobi函数中的电压约束的常数部分。
    tmp_const = V_root_sqr - V_min_sqr  # 参考assignment_gurobi函数中的电压约束的常数部分。
    for i in range(len(mpc["diff_V_terminate_nodes"])):  # 4为feeder数量
        tmp_a, tmp_b, tmp_c = mpc["diff_V_terminate_nodes"][i][0], 0, 0
        for j in not_using_aggs:
            tmp_a += mpc["diff_V_terminate_nodes"][i][j + 1] * aggs[j].capacity * eta_init[j] * E(aggs[j])
        for j in range(len(using_aggs)):
            tmp_agg_id = using_aggs[j]
            tmp_1 = mpc["diff_V_terminate_nodes"][i][tmp_agg_id + 1] * aggs[tmp_agg_id].capacity * E(aggs[tmp_agg_id])
            tmp_b += tmp_1 * new_formular_params_list[j][0]
            tmp_c += tmp_1 * new_formular_params_list[j][1]

            # (V_root_sqr - V_min_sqr) * (10 * 10 * 1e6) >= (1e3) * (mpc["diff_V_terminate_nodes"][i][0] + sum(
            #     mpc["diff_V_terminate_nodes"][i][j + 1] * etas[j] * aggs[j].capacity * E(aggs[i])
        new_formular_params_list.append([tmp_b / tmp_const, (tmp_a + tmp_c) / tmp_const])

    formular_params_list = new_formular_params_list
    # print("formular_params_list:", formular_params_list)

    ans_Pset, ans_eta = 0, []  # 下一拐点处 Pset、eta
    # 根据\eta_{si}约束，使用线性规划求解V
    manual_LP = True
    if manual_LP:
        # 手动求解LP
        tmp_pset_up = 9999999999  # tmp_Pset
        tmp_pset_low = -999999999
        tmp_low_bound, last_tmp_low_bound = 0, 0
        for i in range(len(formular_params_list)):
            if i < len(using_aggs):
                tmp_up_bound = (1 - formular_params_list[i][1]) / (
                    formular_params_list[i][0] if formular_params_list[i][0] != 0 else 1e-8)
                tmp_low_bound = (0 - formular_params_list[i][1]) / (
                    formular_params_list[i][0] if formular_params_list[i][0] != 0 else 1e-8)
                if formular_params_list[i][0] < 0:
                    tmp_exchange = tmp_up_bound
                    tmp_up_bound = tmp_low_bound
                    tmp_low_bound = tmp_exchange
                # print("tmp_v_range :", tmp_low_bound, "<= v <=", tmp_up_bound)
                if tmp_pset_up >= tmp_up_bound:
                    tmp_pset_up = tmp_up_bound
                if tmp_pset_low <= tmp_low_bound:
                    if i != len(formular_params_list) - 1:
                        tmp_pset_low = tmp_low_bound
                # print("tmp_v_range :", tmp_low_bound, "<= v <=", tmp_up_bound, tmp_pset_low, tmp_pset_up)
                # print("tmp_pset_low =", tmp_pset_low, "tmp_pset_up =", tmp_pset_up)
            else:
                tmp_up_bound = (1 - formular_params_list[i][1]) / (
                    formular_params_list[i][0] if formular_params_list[i][0] != 0 else 1e-8)
                if tmp_pset_up >= tmp_up_bound:
                    tmp_pset_up = tmp_up_bound
                # print("tmp_pset_up =", tmp_pset_up)

        ans_Pset = tmp_pset_up  # ans_Pset = 9999999999
        if tmp_pset_low > tmp_pset_up or tmp_pset_up < 0:
            ans_Pset = 9999999999  # 如果无解，那么下一个拐点是using_aggs中最先用光容量的那个Pset
    return ans_Pset, ans_eta_V_Pset_const


# 求解using_agg中有多个聚合商时，由于聚合商之间不能维持边际价值相等的状态导致的拐点
def get_equal_aggs_inflection_end(aggs, a, mpc, eta_init, using_aggs, agg_j):
    # 获取多个aggs等价的终点处的Pset
    # agg_j就是using_aggs中触碰上边界eta_{si}=1的那个聚合商,传入的using_aggs已经将其移除了

    global V_root_sqr, V_max_sqr, V_min_sqr
    # print("get_equal_aggs_inflection_end", agg_j)
    # print("\nusing_aggs", using_aggs, agg_j)
    # print("eta_init", eta_init)

    not_using_aggs = []
    for i in range(len(eta_init)):
        if i not in using_aggs:
            not_using_aggs.append(i)
    last_pset = sum(eta_init[i] * aggs[i].capacity for i in range(len(aggs)))
    # print("not_using_aggs:", not_using_aggs)

    formular_params_list = []  # [[a, b, c]]   \eta_{si} = a * V + b * Pset + c
    i_0 = using_aggs[0]
    for i in range(len(using_aggs) + 1):
        tmp_i = -1  # tmp_i is agg id.
        if i == len(using_aggs):  # 在将eta_{s agg_j}=1代入其他偏导数等式的情况下，仍然获取agg_j关于r与Pset的等式关系式
            tmp_i = agg_j
        else:
            tmp_i = using_aggs[i]

        # 0 <= \eta_{s0} <= 1
        if i == 0:
            # tmp_a * \eta_{si} = tmp_b * r + tmp_c * iso.p_set + tmp_const, 此处是指using_aggs中的第1个agg_id, r为price
            # 此处公式是按照Incremental rate推导的, (Incremental rate = -1.0 * \frac{\partial E(Cost)}{\partial c_i\eta_{si}})
            tmp_a = 2.0 * a[2] * aggs[tmp_i].capacity * Var(aggs[tmp_i]) \
                    + 2.0 * a[2] * aggs[tmp_i].capacity * E(aggs[tmp_i]) * E(aggs[tmp_i]) \
                    + 2.0 * a[2] * E(aggs[tmp_i]) \
                    * sum((aggs[k].capacity * E(aggs[k])) * ((aggs[tmp_i].capacity * Var(aggs[tmp_i]) * E(aggs[k])) / (
                    E(aggs[tmp_i]) * aggs[k].capacity * Var(aggs[k]))) for k in using_aggs[1:])

            tmp_const = E(aggs[tmp_i]) * (a[1] - aggs[tmp_i].price_param) \
                        - 2.0 * a[2] * E(aggs[tmp_i]) * \
                        (sum(aggs[k].capacity * E(aggs[k]) *
                             E(aggs[k]) * (aggs[tmp_i].price_param - aggs[k].price_param) /
                             (2 * a[2] * aggs[k].capacity * Var(aggs[k])) for k in
                             using_aggs[1:])
                         + sum(aggs[k].capacity * eta_init[k] * E(aggs[k]) for k in not_using_aggs))

            tmp_b = 1.0 - 2.0 * a[2] * E(aggs[tmp_i]) \
                    * sum((aggs[k].capacity * E(aggs[k])) * (
                    (E(aggs[tmp_i]) - E(aggs[k])) / (
                    2.0 * a[2] * aggs[k].capacity * Var(aggs[k]) * E(aggs[tmp_i]))) for k in using_aggs[1:])

            tmp_c = 2.0 * a[2] * E(aggs[tmp_i])

            # tmp_a * \eta_{si} = tmp_b * r + tmp_c * iso.p_set + tmp_const
            formular_params_list.append([tmp_b / tmp_a, tmp_c / tmp_a, tmp_const / tmp_a])
        else:
            # 0 <= \eta_{si} <= 1
            eta_coefficient = (aggs[i_0].capacity * Var(aggs[i_0]) * E(aggs[tmp_i])) \
                              / (E(aggs[i_0]) * aggs[tmp_i].capacity * Var(aggs[tmp_i]))

            tmp_b = (E(aggs[i_0]) - E(aggs[tmp_i])) \
                    / (2.0 * a[2] * aggs[tmp_i].capacity * Var(aggs[tmp_i]) * E(aggs[i_0])) \
                    + eta_coefficient * formular_params_list[0][0]

            tmp_c = eta_coefficient * formular_params_list[0][1]

            tmp_const = (E(aggs[tmp_i]) * (aggs[i_0].price_param - aggs[tmp_i].price_param)) / (
                    2.0 * a[2] * aggs[tmp_i].capacity * Var(aggs[tmp_i])) \
                        + eta_coefficient * formular_params_list[0][2]

            # tmp_a(=1) * \eta_{si} = tmp_b * r + tmp_c * iso.p_set + tmp_const
            formular_params_list.append([tmp_b, tmp_c, tmp_const])
            # could draw function of \eta_{si}~r by https://www.desmos.com/calculator?lang=zh-CN

    # print("formular_params_list(V,Pset,Const)", formular_params_list)
    # 根据agg_j的偏导=V且eta_{s agg_j}=1,找到V和Pset的关系，形如V = v1 * Pset + v2. 进而代入不等式组，求最大Pset(此处V也就是r)
    tmp_formular_params_list_agg_j = copy.deepcopy(formular_params_list[-1])
    v1 = -tmp_formular_params_list_agg_j[1] / tmp_formular_params_list_agg_j[0]
    v2 = (eta_init[agg_j] - tmp_formular_params_list_agg_j[2]) / tmp_formular_params_list_agg_j[0]
    # print("v1, v2", v1, v2)
    del formular_params_list[-1]
    ans_eta_V_Pset_const = copy.deepcopy(formular_params_list)

    new_formular_params_list = []
    for i in range(len(formular_params_list)):
        new_formular_params_list.append([formular_params_list[i][1] + formular_params_list[i][0] * v1,
                                         formular_params_list[i][2] + formular_params_list[i][0] * v2])

    # print("new_formular_params_list:", new_formular_params_list)
    # sum(c_i*\eta_{si}) <= P_set  ==> tmp_a * P_set + tmp_b + 1 <= 1
    tmp_a = sum(new_formular_params_list[i][0] * aggs[using_aggs[i]].capacity for i in range(len(using_aggs))) - 1
    tmp_b = sum(new_formular_params_list[i][1] * aggs[using_aggs[i]].capacity for i in range(len(using_aggs)))
    tmp_b += sum(eta_init[k] * aggs[k].capacity for k in not_using_aggs) + 1
    new_formular_params_list.append([tmp_a, tmp_b])

    # 电压约束
    # sum(k_i*c_i*\eta_{si}*E_i) + k_0 <= const. 其中k_i为mpc中计算的常数参数,const为电压下限
    # ==> tmp_a * P_set + tmp_b <= 1
    # tmp_const = (V_root_sqr - V_min_sqr) * (10 * 10 * 1e6) / 1e3  # 参考assignment_gurobi函数中的电压约束的常数部分。
    tmp_const = V_root_sqr - V_min_sqr  # 参考assignment_gurobi函数中的电压约束的常数部分。
    for i in range(len(mpc["diff_V_terminate_nodes"])):
        # (V_root_sqr - V_min_sqr) * (10 * 10 * 1e6) >= (1e3) * (mpc["diff_V_terminate_nodes"][i][0] + sum(
        #     mpc["diff_V_terminate_nodes"][i][j + 1] * etas[j] * aggs[j].capacity * E(aggs[i])
        tmp_a, tmp_b, tmp_c = mpc["diff_V_terminate_nodes"][i][0], 0, 0
        for j in not_using_aggs:
            tmp_a += mpc["diff_V_terminate_nodes"][i][j + 1] * aggs[j].capacity * eta_init[j] * E(aggs[j])
        for j in range(len(using_aggs)):
            tmp_agg_id = using_aggs[j]
            tmp_1 = mpc["diff_V_terminate_nodes"][i][tmp_agg_id + 1] * aggs[tmp_agg_id].capacity * E(aggs[tmp_agg_id])
            tmp_b += tmp_1 * new_formular_params_list[j][0]
            tmp_c += tmp_1 * new_formular_params_list[j][1]
        new_formular_params_list.append([tmp_b / tmp_const, (tmp_a + tmp_c) / tmp_const])

    formular_params_list = new_formular_params_list
    # print("formular_params_list:(Pset,Const)", formular_params_list)

    ans_Pset, ans_eta = 0, []  # 下一拐点处 Pset、eta
    # 根据\eta_{si}约束，使用线性规划求解V
    manual_LP = True
    if manual_LP:
        # 手动求解LP
        tmp_pset_up = 9999999999  # tmp_Pset
        tmp_pset_low = -999999999

        for i in range(len(formular_params_list)):
            if i < len(using_aggs):
                tmp_up_bound = (1 - formular_params_list[i][1]) / formular_params_list[i][0]
                tmp_low_bound = (0 - formular_params_list[i][1]) / formular_params_list[i][0]
                if formular_params_list[i][0] < 0:
                    tmp_exchange = tmp_up_bound
                    tmp_up_bound = tmp_low_bound
                    tmp_low_bound = tmp_exchange

                if tmp_pset_up >= tmp_up_bound:
                    tmp_pset_up = tmp_up_bound
                if tmp_pset_low <= tmp_low_bound:
                    if i != len(formular_params_list) - 1:
                        tmp_pset_low = tmp_low_bound
                # print("tmp_v_range :", tmp_low_bound, "<= v <=", tmp_up_bound, tmp_pset_low, tmp_pset_up)
                # print("tmp_pset_low =", tmp_pset_low, "tmp_pset_up =", tmp_pset_up)
            else:
                tmp_up_bound = (1 - formular_params_list[i][1]) / (
                    formular_params_list[i][0] if formular_params_list[i][0] != 0 else 1e-8)
                if tmp_pset_up >= tmp_up_bound:
                    tmp_pset_up = tmp_up_bound
                # print("tmp_pset_up =", tmp_pset_up)

        # print("tmp_pset_low =", tmp_pset_low, "tmp_pset_up =", tmp_pset_up)
        ans_Pset = tmp_pset_up  # 999999999
        if tmp_pset_low > tmp_pset_up or tmp_pset_up < 0:
            ans_Pset = 9999999999  # 如果无解，那么下一个拐点是using_aggs中最先用光容量的那个Pset
        else:
            if tmp_pset_low > last_pset:
                ans_Pset = tmp_pset_low
        # print("tmp_Pset =", ans_Pset)
    # print("\n")
    return ans_Pset, ans_eta_V_Pset_const


# 求解有最优微增量的聚合商被完全调用不得不退出ALA时的Pset
def get_Pset_of_Euqal_agg_end_point(iso, aggs, cost_param, mpc, eta_init, using_aggs):
    # print("input of get_Pset_of_Euqal_agg_end_point:", iso.p_set, eta_init, using_aggs)
    # 200.9999999999174 [0, np.float64(0.9999999999991739), 0, np.float64(0.9999999999999998), 0] [4]

    not_using_aggs, full_called_aggs = [], []
    for i in range(len(eta_init)):
        if i not in using_aggs:
            not_using_aggs.append(i)
        if abs(1 - i) < 1e-4:
            full_called_aggs.append(i)
    # print("----------------------")
    # print("using_aggs", using_aggs)
    # print("eta_init", eta_init, iso.p_set)

    if len(using_aggs) == 1:
        # 获取单个等价agg达到某一个agg的容量上限处的Pset
        # print("\n------finding Pset of single-aggs end------")

        # 认为 \Delta Pset=\Delta \eta_{s0} * c_0 * E(aggs[0]). Pset的增量=LA0的调用增量的期望
        next_Pset = iso.p_set + (1 - eta_init[using_aggs[0]]) * aggs[using_aggs[0]].capacity * max(1,
                                                                                                   E(aggs[using_aggs[
                                                                                                       0]]))
        # print("DEBUG next_Pset:", next_Pset)
        # print("----------------------")
        return next_Pset
    elif len(using_aggs) >= 2:
        # 获取using_agg内部某个成员无法维持最高的边际价值导致拐点的Pset.
        # print("\n------finding Pset of multi-aggs end------")
        next_Psets = []
        for i in range(len(using_aggs)):
            tmp_eta_init = copy.deepcopy(eta_init)
            tmp_eta_init[using_aggs[i]] = 1.0
            tmp_using_aggs = [x for x in using_aggs if x != using_aggs[i]]
            next_Pset, b = get_equal_aggs_inflection_end(aggs, cost_param, mpc, tmp_eta_init, tmp_using_aggs,
                                                         using_aggs[i])
            tmp_iso1 = iso
            tmp_iso1.p_set = next_Pset
            next_Psets.append(next_Pset)
        # print("DEBUG next_Psets:", next_Psets, "len(using_aggs) >= 1")
        # print("------------\n")
        return min(next_Psets)
    else:
        return 999999999


# 求解当前ALA为空集时，NLA中某聚合商的微增量刚刚<0时的解（由于ALAS为空，微增量<0就是最优微增量）
def get_next_value_eq_0(aggs, a, etas, i):
    # 当前eta_init下使得using_aggs为空，求下一次using_aggs不为空时，P_set的值
    a1 = -aggs[i].capacity * E(aggs[i]) * a[1] + 2 * a[2] * aggs[i].capacity * aggs[i].capacity * etas[i] * Var(
        aggs[i]) + 2 * a[2] * aggs[i].capacity * E(aggs[i]) * sum(
        aggs[k].capacity * etas[k] * E(aggs[k]) for k in range(len(aggs))) + aggs[i].price_param * aggs[i].capacity * E(
        aggs[i])
    b1 = 2 * a[2] * aggs[i].capacity * E(aggs[i])
    return a1 / b1


def get_inflection_positions(iso, aggs, cost_param, mpc):
    # Calculate Pset if existing aggs[j] could replace any one of using_aggs
    tmp_iso = copy.deepcopy(iso)
    positions = [[0, [0 for i in range(len(aggs))], []]]
    next_Pset, next_eta_init, next_using_aggs = 0, [], []

    # deadband round
    next_Pset, next_using_aggs = deadband(aggs, cost_param, return_flag=2)
    next_eta_init = [0 for i in range(len(aggs))]
    positions.append([next_Pset, next_eta_init, next_using_aggs])
    tmp_iso.p_set = next_Pset

    # dispatch rounds
    first_full_called_flag = True
    # while tmp_iso.p_set<=iso.p_set:
    # set current state
    not_using_aggs = []
    for j in range(len(next_eta_init)):
        if j not in next_using_aggs:
            not_using_aggs.append(j)
    # print("init state:\nPset", next_Pset)
    # print("eta_init", next_eta_init)
    # print("using_aggs", next_using_aggs)
    # print("not_using_aggs", not_using_aggs)
    deadbands = deadband(aggs, cost_param, return_flag=1)
    empty_step_flag = False
    first_after_skip_area = True
    last_skip_Pset_range = positions[-1][0] - positions[-2][0]

    # next_Pset, eta_V_Pset_const = get_equal_value_inflection(aggs, cost_param, mpc, next_eta_init, next_using_aggs, 2)
    # print("------next_Pset", next_Pset, deadbands[2])
    wjc_it = 0
    while next_Pset < iso.p_set:
        # if wjc_it > 20: break
        wjc_it += 1
        # 如果这一步没有using_aggs，那么直接在本轮中计算并存储下次有变化的critical point.
        last_pset = next_Pset
        if empty_step_flag:
            empty_step_flag = False
            skip_to_Pset = []
            for i in range(len(aggs)):
                if abs(next_eta_init[i] - 1) < 1e-6:
                    skip_to_Pset.append(99999999)
                else:
                    skip_to_Pset.append(get_next_value_eq_0(aggs, cost_param, positions[-1][1], i))
            next_Pset, next_using_aggs = min(skip_to_Pset), []
            for i in range(len(aggs)):
                if abs(next_Pset - skip_to_Pset[i]) < 1e-6:
                    next_using_aggs.append(i)

            tmp_iso.p_set = next_Pset
            positions.append([next_Pset, copy.deepcopy(positions[-1][1]), next_using_aggs])
            first_after_skip_area = True
            last_skip_Pset_range = positions[-1][0] - positions[-2][0]

        else:
            # 当前using_aggs不为空
            # print("\n---iter", wjc_it, "current_Pset", next_Pset, next_using_aggs, "--------")
            # set current state
            not_using_aggs, full_used_aggs = [], []
            for j in range(len(next_eta_init)):
                if j not in next_using_aggs:
                    not_using_aggs.append(j)
                if abs(next_eta_init[j] - 1) < 1e-5:
                    full_used_aggs.append(j)
            current_using_aggs = copy.deepcopy(next_using_aggs)
            # current_Pset = next_Pset
            # print("not_using_aggs:", not_using_aggs)
            # update Pset
            next_Psets = []
            for i in range(len(aggs)):
                if i in not_using_aggs:
                    # print(i, next_eta_init, next_using_aggs)
                    next_Pset, _ = get_equal_value_inflection(aggs, cost_param, mpc, next_eta_init, next_using_aggs, i)
                    # print("tmp_next_Pset1", next_Pset)
                    # verify Value of aggs on the end position to avoid numbercial computation accuracy
                    if next_Pset < 999999:
                        tmp_iso1 = copy.deepcopy(iso)
                        tmp_iso1.p_set = next_Pset
                        _, tmp_next_eta_init = dispatch_with_inital_state(tmp_iso1, aggs, cost_param, mpc,
                                                                          next_eta_init,
                                                                          next_using_aggs)
                        # print("tmp_next_eta_init", tmp_next_eta_init)
                        tmp_next_V = []
                        for j in range(len(aggs)):
                            if abs(1 - tmp_next_eta_init[j]) < 1e-6:
                                tmp_next_V.append(-1)
                            else:
                                tmp_next_V.append(
                                    value_estimate_function(tmp_iso1, aggs, cost_param, tmp_next_eta_init, j))
                        # print("tmp_next_V", tmp_next_V, next_using_aggs[0], i)
                        # if abs(tmp_next_V[next_using_aggs[0]] - tmp_next_V[i]) > 1e-3:
                        #     next_Psets.append(9999999999)
                        if tmp_next_V[i] < 1e-3:  # 和上面2行if注释选用一个
                            next_Psets.append(9999999999)
                        else:
                            next_Psets.append(next_Pset)
                    else:
                        next_Psets.append(8888888888)
                else:
                    next_Psets.append(7777777777)
            # print("next_Psets(0):", next_Psets)
            # print("Euqal_agg_end_point:")
            next_Psets.append(
                get_Pset_of_Euqal_agg_end_point(tmp_iso, aggs, cost_param, mpc, next_eta_init, next_using_aggs))
            # print("next_Psets(1):", next_Psets)
            # print("first_after_skip_area", first_after_skip_area)
            if first_after_skip_area:
                next_Psets[-1] -= last_skip_Pset_range
                first_after_skip_area = False

            find_inflection_Pset_flag = False
            for i in range(len(next_Psets)):  # next_Pset must be greater than last_Pset
                if next_Psets[i] <= last_pset:
                    next_Psets[i] = 99999999
            # print(last_pset, "`s next_Psets", next_Psets)
            while not find_inflection_Pset_flag:
                next_Pset = min(next_Psets)
                # print("next_position:", next_Pset)

                # update eta_init
                tmp_iso.p_set = next_Pset
                _, next_eta_init = dispatch_with_inital_state(tmp_iso, aggs, cost_param, mpc, next_eta_init,
                                                              next_using_aggs)
                # print("next_eta_init:", next_eta_init)

                # update using_aggs
                next_V = []
                for i in range(len(aggs)):
                    if abs(1 - next_eta_init[i]) < 1e-6:
                        next_V.append(-1)
                    else:
                        next_V.append(value_estimate_function(tmp_iso, aggs, cost_param, next_eta_init, i))
                # print("next_V", next_V)

                max_tmp_V, next_using_aggs = max(next_V), []
                if max_tmp_V <= -1e-8:
                    empty_step_flag = True
                    next_using_aggs = []
                else:
                    for i in range(len(aggs)):
                        if abs(max_tmp_V - next_V[i]) < 1e-3:
                            next_using_aggs.append(i)
                # print("next_using_aggs", next_using_aggs)
                if next_using_aggs != current_using_aggs or last_pset != next_Pset:
                    find_inflection_Pset_flag = True
                else:
                    next_Psets.remove(next_Pset)

        if abs(next_Pset - positions[-1][0]) > 1e-4:
            positions.append([next_Pset, next_eta_init, next_using_aggs])

        # 所有agg容量已经都被完全调用的话，不会再改变了
        if abs(sum(next_eta_init) - len(aggs)) < 1e-6:
            break

    return positions


# ------ For get_next_inflection End------

def draw_two_stage_dynamic_dispatch_algorithm(iso, aggs, cost_param, mpc, sample_times=600, draw=False):
    draw_p_set, draw_etas, draw_costs, draw_costs_PAM = [], [], [], []
    for i in range(len(aggs)):
        draw_etas.append([])

    # 按顺序求所有分配规则变化的节点
    tmp_iso = copy.deepcopy(iso)
    positions = get_inflection_positions(iso, aggs, cost_param, mpc)
    positions.append([99999999999, [], []])
    print("Turning positions:")
    for i in range(len(positions)):
        print(positions[i])
    print("")
    # draw_critical_points(iso, aggs, positions)

    # 求解Pset_list中的每一个Pset下的最优功率分配
    Pset_list = [i / sample_times * iso.p_set for i in range(sample_times)]
    for i in range(len(Pset_list)):
        tmp_Pset, tmp_eta_init, tmp_using_aggs = Pset_list[i], [0 for k in range(len(aggs))], []
        stage_id = -1
        for j in range(1, len(positions)):
            if tmp_Pset <= positions[j][0]:
                stage_id = j - 1
                break
        # print(Pset_list[i], "MW stage_id", stage_id, tmp_Pset, positions[stage_id][0])
        cost, etas = 0, []
        if stage_id == len(positions) - 2:
            tmp_iso.p_set = tmp_Pset
            etas = positions[-2][1]
            cost = E_cost_function(aggs, tmp_iso, cost_param, etas)
        else:
            tmp_iso.p_set = tmp_Pset
            tmp_eta_init, tmp_using_aggs = positions[stage_id][1], positions[stage_id][2]
            cost, etas = dispatch_with_inital_state(tmp_iso, aggs, cost_param, mpc, tmp_eta_init, tmp_using_aggs)

        draw_p_set.append(Pset_list[i])
        for j in range(len(aggs)):
            draw_etas[j].append(etas[j])
        draw_costs.append(cost)
        # print(np.round(Pset_list[i], 2), "MW stage_id", stage_id, np.round(cost, 2), np.round(etas, 2))

    # 设置图形大小和线条颜色
    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot(111)
    lines = []

    # draw distribution result curves
    for i in range(len(aggs)):
        line, = ax1.plot(draw_p_set, draw_etas[i], label=f"$LA_{{{i + 1}}}$")
        lines.append(line)
    # ax1.set_title("Demand response assignment for " + str(len(aggs)) + " aggregator. 2-Stage")
    ax1.set_xlabel(r"$P_{AGC}$ (kW)")
    ax1.set_ylabel(r"Usage Proportions")
    ax1.set_xlim([0, iso.p_set + 1])
    xticks_pos1 = [i * 100 for i in range(int(iso.p_set / 100 + 1))]
    xticks_labels1 = [f"{int(tick)}" for tick in xticks_pos1]
    plt.xticks(xticks_pos1, xticks_labels1, rotation=0)

    ax1.set_ylim([0, 1])
    yticks_pos1 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    yticks_labels1 = [f"{int(tick * 100)}\%" for tick in yticks_pos1]
    ax1.set_yticks(yticks_pos1)
    ax1.set_yticklabels(yticks_labels1)

    # ax1.legend(handles=lines, loc='upper left', bbox_to_anchor=(1.1, 1), borderaxespad=0.)
    ax1.legend(handles=lines, loc='upper left')
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.subplots_adjust(right=0.8)

    # print("Costs:", costs)
    # print("Dispatching etas:", etas)
    name = []
    for i in range(len(aggs)):
        name.append("agg" + str(i + 1))
    test = pd.DataFrame(data=draw_etas).transpose()
    test.columns = name
    test.to_csv('./2stage_etas.csv', encoding='gbk')
    test = pd.DataFrame(data=draw_costs)
    test.columns = ["cost"]
    test.to_csv('./2stage_cost.csv', encoding='gbk')

    plt.savefig("./2stage-assignment.pdf", format="pdf")
    plt.show()


def two_stage_dynamic_dispatch_algorithm_time(iso, aggs, cost_param, mpc, sample_times=600):
    # 按顺序求所有分配规则变化的节点
    start_time = time.time()
    tmp_iso = copy.deepcopy(iso)
    positions = get_inflection_positions(iso, aggs, cost_param, mpc)
    positions.append([99999999999, [], []])
    end_time = time.time()
    stage1_time = end_time - start_time
    print("TEIRP-PDA (stage 1) cost", stage1_time, " s")

    # 求解Pset_list中的每一个Pset下的最优功率分配
    # Pset_list = [16, 410, 520, 560, 800, 900, 1200]
    Pset_list = [i / sample_times * iso.p_set for i in range(sample_times + 1)]
    start_time = time.time()
    costs = []
    for i in range(len(Pset_list)):
        tmp_Pset, tmp_eta_init, tmp_using_aggs = Pset_list[i], [0 for k in range(len(aggs))], []
        stage_id = -1
        for j in range(1, len(positions)):
            if tmp_Pset <= positions[j][0]:
                stage_id = j - 1
                break
        # print("stage_id", stage_id, tmp_Pset, positions[stage_id][0])
        cost, etas = 0, []
        if stage_id == len(positions) - 2:
            tmp_iso.p_set = tmp_Pset
            etas = positions[-2][1]
            cost = E_cost_function(aggs, tmp_iso, cost_param, etas)
        else:
            tmp_iso.p_set = tmp_Pset
            tmp_eta_init, tmp_using_aggs = positions[stage_id][1], positions[stage_id][2]
            cost, etas = dispatch_with_inital_state(tmp_iso, aggs, cost_param, mpc, tmp_eta_init, tmp_using_aggs)
        costs.append(float(cost))
    end_time = time.time()
    stage2_time = end_time - start_time
    print("TEIRP-PDA (stage 2) cost", stage2_time, " s")
    return stage1_time, stage2_time, stage1_time + stage2_time, costs


def test_gurobi_time(iso, aggs, cost_param, mpc, sample_times=600, method="gurobi"):
    Pset_list = [i / sample_times * iso.p_set for i in range(sample_times + 1)]
    tmp_iso = copy.deepcopy(iso)
    costs, total_t = [], 0
    for i in range(len(Pset_list)):
        tmp_Pset = Pset_list[i]
        tmp_iso.p_set = tmp_Pset
        cost, etas, t = assignment_gurobi(tmp_iso, aggs, cost_param, mpc, method)
        total_t = total_t + t
        costs.append(cost)
    return total_t, costs


def compare_time(iso, aggs, cost_param, mpc, my_sample_times, method, system_name="IEEE33"):
    # method参数可选值: "TEIRP-PDA" "gurobi" "cplex"
    print("Optimizer:", method, "System:", system_name)
    gurobi_times, gurobi_SOCP_times = [], []
    cplex_times, cplex_SOCP_times = [], []
    two_stage_times = []
    stage1_times, stage2_times = [], []

    if method == "TEIRP-PDA":
        stage1_time, stage2_time, two_stage_time, costs = two_stage_dynamic_dispatch_algorithm_time(iso, aggs,
                                                                                                    cost_param, mpc,
                                                                                                    sample_times=my_sample_times)
        two_stage_times.append(two_stage_time)
        stage1_times.append(stage1_time)
        stage2_times.append(stage2_time)
        print("TEIRP-PDA time: ", "Preparation stage: ", stage1_time, " Real-time stage: ", stage2_time,
              two_stage_time, "s")
        print("Costs", costs)
    elif method == "gurobi":
        gurobi_time, costs = test_gurobi_time(iso, aggs, cost_param, mpc, sample_times=my_sample_times,
                                              method=method)
        gurobi_times.append(gurobi_time)
        print("Gurobi time:", gurobi_time, "s")
        print("Costs", costs)
    elif method == "cplex":
        cplex_time, costs = test_gurobi_time(iso, aggs, cost_param, mpc, sample_times=my_sample_times,
                                             method=method)
        cplex_times.append(cplex_time)
        print("Cplex time:", cplex_time, "s")
        print("Costs", costs)


def two_aggregator():
    ##### ---  Set parameters ---
    cost_param = [100, 0, 0.1]
    mpc = case33()
    aggs, total_power = ([AGG(250, 19, 2, 0.8),
                          AGG(400, 16, 16, 0.25),
                          AGG(350, 4, 36, 0.6),
                          AGG(200, 4, 35, 0.2),
                          # AGG(200, 18, 2, 0.5)
                          ],
                         1200)
    iso = ISO(total_power, len(aggs))
    define_mpc(mpc, aggs, nid2aid={3: [0], 21: [1], 17: [2], 28: [3]})

    ##### --- Compare distribution result ---
    draw_assignment(aggs, total_power, iso, cost_param, mpc, sample_time=600, method="gurobi")
    draw_value_estimate_function(aggs, total_power, iso, cost_param, mpc, 600, "gurobi")

    # draw power assignment results calculated by TEIRP-PDA
    draw_two_stage_dynamic_dispatch_algorithm(iso, aggs, cost_param, mpc, sample_times=600)
    draw_value_estimate_function(aggs, total_power, iso, cost_param, mpc, 600, "TEIRP")

    ##### --- Compare calculation time ---
    my_sample_times = 150
    compare_time(iso, aggs, cost_param, mpc, my_sample_times, method="TEIRP-PDA", system_name="IEEE33")
    compare_time(iso, aggs, cost_param, mpc, my_sample_times, method="gurobi", system_name="IEEE33")
    compare_time(iso, aggs, cost_param, mpc, my_sample_times, method="cplex", system_name="IEEE33")


if __name__ == "__main__":
    # 全局默认设置
    os.environ['PATH'] += os.pathsep + '/path/to/latex'
    # np.random.seed(0)
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['font.size'] = 9  # 设置全局字体大小为9磅
    plt.rcParams['font.family'] = 'Times New Roman'  # 设置全局字体为新罗马
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams.update({'text.usetex': True})
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',  # 使用颜色编码定义颜色
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    # 常数参数设置
    over_response_ratio_low = -0.15
    over_response_ratio_up = 0.15
    V_root_sqr = 1.0  # root node的电压为固定的1.0pu,端点处电压>0.95pu
    V_max_sqr, V_min_sqr = 1.21, 0.81  # 各节点电压范围为1.05pu--0.95pu

    two_aggregator()
