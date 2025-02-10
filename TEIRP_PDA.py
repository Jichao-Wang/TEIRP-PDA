import numpy as np
from scipy.linalg import solve
from scipy.optimize import minimize
from scipy.special import gamma
from scipy.stats import beta

from gurobipy import *
import matlab  # matlab用于数据转换等工作
import matlab.engine  # engine启动调用的matlab函数

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.legend_handler import HandlerTuple
import matplotlib.patches as mpatches
import pandas as pd
import random
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


def assignment_gurobi(iso, aggs, cost_param):
    # 求ISO对每个AGG的最优功率分配策略

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

    # Add equality constraint
    m.addConstr(iso.p_set >= sum(etas[i] * aggs[i].capacity for i in range(len(aggs))), "c_e")

    # Add inequality constraints:
    for i in range(len(etas)):
        m.addConstr(etas[i] >= 0, "c_m" + str(i + 1))
        m.addConstr(etas[i] <= 1, "c_l" + str(i + 1))

    m.optimize()

    if m.status == 2:
        eta, cost = [], m.objVal
        for v in m.getVars():
            if abs(v.x) < 1e-6:
                eta.append(0)
            else:
                eta.append(v.x)
        # print("eta:", eta)
        # print('cost:', m.objVal)
        return cost, eta
    else:
        print("Not optimal solution from Gurobi.")
        return -1, []


def assignment_ito_gurobi(iso, aggs, cost_param):
    # 在方差随调用量线性变化的情况下，求ISO对每个AGG的最优功率分配策略.

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
        [aggs[i].capacity * aggs[i].alpha / (aggs[i].alpha + aggs[i].beta) * etas[i] for i in
         range(len(aggs))])) +
                   cost_param[2] * (pow(iso.p_set - sum(
        [aggs[i].capacity * aggs[i].alpha / (aggs[i].alpha + aggs[i].beta) * etas[i] for i in range(len(aggs))]),
                                        2) +
                                    sum(aggs[i].capacity * aggs[i].capacity * etas[i] * etas[i] * Var(aggs[i]) for i in
                                        range(len(aggs)))
                                    ) +

                   sum(aggs[i].price_param * aggs[i].capacity * etas[i] * aggs[i].alpha / (
                           aggs[i].alpha + aggs[i].beta) for i in range(len(aggs)))

                   , GRB.MINIMIZE)

    # Add equality constraint
    m.addConstr(iso.p_set >= sum(etas[i] * aggs[i].capacity for i in range(len(aggs))), "c_e")

    # Add inequality constraints:
    for i in range(len(etas)):
        m.addConstr(etas[i] >= 0, "c_m" + str(i + 1))
        m.addConstr(etas[i] <= 1, "c_l" + str(i + 1))

    m.optimize()

    if m.status == 2:
        eta, cost = [], m.objVal
        for v in m.getVars():
            if abs(v.x) < 1e-6:
                eta.append(0)
            else:
                eta.append(v.x)
        # print("eta:", eta)
        # print('cost:', m.objVal)
        return cost, eta
    else:
        print("Not optimal solution from Gurobi.")
        return -1, []


def draw_assignment(aggs, total_power, iso, cost_param, sample_time=200):
    p_set, etas, costs = [], [], []
    for i in range(len(aggs)):
        etas.append([])

    # dead_band = deadband(aggs, cost_param)
    # print("dead-band is at p_set =", dead_band, "MW")  # 死区：最低成本策略是不需求响应的区域

    for i in range(sample_time + 1):
        iso.set_total_power(i / sample_time * total_power)
        cost, eta = assignment_gurobi(iso, aggs, cost_param)
        p_set.append(i / sample_time * total_power)
        for j in range(len(aggs)):
            etas[j].append(eta[j])
        costs.append(cost)

    # 设置图形大小和线条颜色
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    lines = []
    for i in range(len(aggs)):
        line, = ax1.plot(p_set, etas[i], label="LA" + str(i + 1))
        lines.append(line)
    # ax1.set_title("Demand response assignment for " + str(len(aggs)) + " aggregator.")
    ax1.set_xlabel(r"$P_{AGC}$ (MW)")
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
    delta_y = 150
    ax2.set_ylim([0, 5 * delta_y])
    yticks_pos2 = [0, 1 * delta_y, 2 * delta_y, 3 * delta_y, 4 * delta_y, 5 * delta_y]
    yticks_labels2 = [f'{int(tick)}' for tick in yticks_pos2]
    ax2.set_yticks(yticks_pos2)
    ax2.set_yticklabels(yticks_labels2)

    ax1.legend(handles=lines, loc='upper left', bbox_to_anchor=(1.1, 1), borderaxespad=0.)
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.savefig("./assignment.svg", format="svg")
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


def simulate_cost_function(aggs, iso, cost_param, eta):
    # 采样法在某个随机变量取值下，计算实际ISO总成本
    # for Mentle Carlo 仿真
    delta_p, price = iso.p_set, 0
    for i in range(len(aggs)):
        samples = np.random.beta(aggs[i].alpha, aggs[i].beta, size=1)
        pi = aggs[i].capacity * eta[i] * (samples[0] * 0.2 + (1 - 0.10))
        delta_p -= pi
        price += pi * aggs[i].price_param

    cost = cost_param[0] + cost_param[1] * delta_p + cost_param[2] * delta_p * delta_p + price
    return cost


def E_cost_function(aggs, iso, cost_param, eta):
    # 已知所有参数和分配比例情况下，求E(Cost)
    delta_p = iso.p_set - sum(aggs[i].capacity * eta[i] * E(aggs[i]) for i in range(len(aggs)))
    delta_p2 = sum(aggs[i].capacity * aggs[i].capacity * eta[i] * eta[i] * Var(aggs[i]) for i in
                   range(len(aggs))) + delta_p * delta_p
    price = sum(aggs[i].price_param * aggs[i].capacity * eta[i] * E(aggs[i]) for i in range(len(aggs)))
    cost = cost_param[0] + cost_param[1] * delta_p + cost_param[2] * delta_p2 + price
    return cost


def draw_MonteCarlo_assignment(aggs, total_power, iso, cost_param):
    p_set, etas, costs = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000], [], []
    theory_costs = []
    for i in range(len(p_set)):
        costs.append([])

    for i in range(len(p_set)):
        iso.set_total_power(p_set[i])
        theory_cost, eta = assignment_gurobi(iso, aggs, cost_param)
        theory_costs.append(theory_cost)
        for j in range(10000):
            costs[i].append(simulate_cost_function(aggs, iso, cost_param, eta))
    mean_cost = [sum(costs[i]) / 10000 for i in range(len(p_set))]

    fig, ax = plt.subplots()

    h1 = ax.violinplot(costs, positions=p_set, showextrema=False, showmeans=True, widths=60)
    for pc in h1['bodies']:
        pc.set_facecolor('salmon')
        pc.set_edgecolor('black')
        pc.set_alpha(0.5)
    h1['cmeans'].set_color('tomato')
    h1_color = h1['bodies'][0].get_facecolor().flatten()
    h1_patch = mpatches.Patch(color=h1_color, label="Monte Carlo simulation")

    ax.scatter(p_set, mean_cost, color=h1_color, marker="|", zorder=2)
    h2 = ax.scatter(p_set, theory_costs, color="#17becf", marker="_",
                    label="Expectation of cost calculated by TEIRP-PDA",
                    zorder=3)
    # h2_patch = mpatches.Patch(color="r", label="Monte Carlo")
    h3 = ax.scatter(None, None, color="tomato", marker="_", label="Average cost of simulations")

    plt.ylim(0, 700)
    plt.xlabel(r"$P_{AGC}$ (MW)")
    plt.ylabel(r"Expectation of Cost (\$)")

    plt.xlim([-100, 1100])
    xticks_pos1 = [i * 100 for i in range(int(iso.p_set / 100 + 1))]
    xticks_labels1 = [f"{int(tick)}" for tick in xticks_pos1]
    plt.xticks(xticks_pos1, xticks_labels1, rotation=0)

    # ax.legend(handles=[h1_patch, h3, h2], loc="upper left")
    ax.legend(handles=[h3, h2], loc="upper left")
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.savefig("./Monte_Carlo_assignment.svg", format="svg")
    plt.show()
    return 0


def draw_assignment_pset(aggs, total_power, iso, cost_param, sample_time=200):
    ax_x, etas, costs = [], [], []
    for i in range(len(aggs)):
        etas.append([])

    origin_alpha0 = aggs[0].alpha
    for i in range(sample_time):
        aggs[0].set_beta(i * 0.03 * origin_alpha0)
        cost, eta = assignment_gurobi(iso, aggs, cost_param)
        ax_x.append(i * 0.03 * origin_alpha0)
        for j in range(len(aggs)):
            etas[j].append(eta[j])
        costs.append(cost)

    # 设置图形大小和线条颜色
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    lines = []
    for i in range(len(aggs)):
        line, = ax1.plot(ax_x, etas[i], label="agg" + str(i + 1))
        lines.append(line)
    ax1.set_title("Demand response assignment of total power as " + str(iso.p_set) + "MW with various $\\beta_0$ .")
    ax1.set_xlabel("Beta of the aggregator 1")
    ax1.set_ylabel("The proportion of each aggregator capacity called by ISO. (%)")
    ax1.set_xlim([0, int(ax_x[-1]) + 1])
    ax1.set_ylim([0, 1.1])

    ax2 = ax1.twinx()
    line, = ax2.plot(ax_x, costs, color='r', linestyle='--', label="cost")
    lines.append(line)
    ax2.set_ylabel(r"ISO cost expectation (\$)")
    ax2.set_ylim([0, 5500])

    ax1.legend(handles=lines, loc='lower right')
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.savefig("./assignment_pset" + str(iso.p_set) + ".svg", format="svg")
    plt.show()


def value_estimate_function(iso, aggs, cost_param, eta, i):
    # 边际价值
    # 衡量在不同总功率分配下的第i个聚合商的价值函数，即  E(Cost)/(eta[i]) 的偏导数
    tmp_b = aggs[i].capacity * eta[i] * Var(aggs[i]) - E(aggs[i]) * iso.p_set
    for j in range(len(aggs)):
        tmp_b += aggs[j].capacity * eta[j] * E(aggs[j]) * E(aggs[i])
    cost = E(aggs[i]) * (aggs[i].price_param - cost_param[1]) + cost_param[2] * 2 * tmp_b
    # cost *= aggs[i].capacity

    # another implementation for E(Cost)/(c[i]*eta[i]) 的偏导数
    # cost = 2 * cost_param[2] * (
    #         aggs[i].capacity * aggs[i].capacity * eta[i] * Var(aggs[i]) - aggs[i].capacity * E(aggs[i]) * (
    #         iso.p_set - sum(aggs[k].capacity * eta[k] * E(aggs[k]) for k in range(len(aggs))))) \
    #        + aggs[i].price_param * aggs[i].capacity * E(aggs[i]) - aggs[i].capacity * E(aggs[i]) * cost_param[1]

    # val为聚合商帮助ISO避免的损失，所以E(loss)函数对聚合商调用量的偏导是负数，意味着聚合商用得越多，ISO损失越小.聚合商的价值是正的。
    val = 0 - cost
    return val


def draw_value_estimate_function_local(aggs, total_power, iso, cost_param, sample_time=200):
    p_set, prices = [], []
    for i in range(len(aggs)):
        prices.append([])

    for i in range(sample_time + 1):
        iso.set_total_power(i / sample_time * total_power)
        # iso.set_total_power(313*2 + i / sample_time * total_power)
        assign_cost, eta = assignment_gurobi(iso, aggs, cost_param)
        p_set.append(i / sample_time * total_power)
        # p_set.append(313*2 + i / sample_time * total_power)
        # print(i, iso.p_set, eta, end="; ")
        for j in range(len(aggs)):
            prices[j].append(value_estimate_function(iso, aggs, cost_param, eta, j))
        #     print(value_estimate_function(iso, aggs, cost_param, eta, j), end=' ')
        # print("")
    # 设置图形大小和线条颜色
    fig, lines = plt.figure(), []
    ax1 = fig.add_subplot(111)
    for i in range(len(aggs)):
        line, = ax1.plot(p_set, prices[i], label="LA" + str(i + 1))
        lines.append(line)
    # ax1.set_title("Incremental rate of aggregators.")
    ax1.set_xlabel(r"$P_{AGC}$ (MW)")
    ax1.set_ylabel(r"Incremental rate")

    # y_bottom, y_top = -8, 3
    # yticks_pos1 = [i for i in range(y_bottom, y_top, 1)]
    # yticks_pos1.append(y_top)
    # yticks_labels1 = [f"{int(tick)}" for tick in yticks_pos1]
    # ax1.set_yticks(yticks_pos1)
    # ax1.set_yticklabels(yticks_labels1)
    # ax1.set_ylim([y_bottom, y_top])
    #
    # ax1.set_xlim([0, iso.p_set + 1])
    # xticks_pos1 = [i * 100 for i in range(int(iso.p_set / 100 + 1))]
    # xticks_labels1 = [f"{int(tick)}" for tick in xticks_pos1]
    # plt.xticks(xticks_pos1, xticks_labels1, rotation=0)

    ax1.legend(handles=lines, loc='lower right')
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.savefig("./prices.svg", format="svg")
    plt.show()

    # name = []
    # for i in range(len(aggs)):
    #     name.append("LA" + str(i + 1))
    # test = pd.DataFrame(data=prices).transpose()
    # test.columns = name
    # test.to_csv('./incremental_rate.csv', encoding='gbk')


def draw_value_estimate_function(aggs, total_power, iso, cost_param, sample_time=200):
    p_set, prices = [], []
    for i in range(len(aggs)):
        prices.append([])

    for i in range(sample_time + 1):
        iso.set_total_power(i / sample_time * total_power)
        assign_cost, eta = assignment_gurobi(iso, aggs, cost_param)
        p_set.append(i / sample_time * total_power)
        # print(i, iso.p_set, eta, end="; ")
        for j in range(len(aggs)):
            prices[j].append(-1.0 * value_estimate_function(iso, aggs, cost_param, eta, j))
        #     print(value_estimate_function(iso, aggs, cost_param, eta, j), end=' ')
        # print("")
    # 设置图形大小和线条颜色
    fig, lines = plt.figure(), []
    ax1 = fig.add_subplot(111)
    for i in range(len(aggs)):
        line, = ax1.plot(p_set, prices[i], label="LA" + str(i + 1))
        lines.append(line)
    # ax1.set_title("Incremental rate of aggregators.")
    ax1.set_xlabel(r"$P_{AGC}$ (MW)")
    ax1.set_ylabel(r"Incremental Rate")
    ax1.set_xlim([0, total_power])
    y_bottom, y_top = -40, 7
    yticks_pos1 = [np.round(i * 0.1, 2) for i in range(y_bottom, y_top, 5)]
    yticks_pos1.append(y_top * 0.1)
    yticks_labels1 = [f"{np.round(tick, 2)}" for tick in yticks_pos1]
    ax1.set_yticks(yticks_pos1)
    ax1.set_yticklabels(yticks_labels1)
    ax1.set_ylim([min(yticks_pos1), max(yticks_pos1)])

    ax1.set_xlim([0, iso.p_set + 1])
    xticks_pos1 = [i * 100 for i in range(int(iso.p_set / 100 + 1))]
    xticks_labels1 = [f"{int(tick)}" for tick in xticks_pos1]
    plt.xticks(xticks_pos1, xticks_labels1, rotation=0)

    ax1.legend(handles=lines, loc='upper right')
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.savefig("./prices.svg", format="svg")
    plt.show()

    name = []
    for i in range(len(aggs)):
        name.append("LA" + str(i + 1))
    test = pd.DataFrame(data=prices).transpose()
    test.columns = name
    test.to_csv('./incremental_rate.csv', encoding='gbk')


def draw_beta_distribution(aggs):
    ab_pairs = []
    for i in range(len(aggs)):
        ab_pairs.append((aggs[i].alpha, aggs[i].beta))

    x = np.linspace(0, 1, 1002)[1:-1]
    i = 0
    fig = plt.figure()
    for a, b in ab_pairs:
        dist = beta(a, b)
        y = dist.pdf(x)
        plt.plot(x, y, label=r'LA%d\ $\alpha=%.1f,\ \beta=%.1f$' % (i + 1, a, b))
        i += 1

    # 设置标题
    # plt.title(r'Probability density function of $\eta_{r}$.')
    # 设置 x,y 轴取值范围
    plt.xlim(0, 1)
    plt.ylim(0, 10)
    plt.ylabel("Probability density")
    plt.xlabel(r"$\eta_{r}$")
    plt.legend(loc="upper center")
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.savefig("./beta.svg", format="svg")
    plt.show()


def draw_uncertain_power(aggs):
    x_max = -1
    for i in range(len(aggs)):
        x_max = max(x_max, aggs[i].capacity)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    lines = []
    legends, legend_names = [], []

    for i in range(len(aggs)):
        tmp_agg = aggs[i]
        tmp_x, tmp_y, tmp_up, tmp_bottom = [], [], [], []
        for j in range(int(tmp_agg.capacity) + 1):
            interval = beta.interval(0.99, tmp_agg.alpha, tmp_agg.beta)
            tmp_x.append(j)  # j 为ISO设置的被调用容量
            tmp_y.append(j * (beta.mean(tmp_agg.alpha, tmp_agg.beta) * 0.2 + 0.9))
            tmp_bottom.append(j * (interval[0] * 0.2 + 0.9))
            tmp_up.append(j * (interval[1] * 0.2 + 0.9))

        line, = ax1.plot(tmp_x, tmp_y, linewidth=1, label="AGG" + str(i + 1))
        area = ax1.fill_between(tmp_x, tmp_bottom, tmp_up, alpha=0.3, label="AGG" + str(i + 1))
        legends.append((line, area))
        legend_names.append("LA" + str(i + 1))
        lines.append(line)

    # 设置 x,y 轴取值范围
    ax1.set_xlim(0, x_max)
    ax1.set_ylim(0, 350)
    ax1.set_xlabel("set response capacity")
    ax1.set_ylabel("real response capacity")
    ax1.legend(legends, legend_names, loc='upper left', handler_map={tuple: HandlerTuple(ndivide=None)})

    ax1.grid(color='gray', linestyle='--', linewidth=0.5)
    ax1.legend(legends, legend_names)

    plt.savefig("./called_power_uncertainty.svg", format="svg")
    plt.show()


def draw_beta_value_function_ISO(cost_param, p_sets, sample_time=200):
    # 每一个p_set绘制一张图,每个图展示在某一个p_set的最优分配下，每个聚合商的不确定性系数对E的影响
    # 即 在p_set固定的情况下，E/beta[i] 的偏导数

    for p in range(len(p_sets)):
        aggs, total_power = [AGG(200, 18, 2, 0.2), AGG(400, 36, 4, 0.3), AGG(400, 2, 2, 0.3)], p_sets[p]
        iso = ISO(total_power, len(aggs))
        eta, cost = assignment_gurobi(iso, aggs, cost_param)
        for i in range(1, 200):
            for k in range(len(aggs)):
                aggs[k].set_beta_keep_mu((p + 1) / len(p_sets))

    return 0


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
def get_equal_value_dispatch(iso, aggs, a, eta_init, using_aggs):
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

    # 根据\eta_{si}约束，使用线性规划求解r
    manual_LP, r = True, 0
    if manual_LP:
        # 手动求解LP
        tmp_v = 9999999999
        for i in range(len(formular_params_list)):
            tmp_up_bound = (1 - formular_params_list[i][1]) / formular_params_list[i][0]
            tmp_low_bound = (0 - formular_params_list[i][1]) / formular_params_list[i][0]
            if formular_params_list[i][0] < 0:
                tmp_exchange = tmp_up_bound
                tmp_up_bound = tmp_low_bound
                tmp_low_bound = tmp_exchange
            # print("tmp_r_range :", tmp_low_bound, "<= r <=", tmp_up_bound)
            if tmp_v >= tmp_up_bound:
                tmp_v = tmp_up_bound
        # print("tmp_r = ", tmp_v)
        r = min(tmp_v, 0)  # 此处的变量r为边际价值,即-r.限制边际价值(-r)的下限,进而避免\eta超过1.

    # 如果多聚合商等价的话，计算协同使用的比例
    Equivalent_distribution_ratios = []
    # print("r=", r)
    return r, Equivalent_distribution_ratios


def dispatch_with_inital_state(iso, aggs, cost_param, eta_init, using_aggs):
    '''
    思路：计算无约束情况下的最有极值，然后如果突破边界条件的话，考虑讲最靠谱的agg剔除，然后再优化
    '''
    N, etas = len(aggs), []

    # 挑选需要改变的几个agg求解
    # 求边际价格k
    pi, Equivalent_distribution_ratios = 0, []
    if len(using_aggs) >= 1:
        pi, Equivalent_distribution_ratios = get_equal_value_dispatch(iso, aggs, cost_param, eta_init, using_aggs)
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
def get_equal_value_inflection(aggs, a, eta_init, using_aggs, agg_j):
    # print("eta_init", eta_init)
    # print("using_aggs", using_aggs)
    # print("agg_j", agg_j)

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
    v1 = -tmp_formular_params_list_agg_j[1] / (tmp_formular_params_list_agg_j[0])  # right version
    # v1 = -tmp_formular_params_list_agg_j[1] / (tmp_formular_params_list_agg_j[0] + 1e-9)  # add 1e-9 for not 0, only for time
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

        ans_Pset = 9999999999
        if tmp_pset_low > tmp_pset_up or tmp_pset_up < 0:
            ans_Pset = 9999999999  # 如果无解，那么下一个拐点是using_aggs中最先用光容量的那个Pset
        else:
            # 手动找到能够让r相等的那个边界
            tmp_iso1 = ISO(tmp_pset_up, len(aggs))
            tmp_iso2 = ISO(tmp_low_bound, len(aggs))
            cost1, etas1 = dispatch_with_inital_state(tmp_iso1, aggs, a, eta_init, using_aggs)
            cost2, etas2 = dispatch_with_inital_state(tmp_iso2, aggs, a, eta_init, using_aggs)
            tmp_r1, tmp_r2 = [], []
            tmp_deltar1, tmp_deltar2 = [], []
            for i in range(len(aggs)):
                tmp_r1.append(value_estimate_function(tmp_iso1, aggs, a, etas1, i))
                tmp_r2.append(value_estimate_function(tmp_iso2, aggs, a, etas2, i))
            for i in range(len(aggs)):
                if i == agg_j:
                    continue
                else:
                    tmp_deltar1.append(abs(tmp_r1[i] - tmp_r1[agg_j]))
                    tmp_deltar2.append(abs(tmp_r2[i] - tmp_r2[agg_j]))
            if min(tmp_deltar1) > min(tmp_deltar2):
                ans_Pset = tmp_low_bound
            else:
                ans_Pset = tmp_pset_up

        # print("tmp_Pset =", tmp_pset_up)
        # \eta_{si} = a * r + b * Pset + c
        # max_r = (eta_init[agg_j] - tmp_formular_params_list_agg_j[1] * ans_Pset - tmp_formular_params_list_agg_j[2]) / \
        #         tmp_formular_params_list_agg_j[0]
        # print("max_r on ans_Pset", ans_Pset, " is", max_r)
    # print("\n")
    return ans_Pset, ans_eta_V_Pset_const


# (没有使用。目前可以解决r=0且\eta_{si}=1导致的turning point，还解决不了r<0的情况)
# 当前using_agg只有一个agg时，估计其被完全调用时(\eta_{s0}=1)所处的Pset
def get_one_agg_final_point(aggs, a, eta_init, using_aggs):
    tmp_i = using_aggs[0]
    eta_init[tmp_i] = 1.0

    # r = k1 * Pset + k2
    k1 = -2 * a[2] * E(aggs[tmp_i])
    k2 = (E(aggs[tmp_i]) * (aggs[tmp_i].price_param - a[1])
          + 2 * a[2] * aggs[tmp_i].capacity * 1 * Var(aggs[tmp_i]))
    for j in range(len(aggs)):
        k2 += 2 * a[2] * aggs[j].capacity * eta_init[j] * (E(aggs[j]) * E(aggs[tmp_i]))

    not_using_aggs = []
    for i in range(len(aggs)):
        if i != tmp_i:
            not_using_aggs.append(i)
    print("NLAS:", not_using_aggs)
    print("ALAS:", using_aggs)
    print("eta_init:", eta_init)

    # tmp_a * \eta_{si} = tmp_b * r + tmp_c * iso.p_set + tmp_const, 此处是指using_aggs中的第1个agg_id, r为price
    # 此处公式是按照Incremental rate推导的, (Incremental rate = -1.0 * \frac{\partial E(Cost)}{\partial c_i\eta_{si}})
    tmp_a = 2.0 * a[2] * aggs[tmp_i].capacity * Var(aggs[tmp_i]) \
            + 2.0 * a[2] * aggs[tmp_i].capacity * E(aggs[tmp_i]) * E(aggs[tmp_i])

    tmp_const = E(aggs[tmp_i]) * (a[1] - aggs[tmp_i].price_param) \
                - 2.0 * a[2] * E(aggs[tmp_i]) * \
                (sum(aggs[k].capacity * eta_init[k] * E(aggs[k]) for k in not_using_aggs))

    tmp_b = 1.0

    tmp_c = 2.0 * a[2] * E(aggs[tmp_i])

    # eta_si = t1 * Pset + t2
    t1 = (-tmp_b * k1 + tmp_c) / tmp_a
    t2 = (-tmp_b * k2 + tmp_const) / tmp_a

    # s.t. r<=0
    c1 = -k2 / k1
    print(c1, "<= Pset")

    # s.t. 0<=\eta_i<=1
    c2 = (1 - t2) / t1
    c3 = -t2 / t1
    print(c3, "<= Pset <=", c2)

    # s.t. sum(c_i\eta_{si}<=Pset)
    others = 0
    for i in range(len(aggs)):
        if i != tmp_i:
            others += aggs[i].capacity * eta_init[i]
    # (t1*Pset+t2)*aggs[tmp_i]<=Pset-others
    # (t1-1/aggs[tmp_i])Pset<=-1/aggs[tmp_i]*(others)-t2
    c4 = (-1 / aggs[tmp_i].capacity * (others) - t2) / (t1 - 1 / aggs[tmp_i].capacity)
    print("Pset <=", c4)

    # objective: max r   =>   min Pset
    next_Pset = min(c2, c4)  # 在允许的范围内最小的Pset就是LA_i容量用尽时，下一个拐点的位置

    # debug. verify increment at a certain turning point.
    print("\nCalculated next_Pset:", next_Pset)
    p_set = next_Pset
    print("r =", k1 * p_set + k2)
    eta_i = (tmp_b * (k1 * p_set + k2) + tmp_c * p_set + tmp_const) / tmp_a
    print("eta_i =", eta_i, "(should be 1.0)")

    p_set = 1200
    print("\nIf p_set=", p_set)
    r = k1 * p_set + k2
    print("r =", r)
    eta_i = (tmp_b * (k1 * p_set + k2) + tmp_c * p_set + tmp_const) / tmp_a
    # eta_i = (tmp_b * r + tmp_c * p_set + tmp_const) / tmp_a
    print("eta_i =", eta_i, "(should be 1.0)")

    return next_Pset


# 求解using_agg中有多个聚合商时，由于聚合商之间不能维持边际价值相等的状态导致的拐点
def get_equal_aggs_inflection_end(aggs, a, eta_init, using_aggs, agg_j):
    # 获取多个aggs等价的终点处的Pset
    # agg_j就是using_aggs中触碰上边界eta_{si}=1的那个聚合商,传入的using_aggs已经将其移除了

    # print("\nusing_aggs", using_aggs, agg_j)
    # print("eta_init", eta_init)

    not_using_aggs = []
    for i in range(len(eta_init)):
        if i not in using_aggs:
            not_using_aggs.append(i)
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
            # print("tmp_Pset_range :", tmp_pset_low, "<= r <=", tmp_pset_up, tmp_low_bound, tmp_up_bound)

        ans_Pset = 999999999
        if tmp_pset_low > tmp_pset_up or tmp_pset_up < 0:
            ans_Pset = 9999999999  # 如果无解，那么下一个拐点是using_aggs中最先用光容量的那个Pset
        else:
            # 这么做的原因：由于期望影响\eta_{si}超调欠调,
            # sum(c_i_eta_{si}, i\in LAS)的期望可能大于1也可能小于1，因此sum-Pset的符号不固定
            # 进而导致最有一个不等式的不等号的方向也不固定可能是up<Pset<low, 也可能是low<pset<up

            # 手动找到能够让r相等的那个边界
            tmp_iso1 = ISO(tmp_pset_up, len(aggs))
            tmp_iso2 = ISO(tmp_low_bound, len(aggs))
            cost1, etas1 = dispatch_with_inital_state(tmp_iso1, aggs, a, eta_init, using_aggs)
            cost2, etas2 = dispatch_with_inital_state(tmp_iso2, aggs, a, eta_init, using_aggs)
            tmp_r1, tmp_r2 = [], []
            tmp_deltar1, tmp_deltar2 = [], []
            for i in range(len(aggs)):
                tmp_r1.append(value_estimate_function(tmp_iso1, aggs, a, etas1, i))
                tmp_r2.append(value_estimate_function(tmp_iso2, aggs, a, etas2, i))
            for i in range(len(aggs)):
                if i == agg_j:
                    continue
                else:
                    tmp_deltar1.append(abs(tmp_r1[i] - tmp_r1[agg_j]))
                    tmp_deltar2.append(abs(tmp_r2[i] - tmp_r2[agg_j]))
            if min(tmp_deltar1) > min(tmp_deltar2):
                ans_Pset = tmp_low_bound
            else:
                ans_Pset = tmp_pset_up

        # print("tmp_Pset =", tmp_pset_up)
    # print("\n")
    return ans_Pset, ans_eta_V_Pset_const


# 求解有最优微增量的聚合商被完全调用不得不推出ALA时的Pset
def get_Pset_of_Euqal_agg_end_point(iso, aggs, cost_param, eta_init, using_aggs):
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
        # next_Pset = get_one_agg_final_point(aggs, cost_param, eta_init, using_aggs)
        # print("DEBUG next_Pset:", next_Pset)
        # print("----------------------")
        return next_Pset
    elif len(using_aggs) >= 1:
        # 获取多个等价agg达到某一个agg的容量上限处的Pset.(using_agg内部某个成员无法维持最高的边际价值导致拐点)
        # print("\n------finding Pset of multi-aggs end------")
        next_Psets = []
        for i in range(len(using_aggs)):
            tmp_eta_init = copy.deepcopy(eta_init)
            tmp_eta_init[using_aggs[i]] = 1.0
            tmp_using_aggs = [x for x in using_aggs if x != using_aggs[i]]
            next_Pset, b = get_equal_aggs_inflection_end(aggs, cost_param, tmp_eta_init, tmp_using_aggs,
                                                         using_aggs[i])
            tmp_iso1 = iso
            tmp_iso1.p_set = next_Pset
            next_Psets.append(next_Pset)
        # print("DEBUG next_Psets:", next_Psets, "len(using_aggs) >= 1")
        # print("------------\n")
        return min(next_Psets)
    else:
        return 999999999


# 求解当前ALA额外i空集时，NLA中某聚合商的微增量刚刚<0时的解（由于ALA为空，微增量<0就是最优微增量）
def get_next_value_eq_0(aggs, a, etas, i):
    # 当前eta_init下使得using_aggs为空，求下一次using_aggs不为空时，P_set的值
    a1 = -aggs[i].capacity * E(aggs[i]) * a[1] + 2 * a[2] * aggs[i].capacity * aggs[i].capacity * etas[i] * Var(
        aggs[i]) + 2 * a[2] * aggs[i].capacity * E(aggs[i]) * sum(
        aggs[k].capacity * etas[k] * E(aggs[k]) for k in range(len(aggs))) + aggs[i].price_param * aggs[i].capacity * E(
        aggs[i])
    b1 = 2 * a[2] * aggs[i].capacity * E(aggs[i])
    return a1 / b1


def get_inflection_positions(iso, aggs, cost_param):
    # Calculate Pset if existing aggs[j] could replace any one of using_aggs
    tmp_iso = copy.deepcopy(iso)
    positions = [[0, [0 for i in range(len(aggs))], []]]
    next_Pset, next_eta_init, next_using_aggs = 0, [], []

    # deadband round
    next_Pset, next_using_aggs = deadband(aggs, cost_param, return_flag=2)
    next_eta_init = [0 for i in range(len(aggs))]
    positions.append([next_Pset, next_eta_init, next_using_aggs])
    tmp_iso.p_set = next_Pset

    # --- debug one state(could skip deadbound or error rounds) ---
    # next_Pset, next_eta_init, next_using_aggs = 10, [0, 0.0, 0.0], [0]
    # next_Pset, next_eta_init, next_using_aggs = 100.68493150684928, [0.25057678442682035, 0.0, 0.0], [0]
    # next_Pset, next_eta_init, next_using_aggs = 120.68493150684928, [0.30, 0.0, 0.0], [0]
    # next_Pset, next_eta_init, next_using_aggs = 400, [1.0, 0, 0], [2]
    # next_Pset, next_eta_init, next_using_aggs = 540.1380670611471, [1.0, 0, 0.4671268902038265], [1, 2]
    # tmp_iso.p_set = next_Pset
    # ---------

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

    # next_Pset, eta_V_Pset_const = get_equal_value_inflection(aggs, cost_param, next_eta_init, next_using_aggs, 2)
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

            # --- for print ---
            # print("next_position", next_Pset)
            # print("next_eta_init", positions[-1][1])
            # next_V = []
            # for i in range(len(aggs)):
            #     if abs(1 - next_eta_init[i]) < 1e-6:
            #         next_V.append(-1)
            #     else:
            #         next_V.append(value_estimate_function(tmp_iso, aggs, cost_param, next_eta_init, i))
            # print("next_V", next_V)
            # print("next_using_aggs", next_using_aggs, "\n")
            # ------
        else:
            # 当前using_aggs不为空
            # print("\n------current_Pset", next_Pset, wjc_it, "--------")
            # set current state
            not_using_aggs, full_used_aggs = [], []
            for j in range(len(next_eta_init)):
                if j not in next_using_aggs:
                    not_using_aggs.append(j)
                if abs(next_eta_init[j] - 1) < 1e-5:
                    full_used_aggs.append(j)
            current_using_aggs = copy.deepcopy(next_using_aggs)
            # print("not_using_aggs:", not_using_aggs)
            # update Pset
            next_Psets = []
            for i in range(len(aggs)):
                if i in not_using_aggs:
                    # print(next_eta_init, next_using_aggs, i)
                    next_Pset, _ = get_equal_value_inflection(aggs, cost_param, next_eta_init, next_using_aggs, i)
                    # print("tmp_next_Pset1", next_Pset)
                    # verify Value of aggs on the end position to avoid numbercial computation accuracy
                    if next_Pset < 999999:
                        tmp_iso1 = copy.deepcopy(iso)
                        tmp_iso1.p_set = next_Pset
                        _, tmp_next_eta_init = dispatch_with_inital_state(tmp_iso1, aggs, cost_param, next_eta_init,
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
                        if abs(tmp_next_V[next_using_aggs[0]] - tmp_next_V[i]) > 1e-3:
                            next_Psets.append(9999999999)
                        else:
                            next_Psets.append(next_Pset)
                    else:
                        next_Psets.append(8888888888)
                else:
                    next_Psets.append(7777777777)
            # print("next_Psets(0):", next_Psets)
            next_Psets.append(
                get_Pset_of_Euqal_agg_end_point(tmp_iso, aggs, cost_param, next_eta_init, next_using_aggs))
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
                _, next_eta_init = dispatch_with_inital_state(tmp_iso, aggs, cost_param, next_eta_init, next_using_aggs)
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

def draw_two_stage_dynamic_dispatch_algorithm(iso, aggs, cost_param, sample_times=600, draw=False):
    draw_p_set, draw_etas, draw_costs = [], [], []
    for i in range(len(aggs)):
        draw_etas.append([])

    # 按顺序求所有分配规则变化的节点
    tmp_iso = copy.deepcopy(iso)
    positions = get_inflection_positions(iso, aggs, cost_param)
    # positions = [[0, [0, 0, 0], []],
    #              [10.0, [0, 0, 0], [0]],
    #              [400.0, [1, 0, 0], [2]],
    #              [540, [1, 0, 0.35], [1, 2]],
    #              [720, [1, 1, 0.3], [2]]]
    positions.append([99999999999, [], []])
    print("Critical positions:")
    for i in range(len(positions)):
        print(positions[i])
    print("")
    # draw_critical_points(iso, aggs, positions)

    # # example:手动求解优化问题实现一次功率分配
    # print("------Test one-time dispatch------")
    # iso.p_set, eta_init, using_aggs = 544, [1.0, 0, 0.4671268902038232], [1, 2]
    # cost, etas = dispatch_with_inital_state(iso, aggs, cost_param, eta_init, using_aggs)
    # print(iso.p_set, cost, etas)
    # print("------------")

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
            cost, etas = dispatch_with_inital_state(tmp_iso, aggs, cost_param, tmp_eta_init, tmp_using_aggs)
            # if abs(tmp_Pset - 700) < 1:
            #     print("debug", tmp_iso.p_set, tmp_eta_init, tmp_using_aggs)
        # print(etas)
        # print("MyCode", i, cost, etas)
        # cost, etas = assignment_gurobi(tmp_iso, aggs, cost_param)
        # print("Gurobi:", cost, etas)

        draw_p_set.append(Pset_list[i])
        for j in range(len(aggs)):
            draw_etas[j].append(etas[j])
        draw_costs.append(cost)
        # print(np.round(Pset_list[i], 2), "MW stage_id", stage_id, np.round(cost, 2), np.round(etas, 2))

    # 设置图形大小和线条颜色
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    lines = []
    for i in range(len(aggs)):
        line, = ax1.plot(draw_p_set, draw_etas[i], label="LA" + str(i + 1))
        lines.append(line)
    # ax1.set_title("Demand response assignment for " + str(len(aggs)) + " aggregator. 2-Stage")
    ax1.set_xlabel(r"$P_{AGC}$ (MW)")
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

    ax2 = ax1.twinx()
    line, = ax2.plot(draw_p_set, draw_costs, color='#17becf', linestyle='--', label="cost")
    lines.append(line)
    ax2.set_ylabel(r"Expectation of Cost (\$)   ")
    delta_y = 150
    ax2.set_ylim([0, 5 * delta_y])
    yticks_pos2 = [0, 1 * delta_y, 2 * delta_y, 3 * delta_y, 4 * delta_y, 5 * delta_y]
    yticks_labels2 = [f'{int(tick)}' for tick in yticks_pos2]
    ax2.set_yticks(yticks_pos2)
    ax2.set_yticklabels(yticks_labels2)

    # ax1.legend(handles=lines, loc='upper left', bbox_to_anchor=(1.1, 1), borderaxespad=0.)
    ax1.legend(handles=lines, loc='upper left')
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.subplots_adjust(right=0.8)

    plt.savefig("./2stage-assignment.svg", format="svg")
    plt.show()

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


def two_stage_dynamic_dispatch_algorithm_time(iso, aggs, cost_param, sample_times=600):
    # 按顺序求所有分配规则变化的节点
    start_time = time.time()
    tmp_iso = copy.deepcopy(iso)
    positions = get_inflection_positions(iso, aggs, cost_param)
    positions.append([99999999999, [], []])
    end_time = time.time()
    stage1_time = end_time - start_time
    print("two_stage_dynamic_dispatch_algorithm get_inflection_positions cost", end_time - start_time, " s")
    # print("positions:")
    # for i in range(len(positions)):
    #     print(positions[i])

    # 手动求解优化问题实现一次功率分配
    # example:
    # iso.p_set, eta_init, using_aggs = 304 * 2, [1.0, 0, 0.4666565162244413], [1, 2]
    # cost, etas = dispatch_with_inital_state(iso, aggs, cost_param, eta_init, using_aggs)

    # 求解Pset_list中的每一个Pset下的最优功率分配
    # Pset_list = [16, 410, 520, 560, 800, 900, 1200]
    Pset_list = [i / sample_times * iso.p_set for i in range(sample_times)]
    start_time = time.time()
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
            cost, etas = dispatch_with_inital_state(tmp_iso, aggs, cost_param, tmp_eta_init, tmp_using_aggs)

        # print("MyCode", cost, etas)
        # cost, etas = assignment_gurobi(tmp_iso, aggs, cost_param)
        # print("Gurobi:", cost, E_cost_function(aggs, tmp_iso, cost_param, etas), etas)
    end_time = time.time()
    stage2_time = end_time - start_time
    print("two_stage_dynamic_dispatch_algorithm stage 2 cost", end_time - start_time, " s")
    return stage1_time, stage2_time


def two_stage_dynamic_dispatch_algorithm_one_time_test(iso, aggs, cost_param):
    # 按顺序求所有分配规则变化的节点
    tmp_iso = copy.deepcopy(iso)
    positions = get_inflection_positions(iso, aggs, cost_param)
    positions.append([99999999999, [], []])

    # print("positions:")
    # for i in range(len(positions)):
    #     print(positions[i])

    # 手动求解优化问题实现一次功率分配
    tmp_Pset, tmp_eta_init, tmp_using_aggs = iso.p_set, [0 for k in range(len(aggs))], []
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
        cost, etas = dispatch_with_inital_state(tmp_iso, aggs, cost_param, tmp_eta_init, tmp_using_aggs)
    return cost, etas


def test_gurobi_time(iso, aggs, cost_param, sample_times=600):
    Pset_list = [i / sample_times * iso.p_set for i in range(0, sample_times, 2)]
    tmp_iso = copy.deepcopy(iso)
    for i in range(len(Pset_list)):
        tmp_Pset = Pset_list[i]
        tmp_iso.p_set = tmp_Pset
        cost, etas = assignment_gurobi(tmp_iso, aggs, cost_param)
    return 0


def generate_mass_aggs(N):
    aggs = []
    prices = np.random.uniform(0.1, 10.0, N)
    k = 0
    for i in range(10):
        for j in range(int(N / 10)):
            aggs.append(AGG(10, i + 1, j + 1, prices[k]))
            k += 1
    return aggs


def draw_time(cost_param, my_sample_times):
    gurobi_times = []
    two_stage_times = []
    stage1_times, stage2_times = [], []

    for i in range(1, 9):
        print("Number of AGGs:", int(i * 10))
        aggs = generate_mass_aggs(int(i * 10))

        total_power = i * 10 * 10
        iso = ISO(total_power, len(aggs))
        start_time = time.time()
        stage1_time, stage2_time = two_stage_dynamic_dispatch_algorithm_time(iso, aggs, cost_param,
                                                                             sample_times=my_sample_times)
        end_time = time.time()
        two_stage_time = end_time - start_time
        start_time = time.time()
        test_gurobi_time(iso, aggs, cost_param, sample_times=my_sample_times)
        end_time = time.time()
        gurobi_time = end_time - start_time

        gurobi_times.append(gurobi_time)
        two_stage_times.append(two_stage_time)
        stage1_times.append(stage1_time)
        stage2_times.append(stage2_time)
        print(i, stage1_time, stage2_time, gurobi_time, "s")
    print(stage1_times, stage2_times, gurobi_times)

    a, b, c = gurobi_times, stage1_times, stage2_times
    # 设置 x 轴标签
    labels = [f'{(i + 1) * 10}' for i in range(len(a))]

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    plt.figure()

    # 绘制 a 的柱子
    plt.bar(x - width / 2, a, label='Interior point method', width=0.3)

    # 绘制 b 和 c 组合的柱子
    bottom_c = [0] * len(c)
    for i in range(len(c)):
        if i == 0:
            plt.bar(x[i] + width / 2, c[i], bottom=bottom_c[i], label='TEIRP-PDA', color='red',
                    width=0.3)
            bottom_c[i] += c[i]
        else:
            plt.bar(x[i] + width / 2, c[i], bottom=bottom_c[i], color='red', width=0.3)
            bottom_c[i] += c[i]
    # for i in range(len(b)):
    #     if i == 0:
    #         plt.bar(x[i] + width / 2, b[i], bottom=bottom_c[i], label='Stage 1', color='green', width=0.3)
    #     else:
    #         plt.bar(x[i] + width / 2, b[i], bottom=bottom_c[i], color='green', width=0.3)

    # 设置 x 轴标签和标题
    plt.xticks(range(len(a)), labels)
    plt.xlabel('Number of Load Aggregators')
    plt.ylabel('Execution Time of Determining Real-time Power Distribution (s)')

    # 添加图例
    plt.legend()

    # 显示图形
    plt.savefig("./time.svg", format="svg")
    plt.show()
    return 0


def continue_use_id_cost(p_set, aggs, cost_param):
    iso = ISO(p_set, len(aggs))
    cost_init, eta_init = two_stage_dynamic_dispatch_algorithm_one_time_test(iso, aggs, cost_param)
    # print("pset_init", iso.p_set, "cost_init", cost_init, "eta_init", eta_init)

    # 计算强制调用某一agg的话，Cost的变化
    cost1, cost3, cost_optimal = [], [], []
    for i in range(11):
        p = i * 0.1
        iso.p_set = p_set + p
        tmp_cost, _ = two_stage_dynamic_dispatch_algorithm_one_time_test(iso, aggs, cost_param)
        cost_optimal.append(tmp_cost)
    eta = copy.deepcopy(eta_init)
    for i in range(11):
        p = i * 0.1
        # eta[1] = eta_init[1] + p / aggs[1].capacity
        iso.p_set = p_set + p
        # cost1.append(E_cost_function(aggs, iso, cost_param, eta) - cost_optimal[i])
        cost, eta = dispatch_with_inital_state(iso, aggs, cost_param, eta, [1])
        cost1.append(cost - cost_optimal[i])
    eta = copy.deepcopy(eta_init)
    for i in range(11):
        p = i * 0.1
        iso.p_set = p_set + p
        cost, eta = dispatch_with_inital_state(iso, aggs, cost_param, eta, [3])
        cost3.append(cost - cost_optimal[i])

    # print(np.round(cost1, 8))
    # print(np.round(cost3, 8))
    # print(np.round(cost_optimal, 5))

    cost_optimal = [0 for i in range(11)]

    # 设置图形大小和线条颜色
    fig, lines = plt.figure(), []
    ax2 = fig.add_subplot(111)
    p_sets = np.round([p_set + i * 0.1 for i in range(11)], decimals=2)

    line, = ax2.plot(p_sets, cost1, label="Still use LA" + str(2), marker='d', color="#ff7f0e")
    lines.append(line)
    ax2.set_xlabel(r"$P_{AGC}$ (MW)")
    ax2.set_ylabel(r"Expectation of cost (\$)   ")
    ax2.set_xlim([p_set, p_set + 0.1])
    delta_y = 1e-7
    ax2.set_ylim([0, 14 * delta_y])
    yticks_pos2 = [0, 2 * delta_y, 4 * delta_y, 6 * delta_y, 8 * delta_y, 10 * delta_y, 12 * delta_y, 14 * delta_y]
    yticks_labels2 = [f'{int(tick)}' for tick in yticks_pos2]
    ax2.set_yticks(yticks_pos2)
    ax2.set_yticklabels(yticks_labels2)

    ax1 = ax2.twinx()
    # line, = ax1.plot(p_sets, cost1, label="Still use LA" + str(2), marker='d')
    # lines.append(line)
    line, = ax1.plot(p_sets, cost3, label="Still use LA" + str(4), marker='s', color="#d62728")
    lines.append(line)

    ax1.set_ylabel(r"Expectation of Cost1 (\$)")
    y_bottom, y_top = 0, 7
    delta_y = 1e-4
    yticks_pos1 = [0, 1 * delta_y, 2 * delta_y, 3 * delta_y, 4 * delta_y, 5 * delta_y, 6 * delta_y, 7 * delta_y]
    yticks_labels1 = [f"{np.round(tick, 6)}" for tick in yticks_pos1]
    ax1.set_yticks(yticks_pos1)
    ax1.set_yticklabels(yticks_labels1)
    ax1.set_ylim([0, 7 * delta_y])

    xticks_pos1 = [p_set + i * 0.1 for i in range(11)]
    xticks_labels1 = [f"{np.round(tick, 9)}" for tick in xticks_pos1]
    plt.xticks(xticks_pos1, xticks_labels1, rotation=0)

    line, = ax2.plot(p_sets, cost_optimal, label="optimal solution", marker='o', color="#17becf")
    lines.append(line)

    # ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax2.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax2.legend(handles=lines, loc='upper left')
    ax2.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.savefig("./local_id_cost.svg", format="svg")
    plt.show()

    return 0


def continue_use_id_marginal_value(p_set, aggs, cost_param):
    iso = ISO(p_set, len(aggs))
    p_sets = np.round([p_set + i * 0.1 for i in range(11)], decimals=2)
    cost_init, eta_init = two_stage_dynamic_dispatch_algorithm_one_time_test(iso, aggs, cost_param)
    # print("pset_init", iso.p_set, "cost_init", cost_init, "eta_init", eta_init)

    # 计算强制调用某一agg的话，所有agg的Incremental rate的变化
    prices = [[] for i in range(len(aggs))]
    for i in range(11):
        iso.p_set = p_set + i * 0.1
        _, eta = two_stage_dynamic_dispatch_algorithm_one_time_test(iso, aggs, cost_param)
        for j in range(len(aggs)):
            prices[j].append(value_estimate_function(iso, aggs, cost_param, eta, j))

    # 设置图形大小和线条颜色
    fig, lines = plt.figure(), []
    ax1 = fig.add_subplot(111)

    # line, = ax1.plot(p_sets, prices[0], label="LA" + str(1), marker='d')
    # lines.append(line)
    line, = ax1.plot(p_sets, prices[1], label="LA" + str(2), marker='d', color="#ff7f0e")
    lines.append(line)
    line, = ax1.plot(p_sets, prices[2], label="LA" + str(3), marker='s', color="#2ca02c")
    lines.append(line)
    # line, = ax1.plot(p_sets, prices[3], label="LA" + str(4), marker='d')
    # lines.append(line)
    # line, = ax1.plot(p_sets, prices[4], label="LA" + str(5), marker='d')
    # lines.append(line)

    # ax1.set_title("Incremental rate of aggregators.")
    ax1.set_xlabel(r"$P_{AGC}$ (MW)")
    ax1.set_ylabel(r"Incremental rate (\$)")

    # y_bottom, y_top = -1, 9
    # yticks_pos1 = [i * 0.001 for i in range(y_bottom, y_top, 1)]
    # yticks_pos1.append(y_top)
    # yticks_labels1 = [f"{np.round(tick, 3)}" for tick in yticks_pos1]
    # ax1.set_yticks(yticks_pos1)
    # ax1.set_yticklabels(yticks_labels1)
    ax1.set_ylim([0.326, 0.33195])

    ax1.set_xlim([p_set, p_set + 0.1])
    xticks_pos1 = [p_set + i * 0.1 for i in range(11)]
    xticks_labels1 = [f"{np.round(tick, 2)}" for tick in xticks_pos1]
    plt.xticks(xticks_pos1, xticks_labels1, rotation=0)

    ax1.legend(handles=lines, loc='lower left')
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.savefig("./local_marginal_value_optimal.svg", format="svg")
    plt.show()
    iso.p_set = p_set

    # 计算强制调用某一agg的话，所有agg的Incremental rate的变化
    prices1 = [[] for i in range(len(aggs))]
    eta = copy.deepcopy(eta_init)
    for i in range(11):
        eta[1] = eta_init[1] + (i * 0.1) / aggs[1].capacity
        for j in range(len(aggs)):
            prices1[j].append(value_estimate_function(iso, aggs, cost_param, eta, j))
    # 设置图形大小和线条颜色
    fig, lines = plt.figure(), []
    ax1 = fig.add_subplot(111)

    # line, = ax1.plot(p_sets, prices[0], label="LA" + str(1), marker='d')
    # lines.append(line)
    line, = ax1.plot(p_sets, prices1[1], label="LA" + str(2), marker='d', color="#ff7f0e")
    lines.append(line)
    line, = ax1.plot(p_sets, prices1[2], label="LA" + str(3), marker='s', color="#2ca02c")
    lines.append(line)
    # line, = ax1.plot(p_sets, prices[3], label="LA" + str(4), marker='d')
    # lines.append(line)
    # line, = ax1.plot(p_sets, prices[4], label="LA" + str(5), marker='d')
    # lines.append(line)

    # ax1.set_title("Incremental rate of aggregators.")
    ax1.set_xlabel(r"$P_{AGC}$ (MW)")
    ax1.set_ylabel(r"Incremental rate (\$)")

    # y_bottom, y_top = -1, 9
    # yticks_pos1 = [i * 0.001 for i in range(y_bottom, y_top, 1)]
    # yticks_pos1.append(y_top)
    # yticks_labels1 = [f"{np.round(tick, 3)}" for tick in yticks_pos1]
    # ax1.set_yticks(yticks_pos1)
    # ax1.set_yticklabels(yticks_labels1)
    ax1.set_ylim([0.326, 0.33195])

    ax1.set_xlim([p_set, p_set + 0.1])
    xticks_pos1 = [p_set + i * 0.1 for i in range(11)]
    xticks_labels1 = [f"{np.round(tick, 2)}" for tick in xticks_pos1]
    plt.xticks(xticks_pos1, xticks_labels1, rotation=0)

    ax1.legend(handles=lines, loc='lower left')
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.savefig("./local_marginal_value_still_use_LA2.svg", format="svg")
    plt.show()

    # 计算强制调用某一agg的话，所有agg的Incremental rate的变化
    prices2 = [[] for i in range(len(aggs))]
    eta = copy.deepcopy(eta_init)
    for i in range(11):
        eta[2] = eta_init[2] + (i * 0.1) / aggs[2].capacity
        for j in range(len(aggs)):
            prices2[j].append(value_estimate_function(iso, aggs, cost_param, eta, j))
    # 设置图形大小和线条颜色
    fig, lines = plt.figure(), []
    ax1 = fig.add_subplot(111)

    # line, = ax1.plot(p_sets, prices[0], label="LA" + str(1), marker='d')
    # lines.append(line)
    line, = ax1.plot(p_sets, prices2[1], label="LA" + str(2), marker='d', color="#ff7f0e")
    lines.append(line)
    line, = ax1.plot(p_sets, prices2[2], label="LA" + str(3), marker='s', color="#2ca02c")
    lines.append(line)
    # line, = ax1.plot(p_sets, prices[3], label="LA" + str(4), marker='d')
    # lines.append(line)
    # line, = ax1.plot(p_sets, prices[4], label="LA" + str(5), marker='d')
    # lines.append(line)

    # ax1.set_title("Incremental rate of aggregators.")
    ax1.set_xlabel(r"$P_{AGC}$ (MW)")
    ax1.set_ylabel(r"Incremental rate (\$)")

    # y_bottom, y_top = -1, 9
    # yticks_pos1 = [i * 0.001 for i in range(y_bottom, y_top, 1)]
    # yticks_pos1.append(y_top)
    # yticks_labels1 = [f"{np.round(tick, 3)}" for tick in yticks_pos1]
    # ax1.set_yticks(yticks_pos1)
    # ax1.set_yticklabels(yticks_labels1)
    ax1.set_ylim([0.326, 0.33195])

    ax1.set_xlim([p_set, p_set + 0.1])
    xticks_pos1 = [p_set + i * 0.1 for i in range(11)]
    xticks_labels1 = [f"{np.round(tick, 2)}" for tick in xticks_pos1]
    plt.xticks(xticks_pos1, xticks_labels1, rotation=0)

    ax1.legend(handles=lines, loc='lower left')
    ax1.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.savefig("./local_marginal_value_still_use_LA3.svg", format="svg")
    plt.show()


# 绘制critical points
def draw_critical_points(iso, aggs, positions):
    global colors
    positions[-1][0] = iso.p_set
    for i in range(len(aggs)):
        label_str = "LA" + str(int(i + 1))
        plt.plot([0], [0], color=colors[i], label=label_str)
        for j in range(len(positions) - 1):
            tmp_x = np.arange(positions[j][0], positions[j + 1][0])
            if i in positions[j][2]:
                tmp_y = np.full_like(tmp_x, i + 1)
            else:
                tmp_y = np.zeros_like(tmp_x)
            plt.plot(tmp_x, tmp_y, color=colors[i], linewidth=3)

    plt.xlabel(r"$P_{ACE}$")
    plt.ylabel("Elements in ALA")
    plt.xlim(0, positions[-1][0])
    plt.ylim(0, len(aggs) + 0.5)

    y_bottom, y_top = 0, int(len(aggs))
    yticks_pos = [int(i + 1) for i in range(0, y_top)]
    yticks_labels = [f"LA{int(tick)}" for tick in yticks_pos]
    yticks_pos.append(len(aggs) + 0.5)
    yticks_labels.append(" ")
    plt.yticks(yticks_pos, yticks_labels)

    grid_x_values = [positions[i][0] for i in range(1, len(positions) - 1)]
    # print("grid_x_values", grid_x_values)
    xticks_pos = grid_x_values
    xticks_labels = [f"{np.round(tick, 0)}" for tick in xticks_pos]
    plt.xticks(xticks_pos, xticks_labels)

    plt.legend()
    for grid_x_value in grid_x_values:
        plt.axvline(x=grid_x_value, color='gray', linestyle='--', linewidth=1)
    plt.grid(axis="y", color='gray', linestyle='--', linewidth=1)
    plt.savefig("./critical_points.svg", format="svg")
    plt.show()

    positions[-1][0] = 99999999999


def frequency_regulation_simulator1():
    eng = matlab.engine.start_matlab()  # 启动
    eng.addpath(r'E:\wjc_code\matlab\aggregator', nargout=0)  # 确保 MATLAB 能找到函数文件所在路径

    # 仿真参数
    time_delay_gurobi = 1.0
    time_delay_TEIRP_PDA = 0.05
    total_FR_steps = 150  # 负荷波动次数

    # 生成仿真输入序列
    np.random.seed(0)
    power_gap = matlab.double(np.random.uniform(-300, 300, total_FR_steps).tolist())
    # print(power_gap)

    frequency_gurobi, time_out_gurobi = eng.simulate_secondary_frequency(power_gap, time_delay_gurobi,
                                                                         nargout=2)  # 调用 MATLAB 函数
    frequency_TEIRP_PDA, time_out_TEIRP_PDA = eng.simulate_secondary_frequency(power_gap, time_delay_TEIRP_PDA, nargout=2)  # 调用 MATLAB 函数

    frequency_gurobi = np.array(frequency_gurobi).transpose().tolist()[0]
    frequency_TEIRP_PDA = np.array(frequency_TEIRP_PDA).transpose().tolist()[0]
    # print(frequency_gurobi)
    print(len(time_out_TEIRP_PDA))
    fd_gurobi, fd_TEIRP_PDA = 0, 0
    for i in range(len(frequency_gurobi)):
        fd_gurobi += abs(frequency_gurobi[i])
        fd_TEIRP_PDA += abs(frequency_TEIRP_PDA[i])
    fd_gurobi /= len(frequency_gurobi)
    fd_TEIRP_PDA /= len(frequency_TEIRP_PDA)
    print("Frequency divation (gurobi, TEIRP_PDA, difference):", fd_gurobi, fd_TEIRP_PDA, fd_gurobi - fd_TEIRP_PDA)
    print("gurobi", max(frequency_gurobi), min(frequency_gurobi))
    print("TEIRP_PDA", max(frequency_TEIRP_PDA), min(frequency_TEIRP_PDA))

    for i in range(len(frequency_gurobi)):
        frequency_gurobi[i] += 50
        frequency_TEIRP_PDA[i] += 50
    plt.plot(frequency_gurobi, label="Interior point method")
    plt.plot(frequency_TEIRP_PDA, label="TEIRP-PDA")
    plt.xlim(0, 90001)
    plt.ylim(49.90, 50.10)
    # 自定义横坐标刻度和标签
    ticks = [0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000]  # 横坐标刻度
    labels = ['0', '100', '200', '300', '400', '500', '600', '700', '800', '900']  # 横坐标标签
    plt.xticks(ticks, labels)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.legend()
    # plt.legend(handles=lines, loc='upper left', bbox_to_anchor=(1.1, 1), borderaxespad=0.)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.savefig("./primary_frequency.svg", format="svg")
    plt.show()

    eng.quit()


# 验证目标函数的Hessian矩阵是否为正定矩阵(对\eta_{si}、\eta_{sj}求偏导)
def verify_zhengding_matrix1(aggs):
    matrix = np.zeros([len(aggs), len(aggs)])
    for i in range(len(aggs)):
        for j in range(len(aggs)):
            if i != j:
                matrix[i][j] = aggs[i].capacity * aggs[j].capacity * (E(aggs[i]) * E(aggs[j]))
            else:
                matrix[i][j] = aggs[i].capacity * aggs[j].capacity * (E(aggs[i]) * E(aggs[j]) + Var(aggs[i]))
    eigenvalues = np.linalg.eigvals(matrix)
    print("eigenvalues", eigenvalues)


# 验证目标函数的Hessian矩阵是否为正定矩阵(对c_i\eta_{si}、c_j\eta_{sj}求偏导)
def verify_zhengding_matrix2(aggs):
    matrix = np.zeros([len(aggs), len(aggs)])
    for i in range(len(aggs)):
        for j in range(len(aggs)):
            if i != j:
                matrix[i][j] = E(aggs[i]) * E(aggs[j])
            else:
                matrix[i][j] = E(aggs[i]) * E(aggs[j]) + Var(aggs[i])
    eigenvalues = np.linalg.eigvals(matrix)
    print("eigenvalues", eigenvalues)


def two_aggregator():
    # Set parameters
    cost_param = [100, 0, 0.1]
    aggs, total_power = ([AGG(250, 19, 2, 0.8),
                          AGG(400, 16, 16, 0.25),
                          AGG(350, 4, 36, 0.6),
                          AGG(200, 4, 35, 0.2),
                          # AGG(200, 18, 2, 0.5)
                          ],
                         1200)
    iso = ISO(total_power, len(aggs))

    # draw power assignment results calculated by Gurobi
    draw_assignment(aggs, total_power, iso, cost_param, sample_time=600)

    # draw power assignment results calculated by TEIRP-PDA
    draw_two_stage_dynamic_dispatch_algorithm(iso, aggs, cost_param, sample_times=600)
    # draw turning points of stage 1 in TEIRP-PDA
    positions = get_inflection_positions(iso, aggs, cost_param)
    draw_critical_points(iso, aggs, positions)

    # Compare the calculation speed for coordinate mass aggregators by Gurobi or TEIRP-PDA
    my_sample_times = 150
    draw_time(cost_param, my_sample_times)

    # Simulate the frequency
    frequency_regulation_simulator1()


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

    over_response_ratio_low = -0.15
    over_response_ratio_up = 0.15
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',  # 使用颜色编码定义颜色
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    two_aggregator()
