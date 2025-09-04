import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.gridspec import GridSpec
import os


def original_frequency_results():
    data = [9.900780059, 11.22178551, 10.44253303, 10.30340564, 11.27018496, 10.45581931,
            10.5663736, 11.28376409, 9.974117153, 10.71183697, 11.00465759, 10.13680767, 10.68862851, 9.979983408,
            10.441461, 9.567194156, 9.265638742, 10.5401265, 10.21264374, 9.493338916, 9.395470633, 10.1462746,
            10.25484254,
            9.488836133, 9.349020353, 10.41916051, 10.28242659, 9.560491744, 9.646716922, 10.25370528, 10.27088464,
            9.618937228, 9.357610337, 10.38650168, 10.26964124, 10.62086659, 10.623358, 10.30020829, 11.87248361,
            10.23849399, 10.52462743, 10.56951479, 12.4520958, 10.68253358, 10.12266039, 11.01760066, 12.23343309,
            10.21911111, 11.85296062, 10.65185912, 10.99791329, 10.18083782, 9.267705906, 9.459927931, 8.885869445,
            8.730776355, 8.377185396, 8.37459824, 7.053527488, 6.88907463, 6.800127851, 6.788604762, 6.780311834,
            6.773226984, 6.804317317, 6.759192567, 6.701215826, 6.677180942, 6.648686264, 6.594050581, 6.609235122,
            6.591235822, 6.544716436, 6.634445039, 6.733773983, 6.932774024, 7.103967695, 7.228406273, 8.123058055,
            8.333789708, 11.12013179, 8.839586454, 9.897866859, 8.497900117, 9.375699501, 9.029283696,
            10.05620768, 9.374457982, 9.358119952, 11.14429691, 10.00724546, 11.01394128, 9.74399278, 10.13607058,
            11.46583064, 9.984530494]

    for i in range(len(data)):
        if data[i] >= 8:
            data[i] = 200 + data[i]
        else:
            data[i] = data[i] / 8 * 200
    print(data)
    plt.figure(figsize=(12, 6))  # 设置图形大小（宽12英寸，高6英寸）

    # 绘制折线图
    # range(len(data)) 生成X轴（数据点的索引），data 是Y轴的值
    plt.plot(range(len(data)), data,
             marker='o',  # 数据点标记为圆圈
             linestyle='-',  # 线型为实线
             linewidth=1,  # 线宽
             markersize=3,  # 标记点大小
             color='blue',  # 线条颜色
             label='Load (MW)')  # 图例标签

    # 添加标题和标签
    plt.title('Secondary Frequency Control - Load Profile', fontsize=14)
    plt.xlabel('Time Index', fontsize=12)
    plt.ylabel('Load (MW)', fontsize=12)

    # 添加网格线，使读数更方便
    plt.grid(True, linestyle='--', alpha=0.7)

    # 添加图例
    plt.legend()

    # 自动调整布局以避免标签被截断
    plt.tight_layout()

    # 显示图形
    plt.show()


def read_csv_to_list_pandas(filename, a=50):
    # 读取CSV文件，header=None表示没有列标题
    df = pd.read_csv(filename, header=None)
    if a == 50:
        df += a
    else:
        df /= 10
    # 将第一列转换为列表
    return df[0].tolist()


def revised1_frequency_results():
    load = read_csv_to_list_pandas("frequency_results/load_curve.csv", a=0)
    f_TEIRP_1 = read_csv_to_list_pandas("frequency_results/frequency_AGG1000_TEIRP.csv")
    f_TEIRP_2 = read_csv_to_list_pandas("frequency_results/frequency_AGG2000_TEIRP.csv")
    f_TEIRP_3 = read_csv_to_list_pandas("frequency_results/frequency_AGG3000_TEIRP.csv")

    f_Gurobi_1 = read_csv_to_list_pandas("frequency_results/frequency_AGG1000_Gurobi.csv")
    f_Gurobi_2 = read_csv_to_list_pandas("frequency_results/frequency_AGG2000_Gurobi.csv")
    f_Gurobi_3 = read_csv_to_list_pandas("frequency_results/frequency_AGG3000_Gurobi.csv")

    f_Cplex_1 = read_csv_to_list_pandas("frequency_results/frequency_AGG1000_Cplex.csv")
    f_Cplex_2 = read_csv_to_list_pandas("frequency_results/frequency_AGG2000_Cplex.csv")
    f_Cplex_3 = read_csv_to_list_pandas("frequency_results/frequency_AGG3000_Cplex.csv")

    print("Max Frequency Deviation(Hz), Min Frequency Deviation(Hz), Average Aboslute Frequency Deviation(Hz)")
    print("1000 AGGs")
    print(max(f_TEIRP_1), min(f_TEIRP_1), sum(abs(y) for y in f_TEIRP_1) / len(f_TEIRP_1))
    print(max(f_Gurobi_1), min(f_Gurobi_1), sum(abs(y) for y in f_Gurobi_1) / len(f_Gurobi_1))
    print(max(f_Cplex_1), min(f_Cplex_1), sum(abs(y) for y in f_Cplex_1) / len(f_Cplex_1))
    print("2000 AGGs")
    print(max(f_TEIRP_2), min(f_TEIRP_2), sum(abs(y) for y in f_TEIRP_2) / len(f_TEIRP_2))
    print(max(f_Gurobi_2), min(f_Gurobi_2), sum(abs(y) for y in f_Gurobi_2) / len(f_Gurobi_2))
    print(max(f_Cplex_2), min(f_Cplex_2), sum(abs(y) for y in f_Cplex_2) / len(f_Cplex_2))
    print("3000 AGGs")
    print(max(f_TEIRP_3), min(f_TEIRP_3), sum(abs(y) for y in f_TEIRP_3) / len(f_TEIRP_3))
    print(max(f_Gurobi_3), min(f_Gurobi_3), sum(abs(y) for y in f_Gurobi_3) / len(f_Gurobi_3))
    print(max(f_Cplex_3), min(f_Cplex_3), sum(abs(y) for y in f_Cplex_3) / len(f_Cplex_3))

    x = list(range(1, len(f_Cplex_3) + 1))
    for i in range(len(x)):
        x[i] = x[i] / len(x) * 900

    # 使用GridSpec创建3行1列布局
    fig = plt.figure(figsize=(8, 12))
    gs = GridSpec(4, 1, height_ratios=[1, 1, 1, 1])  # 等高

    axes = [fig.add_subplot(gs[i, 0]) for i in range(4)]

    # 设置每个子图的宽高比为6:1
    for ax in axes:
        ax.set_box_aspect(1 / 4)  # 6:1比例

    # 绘制第一张图
    axes[0].plot(x, load, color='blue', linewidth=1)
    # axes[0].set_title('Net Load Curve')
    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('Power (MW)')
    axes[0].set_xlim(0, 900)
    axes[0].set_ylim(0, 25)
    axes[0].grid(color='gray', linestyle='--', linewidth=0.5, dashes=(5, 5))

    axes[1].plot(x, f_TEIRP_1, label='TEIRP-PDA', color='blue', linewidth=1)
    axes[1].plot(x, f_Gurobi_1, label='IPM(Gurobi)', color='red', linewidth=1)
    axes[1].plot(x, f_Cplex_1, label='IPM(Cplex)', color='green', linewidth=1)
    # axes[1].set_title('Coordinating 1000 Aggregators')
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylabel('Frequency (Hz)')
    axes[1].set_xlim(0, 900)
    axes[1].set_ylim(49.90, 50.025)
    axes[1].legend()
    axes[1].grid(color='gray', linestyle='--', linewidth=0.5, dashes=(5, 5))

    axes[2].plot(x, f_TEIRP_2, label='TEIRP-PDA', color='blue', linewidth=1)
    axes[2].plot(x, f_Gurobi_2, label='IPM(Gurobi)', color='red', linewidth=1)
    axes[2].plot(x, f_Cplex_2, label='IPM(Cplex)', color='green', linewidth=1)
    # axes[2].set_title('Coordinating 2000 Aggregators')
    axes[2].set_xlabel('Time (s)')
    axes[2].set_ylabel('Frequency (Hz)')
    axes[2].set_xlim(0, 900)
    axes[2].set_ylim(49.90, 50.025)
    axes[2].legend()
    axes[2].grid(color='gray', linestyle='--', linewidth=0.5, dashes=(5, 5))

    axes[3].plot(x, f_TEIRP_3, label='TEIRP-PDA', color='blue', linewidth=1)
    axes[3].plot(x, f_Gurobi_3, label='IPM(Gurobi)', color='red', linewidth=1)
    axes[3].plot(x, f_Cplex_3, label='IPM(Cplex)', color='green', linewidth=1)
    # axes[3].set_title('Coordinating 3000 Aggregators')
    axes[3].set_xlabel('Time (s)')
    axes[3].set_ylabel('Frequency (Hz)')
    axes[3].set_xlim(0, 900)
    axes[3].set_ylim(49.90, 50.025)
    axes[3].legend()
    axes[3].grid(color='gray', linestyle='--', linewidth=0.5, dashes=(5, 5))

    plt.savefig("./frequency.png", format="png", dpi=300)

    # 分别保存每个子图为PDF
    titles = [
        'Net_Load_Curve',
        'Coordinating_1000_Aggregators',
        'Coordinating_2000_Aggregators',
        'Coordinating_3000_Aggregators'
    ]
    for i, ax in enumerate(axes):
        # 创建新的图形
        fig_single = plt.figure(figsize=(8, 3))
        ax_single = fig_single.add_subplot(111)

        # 复制原图内容
        lines = ax.get_lines()
        for line in lines:
            x_data = line.get_xdata()
            y_data = line.get_ydata()
            label = line.get_label()
            color = line.get_color()
            linewidth = line.get_linewidth()
            ax_single.plot(x_data, y_data, label=label, color=color, linewidth=linewidth)

        # 复制标题和标签
        # ax_single.set_title(ax.get_title())
        ax_single.set_xlabel(ax.get_xlabel())
        ax_single.set_ylabel(ax.get_ylabel())
        ax_single.set_xlim(ax.get_xlim())
        ax_single.set_ylim(ax.get_ylim())

        # 复制图例和网格
        if ax.get_legend() is not None:
            ax_single.legend()
        ax_single.grid(color='gray', linestyle='--', linewidth=0.5, dashes=(5, 5))

        # 设置相同的宽高比
        ax_single.set_box_aspect(1 / 4)

        # 调整布局并保存
        plt.tight_layout()
        plt.savefig(f"./{titles[i]}.pdf", format="pdf", dpi=300, bbox_inches='tight')
        plt.close(fig_single)

    plt.show()


os.environ['PATH'] += os.pathsep + '/path/to/latex'
plt.rcParams['axes.unicode_minus'] = False
# plt.rcParams['font.size'] = 9  # 设置全局字体大小为9磅
plt.rcParams['font.family'] = 'Times New Roman'  # 设置全局字体为新罗马
# Case 13
revised1_frequency_results()
