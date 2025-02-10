import matlab  # matlab用于数据转换等工作
import matlab.engine  # engine启动调用的matlab函数
import numpy as np
import matplotlib.pyplot as plt


def matlab_v(K=300, T=10.0, current_delta_f=0, load=[]):
    # K = 300 # 输入功率
    # T = 6.0 # 时间常数
    # current_delta_f # 当前频率偏差=0Hz

    eng = matlab.engine.start_matlab()  # 启动

    # 定义惯性环节
    num = matlab.double([0, 0.01 * K])  # 惯性环节分子
    den = matlab.double([T, 40])  # 惯性环节分母(Ts+D)
    t_step = 0.001
    t = matlab.double(np.arange(0, 2.0 + t_step, t_step))

    tmp = eng.step(num, den, t, nargout=3)  # c, x, t = eng.step(num, den, t, nargout=3)
    c = np.array(tmp[0]).flatten().tolist()
    t = np.array(tmp[2]).flatten().tolist()
    # 调用 step 函数
    print("c", len(c), c)
    print("t", len(t), t)
    # print(c[int(T / t_step)])  # 获取时间常数T处的响应，此处应为63.21%
    # print(c[int((T + 2.389) / t_step)])  # 获取时间常数T+ 2.389处的响应，此处应为
    # print(c[int((T + 0.648) / t_step)])  # 获取时间常数T+ 0.648处的响应，此处应为

    # 绘制惯性环节阶跃响应结果曲线
    fig, ax = plt.subplots()
    ax.plot(t, c, label='T=' + str(T))
    # ax.plot(t, c2, label='T=10')
    # ax.plot(t, x, label='Input')
    ax.set_title('First Order System')
    ax.legend(loc=0)
    plt.show()

    print(c[int((0.2389) / t_step)])  # 获取时间常数T- 2.389处的响应，此处应为
    print(c[int((0.04) / t_step)])  # 获取时间常数T- 0.648处的响应，此处应为

    eng.quit()  # 停止 或eng.exit()
    return 0


def calculate_delta_f(T=10.0, current_delta_f=0, load=[], time=[]):
    # T = 6.0 # 时间常数
    # current_delta_f # 当前频率偏差=0Hz

    eng = matlab.engine.start_matlab()  # 启动

    matlab_variable_name = 'input_load'
    eng.putvar('base', matlab_variable_name, load)

    # 假设时间向量为t，和上面的python_array长度相同
    eng.eval("input_timeseries = timeseries(input_load, time);", nargout=0)

    model_name = 'your_model_name'
    sim_results = eng.sim(model_name)
    print("sim_results", sim_results)

    # 定义惯性环节
    # num = matlab.double([0, 0.01])  # 惯性环节分子
    # den = matlab.double([T, 40])  # 惯性环节分母(Ts+D)
    # t_step = 0.001
    # t = matlab.double(np.arange(0, 2.0 + t_step, t_step))

    # tmp = eng.step(num, den, t, nargout=3)  # c, x, t = eng.step(num, den, t, nargout=3)
    # c = np.array(tmp[0]).flatten().tolist()
    # t = np.array(tmp[2]).flatten().tolist()

    # 绘制惯性环节阶跃响应结果曲线
    # fig, ax = plt.subplots()
    # ax.plot(t, c, label='T=' + str(T))
    # # ax.plot(t, c2, label='T=10')
    # # ax.plot(t, x, label='Input')
    # ax.set_title('First Order System')
    # ax.legend(loc=0)
    # plt.show()

    eng.quit()  # 停止 或eng.exit()
    return 0


time_step_num = int(6 / 0.0001)
delta_loads, time = [], []
delta_load = 300  # 初始负荷偏差
delta_f = 0  # 初始频率偏差

for i in range(time_step_num):
    time.append(i * 0.0001)
    if i < int(0.012 / 0.0001):
        delta_loads.append(delta_load)
    else:
        delta_loads.append(1 / 6 * delta_load)
calculate_delta_f(T=10.0, current_delta_f=delta_f, load=delta_loads, time=time)
