function [frequency, time_out] = s3(input_power_gap, time_delay1)
    %time = 0:0.1:10;
    %time = transpose(time);
    %sinewave=sin(time);
    %sig1=[time sinewave];
    
    %我自己的输入波形
    my_time = 1:1:size(input_power_gap,2) * 1;
    power_gap = 1:1:size(input_power_gap,2) * 1;
    time_interval = 0.001;
    time_delay = 6;
    sim_time = size(input_power_gap,2) * time_delay;
    
    
    for i=0:1:(size(input_power_gap,2)-1)
        my_time(i+1)=i*time_delay;
        power_gap(i+1)=input_power_gap(i+1);
    end

    my_time = transpose(my_time);
    power_gap= transpose(power_gap);
    %disp("my_time");
    %disp(my_time);
    %disp("power_gap");
    %disp(power_gap);
    sig1=[my_time, power_gap];
    assignin('base', 'sig1', sig1);
    %sig1=timeseries(my_time, power_gap);
    
    %设置模型
    open_system('untitled.slx');
    set_param('untitled','StopTime', num2str(sim_time));
    set_param('untitled', 'FixedStep', num2str(time_interval)); % 设置仿真步长为 0.01 秒
    out = sim('untitled', 'ExternalInput', 'sig1'); % 运行仿真并传入 timeseries 信号
    assignin('base', 'out', out);

    time_out = out.tout;
    loads = out.sim_out;
    frequency = out.delta_f;

    frequency_AGG1000_TEIRP = out.delta_f1;
    frequency_AGG2000_TEIRP = out.delta_f4;
    frequency_AGG3000_TEIRP = out.delta_f7;

    frequency_AGG1000_Gurobi = out.delta_f2;
    frequency_AGG2000_Gurobi = out.delta_f5;
    frequency_AGG3000_Gurobi = out.delta_f8;
    
    frequency_AGG1000_Cplex = out.delta_f3;
    frequency_AGG2000_Cplex = out.delta_f6;
    frequency_AGG3000_Cplex = out.delta_f9;

    csvwrite('load_curve.csv',loads);
    csvwrite('frequency_AGG1000_TEIRP.csv',frequency_AGG1000_TEIRP);
    csvwrite('frequency_AGG2000_TEIRP.csv',frequency_AGG2000_TEIRP);
    csvwrite('frequency_AGG3000_TEIRP.csv',frequency_AGG3000_TEIRP);
    csvwrite('frequency_AGG1000_Gurobi.csv',frequency_AGG1000_Gurobi);
    csvwrite('frequency_AGG2000_Gurobi.csv',frequency_AGG2000_Gurobi);
    csvwrite('frequency_AGG3000_Gurobi.csv',frequency_AGG3000_Gurobi);
    csvwrite('frequency_AGG1000_Cplex.csv',frequency_AGG1000_Cplex);
    csvwrite('frequency_AGG2000_Cplex.csv',frequency_AGG2000_Cplex);
    csvwrite('frequency_AGG3000_Cplex.csv',frequency_AGG3000_Cplex);


    %disp(out.sim_out)
    %disp(frequency)
    %disp("Finish!")
    pause(1);
end