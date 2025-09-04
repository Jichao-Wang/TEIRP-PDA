function [frequency, time_out] = simulate_secondary_frequency(input_power_gap, time_delay1)
    %time = 0:0.1:10;
    %time = transpose(time);
    %sinewave=sin(time);
    %sig1=[time sinewave];
    
    %我自己的输入波形
    my_time = 1:1:size(input_power_gap,2) * 4;
    power_gap = 1:1:size(input_power_gap,2) * 4;
    time_interval = 0.01;
    not_adjust_power_rate = 0.3;
    sim_time = size(input_power_gap,2) * 6;
    %time_delay = 6;
    time_delay = time_delay1;

    %time_delay = 5;
    
    for i=0:1:(size(input_power_gap,2)-1)
        my_time(i*4+1)=i*6;
        power_gap(i*4+1)=input_power_gap(i+1);
        my_time(i*4+2)=i*6+time_delay;
        power_gap(i*4+2)=input_power_gap(i+1);
        
        my_time(i*4+3)=i*6+time_delay;
        power_gap(i*4+3)=not_adjust_power_rate*input_power_gap(i+1);
        my_time(i*4+4)=(i+1)*6;
        power_gap(i*4+4)=not_adjust_power_rate*input_power_gap(i+1);
    end

    my_time = transpose(my_time);
    power_gap= transpose(power_gap);
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
    print_input = out.sim_out;
    frequency = out.delta_f;

    %disp(out.sim_out)
    %disp(frequency)
    %disp("Finish!")
    pause(1);
end