% An example usage of MMP for prediction, general model
%==========================================================================
% Copyright (c) 2019 Hui Xiao
%==========================================================================
% 3-20-2019

% In this example, the MMP will predict future data (k steps)
addpath('functions')

%% define slow and fast sampling speed
clear
Tsf = 0.008;   % fast smapling period
L = 6;
Tss = Tsf*L;     % slow sampling period
Tsc = Tsf/10;  % sampling time for the continuous signal
n_amp = 0.001;   % noise amplitude
t_run = 20;    % simulation running time
%% define siganl frequency distribution
% constant model
A_const = [1 -1];
% linear model
A_linear = [1 -2 1];
% quadratic model
A_quad = [1 -3 3 -1];

Apara = A_quad;
%% Calculate prediction parameters
predic_step = 7;
m = length(Apara)-1;
p = m - 1; % give unique solution.
W = FIR_MMP(Apara, L, p, predic_step);

%% run simulation
sim('MMP_simulink_prediction_general_example');

%% plot results
figure,
plot(signal_c.time, signal_c.signals.values);
hold on,
stairs(measurement.time, measurement.signals.values);
stairs(FIR_output.time, FIR_output.signals.values(:));
legend('true signal','measurement','recovered signal');
xlabel('time (seconds)');
title('FIR-MMP');
