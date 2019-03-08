% An example usage of MMP simulink modual.
%==========================================================================
% Copyright (c) 2019 Hui Xiao
%==========================================================================
% 3-7-2019

%% define slow and fast sampling speed
clear
Tsf = 0.008;   % fast smapling period
L = 6;
Tss = Tsf*L;     % slow sampling period
Tsc = Tsf/10;  % sampling time for the continuous signal
n_amp = 0.005;   % noise amplitude
t_run = 20;    % simulation running time
%% define siganl frequency distribution
f = [1.2 4];
Amp = [1 0.8];
Phi = [0.6 0];

%% Calculate prediction parameters
[B,a] = IIR_MMP(f,L,Tss,0.95);  % IIR parameters

Apara = Apara_prd(f,Tsf); % calculate the signal model
m = length(Apara)-1;
p = m - 1; % give unique solution.
W = FIR_MMP(Apara, L, p, 3);

%% run simulation
sim('MMP_simulink_prediction_example');

%% plot results
figure,
plot(signal_c.time, signal_c.signals.values);
hold on,
stairs(measurement.time, measurement.signals.values);
stairs(FIR_output.time, FIR_output.signals.values);
legend('true signal','measurement','recovered signal');
xlabel('time (seconds)');
title('FIR-MMP');

figure,
plot(signal_c.time, signal_c.signals.values);
hold on,
stairs(measurement.time, measurement.signals.values);
stairs(IIR_output.time, IIR_output.signals.values);
legend('true signal','measurement','recovered signal');
xlabel('time (seconds)');
title('IIR-MMP');