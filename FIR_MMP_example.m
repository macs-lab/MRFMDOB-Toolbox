% This is an example usage of multi-rate model-based predictor. 
clear
addpath(genpath('functions'));
% define parameters
L = 4;
% generate time series
Tsf = 0.1;  % fast sampling period.
Tss = Tsf*L; % slow sampling period.
steps = 100;
end_time = steps*Tsf*L;
d_time = 0:Tsf:end_time;
dL_time = 0:Tss:end_time;
dc_time = 0:Tsf/10:end_time;
%% generate signal
f = [1.2 1.67 2.4]; %Hz
A = [1 0.8 0.3];
phi = [0.6 0.3 1.4];
b = [0 0 0];
noise_amp = 0.05;
[dc,d,dL] = SensorHarmonic(f,A,phi,b,dc_time,d_time,dL_time,noise_amp);

figure,
plot(dc_time, dc);
hold on
stairs(d_time, d);
stairs(dL_time, dL);
legend('continus signal','fast sampled true value','slow sampled signal');

%% solve predicting parameters
Apara = Apara_prd(f,Tsf);
m = length(Apara)-1;
p = m -1; % give unique solution.
W = FIR_MMP(Apara, L, p);
%% recover the signal.
d_est = FIR_MMP_est(dL,W);
%% plot results        
figure,
plot(dc_time,dc,'-.');
hold on
plot(dL_time,dL,'ro','MarkerSize',10,'LineWidth',2);
stairs(d_time(1:steps*L),d_est(1:steps*L),'LineWidth',2,'Color',[0.49 0.18 0.56])
legend('d_c(t)','dL[n]','d[n]')
xlabel('time')

figure,
subplot(2,1,1);
plot(dc_time,dc,'b-.');
hold on
stairs(dL_time,dL,'r','LineWidth',1.5);
legend('d_c(t)','dL[n]')

subplot(2,1,2);
p21 = plot(dc_time,dc,'b-.');
hold on
p22 = stairs(d_time(1:steps*L),d_est(1:steps*L),'Color',[0.49 0.18 0.56],'LineWidth',2);
legend('d_c(t)','d[n]')
xlabel('time')