clear
% define parameters
L = 4;
% generate time series
T = 0.1;  % fast sampling period.
Tu = T*L; % slow sampling period.
steps = 100;
end_time = steps*T*L;
d_time = 0:T:end_time;
dL_time = 0:Tu:end_time;
dc_time = 0:T/10:end_time;
%% generate signal
f = [1.2 1.67 2.4];
A = [1 0.8 0.3];
phi = [0.6 0.3 1.4];
b = [0 0 0];
noise_amp = 0;
[dc,d,dL] = SensorHarmonic(f,A,phi,b,dc_time,d_time,dL_time,noise_amp);

figure,
plot(dc_time, dc);
hold on
plot(d_time, d);
plot(dL_time, dL);

%% solve predicting parameters
Apara = Apara_prd(f,T);
m = length(Apara)-1;
p = m -1; % give unique solution.
W = MMP(Apara, L, p);
%% recover the signal.
d_est = MMP_est(dL,W);
%% plot results        
figure,
plot(dc_time,dc,'-.');
hold on
plot(dL_time,dL,'ro','MarkerSize',10,'LineWidth',2);
stairs(d_time,d_est,'LineWidth',2,'Color',[0.49 0.18 0.56])
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
p22 = stairs(d_time,d_est,'Color',[0.49 0.18 0.56],'LineWidth',2);
legend('d_c(t)','d[n]')
xlabel('time')