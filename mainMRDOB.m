%Given Ts, L and disturbance frequencys, implement MR-FMDOB techniques to reject
%disturbance.
%Created by Hui Xiao, 10-13-2015

clear
%load plant model
load('Tz.mat')
load('Tss.mat')
%% input arguments
Tu = input('Tu = ? seconds, default is 2.5e-5: ');
if isempty(Tu)
    Tu = 2.5e-5;
end
PdL=d2d(Tz,Tu);
L = input('L = ?, default is 2: ');
if isempty(L)
    L = 2;
end
Ts = Tu*L;
Nyquist = 1/Ts/2;
distN = input('how many freqency components does the disturbance has?(default 1)');
if isempty(distN)
    distN=1;
end
freq = zeros(1,distN);
for i=1:distN
    try
        freq(i) = input(['Please enter the ',num2str(i),'th frequency value(Hz), or created a random one: ']);
    catch
        if freq(i) == 0
            freq(i)= rand*2*Nyquist;
        end
    end
end
disp(['the frequeny has ',num2str(distN),' different components:'])
for i=1:distN
    disp([num2str(i),': ',num2str(freq(i)),'Hz'])
end
% give random amplitude and phase value
dist.freq  = zeros(1,distN);
dist.phase = zeros(1,distN);
dist.amp   = zeros(1,distN);
for i=1:distN
    dist.freq(i)  = freq(i);
    dist.phase(i) = rand*pi;
    dist.amp(i)   = 0.5+0.5*rand;
end
%% Forward model Youla Q design >> direct approach
[stdBP,~] = lattice_prod(freq,2,Tu);
phaF = phaseCompFilter_prod(...
    freqresp(Tz,freq*2*pi),...
    freq,Tu);
QFM = stdBP*phaF;
if 1
    xbodeplot(QFM)
end

%% Calculate the pridictor parameters
PRpara = PRpara_prd(freq,Tu,L);

%%
sim('MRFMDOB.slx')
%{
---> Sometimes when choosing a large L, there would be the following errer
messages:
Error using mainMRDOB (line 62)
Invalid setting in 'MRFMDOB/Constant' for parameter 'Value'.
Caused by:
    Error using mainMRDOB (line 62)
    Error evaluating parameter 'Value' in 'MRFMDOB/Constant'        Error
    using mainMRDOB (line 62)
        Undefined function or variable 'PRpara'.
---> try to manually use clear command before running.
%}
%%
%{
figure,plot(output_c.time,output_c.signals.values,'y--');
hold on
stairs(output_fastsampled.time,output_fastsampled.signals.values,'k');
stairs(output_slowsampled.time,output_slowsampled.signals.values,'b');
legend('continuous output','fastsampled output','slowsampled output');
title('System Output');

figure, stairs(measured.time,measured.signals.values);
hold on
plot(distur.time,distur.signals.values); 
stairs(recovered.time,recovered.signals.values,'k--'); 
legend('measured disturbance','real disturbance','reconstruted disturbance');
title('Disturbance Recovery');

figure, stairs(plant_output.time,plant_output.signals.values,'--');
hold on
stairs(control.time,control.signals.values/10);
plot(distur.time,distur.signals.values); 
plot(output_c.time,output_c.signals.values);
legend('plant output','control singals/10','disturbance','continuous output');
title('Plant Output vs Disturbance');

figure,specCale(recovered.signals.values,1/Tu,50000);
hold on
specCale(offset.signals.values,1/Tu,50000);
legend('recovered singals','offset signals');
title('Spectrul Analysis for signals before and after Q filter');

w = [1 100 1000 10000 20000 30000 40000 disturbance.freq1-20000:disturbance.freq1+20000 disturbance.freq1+10000:10:2*Nyqst*2*pi];
figure,bode(filterGain*series(Qfilter,Tz),Qfilter,Tss,Tz,w);
legend('Q*PdL','Q','Tss','Tz');
title('Bode plot for PdL*Q, Q filter and continuous plant');
[mag_QT,phase_QT] = bode(filterGain*series(Qfilter,Tz),disturbance.freq1);
[mag_Q,phase_Q] = bode(Qfilter,disturbance.freq1);
[mag_Tss,phase_Tss] = bode(Tss,disturbance.freq1);
[mag_Tz,phase_Tz] = bode(Tz,disturbance.freq1);

%figure,specCale(output.signals.values,1/Ts,50000);
%title('output spectrum')
%}
