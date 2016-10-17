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
figure,plot(disturbance.time,disturbance.signals.values,'g--');
hold on
stairs(measuredDist.time,measuredDist.signals.values,'k');
stairs(recoveredDist.time,recoveredDist.signals.values,'b');
legend('disturbance','measured disturbance','recovered disturbance');
title('Disturbance recovery');

figure, stairs(fastSampledOutput.time,fastSampledOutput.signals.values);
hold on
legend('system output (fasted sampled)');
title('system output');