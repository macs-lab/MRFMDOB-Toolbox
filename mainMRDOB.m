%Given Ts, L and disturbance frequencys, implement MR-FMDOB techniques to reject disturbance.
%Created by Hui Xiao, 10-13-2015

clear
addpath('C:\Users\Hui Xiao\Documents\MATLAB\library');
%load plant model
simTime = 25; %simulation time
compON = 2.5; %compensation turning on time
distON = 1; %disturbance turning on time
noiseAmp = 0.05;
% Controler parameters
Tu1 = 0.0004;  % plant identifation sampling time
continousSim = 0.0004/10;
Ts = 0.005; % Let the sensor sampling speed limits at 200Hz.
kp=30;
ki=0.2;
kd=0.05;
%% input arguments
L = input('L = ?, default is 2, maximum is 12: ');
if isempty(L)
    L = 2;
end
Tu = Ts/L;
if Tu<0.0004
    error('invalid L')
end
Nyquist = 1/Ts/2;
disp(['Nyquist frequency is: ',num2str(Nyquist)]);
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
            %freq(i)= (rand*0.9+1.1)*Nyquist;   % always beyond Ny
            freq(i)= (rand*2)*Nyquist;          % between 0 and 2*Ny
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
    dist.amp(i)   = 1/distN;
    %dist.amp(i)   = (0.5+0.5*rand)/distN;
end

noiseFlag = input('Do you want to turn on the noise input?(enter 1 for Yes,default is NO)');
if isempty(noiseFlag)
    noiseFlag = 0;
else
    noiseFlag = 1;
end

FIR_ON = input('FIR or IIR predictor? (enter 0 for IIR, default is FIR)');
if isempty(FIR_ON)
    FIR_ON = 1;
else
    FIR_ON = 0;
end
%% model defination
PdL = tf([0.013 0],[1 -1.9819 0.9819],Tu1); % the linear stage model, with sampling time Ts = 0.0004s
Pc = d2c(PdL);
PdL = d2d(PdL,Tu);
C_baseline = kp + ki*Tu*tf([1 0],[1 -1],Tu) + kd/Tu*tf([1 -1],[1 0],Tu); 
% ----> Is PID parameters still good for different Tu? needs to be determined in the future.

%% Forward model Youla Q design >> direct approach
[stdBP,~] = lattice_prod(freq,2,Tu);
phaF = phaseCompFilter_prod(...
    freqresp(PdL,freq*2*pi),...
    freq,Tu);
%QFM = series(stdBP,stdBP)*phaF;
QFM = stdBP*phaF;
if 0
    xbodeplot(QFM)
end
C_wFMDOB = (C_baseline+QFM)/(1-PdL*QFM);
loopShapingCompare(PdL,C_baseline,C_wFMDOB,Tu,'nosave',[],1);
title 'forward-model Youla: direct optimal approach'
figure, xbodeplot(1-PdL*QFM)
title 'forward-model Youla: direct optimal approach'

%% Calculate the pridictor parameters
Apara = Apara_prd(freq,Tu);
W = MMP(Apara,L);
[B,A] = IIR_MMP(freq,L,Ts,0.95);
FIR = tf_W(W,Ts);
IIR = tf_W(B,Ts,A);

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
---> try to manually use 'clear' command before running.
%}

%% Bode plot
for i=1:distN
    if freq(i) >= Nyquist
        f_b(i) = 2*Nyquist-freq(i);
        flag_f(i) = 1;        % is aliased frequency band
    else
        f_b(i) = freq(i);
        flag_f(i) = 0;        % is not aliased frequency band
    end
end
w = 0.1:0.1:Nyquist*2*pi;
[magI,phaseI]=bode(IIR,w);
[magF,phaseF]=bode(FIR,w);
mag_dbI = 20*log10(magI);
mag_dbF = 20*log10(magF);
figure,subplot(2,1,1)
plot(w/2/pi,mag_dbI(:));
hold on
plot(w/2/pi,mag_dbF(:));
xlim([0.1 Nyquist]);
for i=1:distN
    if flag_f(i)==1
        vline(f_b(i),'r--');
    else
        vline(f_b(i),'k--');
    end
end
xlabel('Frequency (Hz)');
ylabel('Magnitude (db)');
title('Bode plot');
legend('IIR predictor','FIR predictor')
subplot(2,1,2)
plot(w/2/pi,phaseI(:))
hold on
plot(w/2/pi,phaseF(:));
xlim([0.1 Nyquist]);
for i=1:distN
    if flag_f(i)==1
        vline(f_b(i),'r--');
    else
        vline(f_b(i),'k--');
    end
end
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
figure,bode(IIR);
hold on
bode(FIR);
legend('IIR design','FIR design')
%%
figure,plot(disturbance.time,disturbance.signals.values,'--');
hold on
stairs(measuredDist.time,measuredDist.signals.values,'k');
stairs(recoveredDist.time,recoveredDist.signals.values,'b');
legend('disturbance','measured disturbance','recovered disturbance');
title(['Disturbance recovery (L=',num2str(L),', distN=',num2str(distN),')']);

figure, stairs(fastSampledOutput.time,fastSampledOutput.signals.values);
hold on
legend('system output (fasted sampled)');
title(['system output (L=',num2str(L),', distN=',num2str(distN),')']);

settlingIndex = round((simTime-1)/Tu);
figure,specCale(disturbance.signals.values,1/Tu*10);
hold on
specCale(recoveredDist.signals.values,1/Tu);
specCale(filterOutput.signals.values,1/Tu);
legend('disturbance','recovered disturbance','filter output');
title(['Spectral analysis ((L=',num2str(L),', distN=',num2str(distN),')']);
