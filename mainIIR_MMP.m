% varify the IIR MMP and compare it with FIR MMP
clear
addpath('C:\Users\Hui Xiao\Documents\MATLAB\library');
%% Define the condition

Ts = 0.0008; %slow sampling time
L  = 2;
Tu = Ts/L;   %fast sampling time
fN = 1/Ts/2; %Nyquist frequency
f = [0.4*fN 1.3*fN];
% 7 frequency bands
 f  = [fN*0.3,fN*0.6,fN*1.3,fN*1.6,fN*1.8,fN*0.12,fN*1.5]; 
% frequency near Nyqyist frequency
 %f = 0.9*fN; 
% f1+f2 near 2*fN
%f = [0.41*fN 1.6*fN];  
n = length(f); % number of frequency bands
Apara = Apara_prd(f,Tu);
A1 = Apara(2);
A2 = Apara(3);
for i=1:n
    if f(i)>=fN
        f_b(i) = 2*fN-f(i);
        flag_f(i) = 1;        % is aliased frequency band
    else
        f_b(i) = f(i);
        flag_f(i) = 0;        % is not aliased frequency band
    end
end
%% FIR MMP
W = MMP(Apara,L);
FIR = tf_W(W,Ts);
%% IIR MMP
[B,a] = IIR_MMP(f,L,Ts,0.95);
IIR = tf_W(B,Ts,a);
%% Bode plot
w = 0.1:0.1:fN*2*pi;
[magI,phaseI]=bode(IIR,w);
[magF,phaseF]=bode(FIR,w);
mag_dbI = 20*log10(magI);
mag_dbF = 20*log10(magF);
figure,subplot(2,1,1)
plot(w/2/pi,mag_dbI(:));
hold on
plot(w/2/pi,mag_dbF(:));
xlim([0.1 fN]);
for i=1:n
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
xlim([0.1 fN]);
for i=1:n
    if flag_f(i)==1
        vline(f_b(i),'r--');
    else
        vline(f_b(i),'k--');
    end
end
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
[mag_fI,phase_fI]=bode(IIR,f*2*pi);
disp(['magnitude at frequency [',num2str(f),'] is: [',num2str(mag_fI(:)'),']']);
disp(['phase at frequency [',num2str(f),'] is: [',num2str(phase_fI(:)'),']']);
[mag_fF,phase_fF]=bode(FIR,f*2*pi);
disp(['magnitude at frequency [',num2str(f),'] is: [',num2str(mag_fI(:)'),']']);
disp(['phase at frequency [',num2str(f),'] is: [',num2str(phase_fI(:)'),']']);
figure,bode(IIR);
hold on
bode(FIR);
legend('IIR design','FIR design')
%% Generate sampled data

phi = pi*rand([1 n]); %phase of disturbance
amp = ones([1 n])/n;
step = 500;
dL        = zeros([1,step*20]);
timeL     = zeros([1,step*20]);
d         = zeros([1,step]);
time      = zeros([1,step]);
for i=1:(step+1)*20
    dL(i)=0;
    for j=1:n
        dL(i) = dL(i)+amp(j)*sin(2*pi*f(j)*Ts/20*i+phi(j)); %fast sampled singal
    end
    timeL(i) = Ts*i/20; 
end
for i=1:step
    d(i)=0;
    for j=1:n
        d(i) = d(i)+amp(j)*sin(2*pi*f(j)*Ts*i+phi(j)); %slow sampled singal
    end
    time(i) = Ts*i;
end
%% Verify the pridictor
d_FIR     = zeros([1,step]);
timeFIR   = zeros([1,step]);
d_IIR     = zeros([1,step]);
timeIIR   = zeros([1,step]);
error_FIR = zeros([1,step]);
error_IIR = zeros([1,step]);
phi_FIR = zeros(length(W),1);
for i=1:step
    % update phi_FIR
    for j=length(W):-1:2
        phi_FIR(j) = phi_FIR(j-1);
    end
    phi_FIR(1) = d(i);
    % do prediction
    d_FIR(i)= W*phi_FIR;
    error_FIR(i) = d_FIR(i)-dL(i*20+10);
    timeFIR(i) = Ts*i+Tu;
end

phi_B = zeros(length(B),1);
phi_a = zeros(length(a),1);
for i=1:step
    % update phi_B
    for j=length(phi_B):-1:2
        phi_B(j)=phi_B(j-1);
    end
    phi_B(1) = d(i);
    % update phi_a
    for j=length(phi_a):-1:2
        phi_a(j)=phi_a(j-1);
    end
    if i>1
        phi_a(1) = d_IIR(i-1);
    else
        phi_B(1) = 0;
    end
    % do prediction
    d_IIR(i) = B*phi_B - a*phi_a;
    error_IIR(i) = d_IIR(i)-dL(i*20+10);
    timeIIR(i) = Ts*i+Tu;
end
figure,
plot(timeL,dL)
hold on
stairs(time,d)
stairs(timeFIR,d_FIR,'ro');
stairs(timeIIR,d_IIR,'kx');
xlabel('time (sec)')
title('Prediction result (without noise)')

%plot the prediction error
figure,
stairs(error_IIR)
hold on
stairs(error_FIR);
legend('IIR preidictor','FIR predictor');
xlabel('Step')
ylabel('Prediction error');
title('Prediction error of IIR/FIR predictor');
%% Robustness to noise
phi = pi*rand([1 n]); %phase of disturbance
amp = ones([1 n])/n;
d_FIR   = zeros([1,step]);
timeFIR = zeros([1,step]);
d_IIR   = zeros([1,step]);
timeIIR = zeros([1,step]);
error_FIR = zeros([1,step]);
error_IIR = zeros([1,step]);
% Add uniforml distributed random noise to slow sampled signal.
d_n = d + 0.1*(2*rand([1 step])-ones(1,step));
figure,
plot(timeL,dL)
hold on
stairs(time,d_n)

phi_FIR = zeros(length(W),1);
for i=1:step
    % update phi_FIR
    for j=length(W):-1:2
        phi_FIR(j) = phi_FIR(j-1);
    end
    phi_FIR(1) = d_n(i);
    % do prediction
    d_FIR(i)= W*phi_FIR;
    error_FIR(i) = d_FIR(i)-dL(i*20+10);
    timeFIR(i) = Ts*i+Tu;
end

phi_B = zeros(length(B),1);
phi_a = zeros(length(a),1);
for i=1:step
    % update phi_B
    for j=length(phi_B):-1:2
        phi_B(j)=phi_B(j-1);
    end
    phi_B(1) = d_n(i);
    % update phi_a
    for j=length(phi_a):-1:2
        phi_a(j)=phi_a(j-1);
    end
    if i>1
        phi_a(1) = d_IIR(i-1);
    else
        phi_B(1) = 0;
    end
    % do prediction
    d_IIR(i) = B*phi_B - a*phi_a;
    error_IIR(i) = d_IIR(i)-dL(i*20+10);
    timeIIR(i) = Ts*i+Tu;
end

stairs(timeFIR,d_FIR,'ro');
stairs(timeIIR,d_IIR,'kx');
xlabel('time (sec)');
title('Prediction result (subject to evenly distribuded noise')

%plot the prediction error
figure,
stairs(error_IIR)
hold on
stairs(error_FIR);
legend('IIR preidictor','FIR predictor');
xlabel('Step')
ylabel('Prediction error');
title('Prediction error of IIR/FIR predictor');