function Apara = Apara_prd(freqinHz,Ts)
% function Apara = Apara_prd(freq,Tu)
% Creat parameters of disturbance model.
% freqinHz: given freqency components, freq = [f1 f2 f3 ... fn]
% Ts: sampling period of the digital signal
% Created by Hui Xiao, 10-12-2016
Apara = 1;
n = length(freqinHz);
for i = 1:n
    Apara = conv(Apara,[1 -2*cos(2*pi*Ts*freqinHz(i)) 1]);
end
disp(['Apara= [',num2str(Apara),']']);
end