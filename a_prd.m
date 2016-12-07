function a = a_prd(freqinHz,Ts,alpha)
% function apara = apara_prd(freqinHz,Ts,alpha)
% freqinHz: given freqency components, freq = [f1 f2 f3 ... fn]
% Created by Hui Xiao, 12-4-2016
a = 1;
n = length(freqinHz);
for i = 1:n
    a = conv(a,[1 -2*alpha*cos(2*pi*Ts*freqinHz(i)) alpha^2]);
end
disp(['apara= [',num2str(a),']']);
end