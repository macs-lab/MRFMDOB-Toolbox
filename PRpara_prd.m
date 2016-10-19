function W = PRpara_prd(freqinHz,Tu,L)
% function W = PRpara_prd(freq,Tu,L)
% Creat predictor parameters by given frequency components.
% freqinHz: given freqency components, freq = [f1 f2 f3 ... fn]
% Tu: fast sampling speed in MR-FMDOB structure
% L = Ts/Tu, Ts is the slow sampling speed.
% Created by Hui Xiao, 10-12-2016
Apara = 1;
n = length(freqinHz);
for i = 1:n
    Apara = conv(Apara,[1 -2*cos(2*pi*Tu*freqinHz(i)) 1]);
end
W = MMP(Apara,L);
end









































