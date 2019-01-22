function [B,a] = IIR_MMP(f,L,Tss,alpha)
%function W = IIR_MMP(f,L,Ts)
%The IIR version of Multirate Model based predictor.
% f : a row vector that contains the signal frequency bands, in Hz
% L : Tu=Ts/L
% Tss: the slow sampling time of the signal.
% make alpha less than but close to 1, usually 0.9 will be good
% enough.
%====================================================
% Created by Hui Xiao
% 12/3/2016
dbstop if error
Tu = Tss/L;
Apara = Apara_prd(f,Tu);
a = a_prd(f,Tss,alpha);
n = length(f);
if Apara(1)~=1
    error('the first coefficient must be 1');
end
m = length(Apara)-1;

    function Mr = Mr_prd(Apara,n,L,r)
        Mr = zeros(2*n*L);
        for i=1:2*n*(L-1)
            for j=1:length(Apara)
                Mr(j+i-1,i) = Apara(j);
            end
        end
        for i=2*n*(L-1)+1:2*n*L
            Mr(r+L*(i-(2*n*(L-1)+1)),i)=1;
        end
    end
    
if(L>1)
    B = zeros(L-1,2*n);
    b = zeros(2*n*L,1);
    for i=1:2*n
        b(i) = -Apara(i+1);
    end
    for i=1:2*n
        b(i*L)=b(i*L)+a(i+1);
    end
    for r = 1:(L-1)
        Mr=Mr_prd(Apara,n,L,r);
        x = Mr\b;
        B(r,:) = x(end-2*n+1:end);
    end
elseif(L==1)
    B = 0;
else
    error('invalid L')
end
a = a(2:end);
end