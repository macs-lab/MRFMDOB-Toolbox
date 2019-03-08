function [B,a] = IIR_MMP(f,L,Tss,alpha,pre_step)
%function W = IIR_MMP(f,L,Ts)
%The IIR version of Multirate Model based predictor.
% f : a row vector that contains the signal frequency bands, in Hz
% L : Tu=Ts/L
% Tss: the slow sampling time of the signal.
% make alpha less than but close to 1, usually 0.9 will be good
% enough.
%==========================================================================
% Copyright (c) 2019 Hui Xiao
%==========================================================================
% Created 12/3/2016
if(nargin == 4)
    pre_step = 0;
elseif(nargin < 4 || nargin > 5)
    error('number of argument is incorrect');
end
    
Tu = Tss/L;
Apara = Apara_prd(f,Tu);
a = a_prd(f,Tss,alpha);
n = length(f);
if Apara(1)~=1
    error('the first coefficient must be 1');
end
m = length(Apara)-1;

function Mr = Mr_prd(Apara,n,L,r)
    dim = max([2*n*L, 2*n*L-L+r]);
    Mr = zeros(dim);
    % Mr = [Mr*, e]
    % Mr* has dimention dim x (dim - 2n)
    
    for i = 1:(dim-2*n)
        for j=1:length(Apara)
            Mr(j+i-1,i) = Apara(j);
        end
    end
    ind = 0;
    for i=((dim-2*n)+1):dim
        Mr(r+ind,i)=1;
        ind = ind + L;
    end
end
    
if(L>1)
    B = zeros(L-1+pre_step,2*n);
    for r = 1:(L-1)+pre_step
        dim_now = max([2*n*L, 2*n*L-L+r]);
        b = zeros(dim_now,1);
        for i=1:2*n
            b(i) = -Apara(i+1);
        end
        for i=1:2*n
            b(i*L)=b(i*L)+a(i+1);
        end
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