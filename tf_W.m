function TF = tf_W(ParaW,Ts,apara)
%function TF = tf_W(ParaW,Ts,apara)
%return the transfer function of MMP.
%Created by Hui Xiao
%12/1/2016
ParaW = ParaW(1,:);
if nargin == 2
    den = [1,zeros(1,length(ParaW)-1)];
    TF = tf(ParaW,den,Ts);
elseif nargin == 3
    num = [ParaW,0];
    den = [1,apara];
    TF = tf(num,den,Ts);
else
    error('nargin error')
end
end