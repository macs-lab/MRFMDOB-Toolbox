function TF = tf_W(ParaW,Ts)
%return the transfer function of MMP
%Created by Hui Xiao
%12/1/2016
n=length(ParaW);
den=[1,zeros(1,n-1)];
TF=tf(ParaW,den,Ts);
end