function bode_W(ParaW,Ts)
%draw the bode plot of the predictor
%Created by Hui Xiao
%10/24/2016
n=length(ParaW);
den=[1,zeros(1,n-1)];
W=tf(ParaW,den,Ts);
bode(W);