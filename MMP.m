function W = MMP(Apara,L)
%function W = MMP(Aparameters,L)
%Implement the procedure to produce the multirate model based predictor.
%Reference 'Multirate Forward-model Disturbance Observer for Feedback
%Regulation beyond Nyquist Frequency, Xu chen and Hui Xiao, systems &
%control letters, 2016'
%===============================================
%Apara: model coefficients of disturbance. 
%   Apara=[1 a1 a2 ... an] such that A(z-1)= 1+a1*z-1+...+an*z-n, and 
%   A(z-1)*d[n]=0
%L is the division parameter.
%W is a matrix containing the predictor parameters.
%  W(i,:)=[w_i_0, w_i_1, ..., w_i_m-1], i=1,2,3...,L-1
%================================================
%Created by Hui Xiao, 10-12-2016
dbstop if error
if Apara(1)~=1
    error('the first coefficient must be 1');
end
m = length(Apara)-1;
W = zeros(L-1,m);
for k = 1:(L-1)
    Mk = zeros(L*(m-1)+k);
    for i = 1:(L*(m-1)-m+k)
        for j = i:i+m
            Mk(j,i)=Apara(j-i+1);
        end
    end
    for i = ((L*(m-1)-m+k)+1):(L*(m-1)+k)
        Mk(k+(i-((L*(m-1)-m+k)+1))*L,i) = 1;
    end
    a = -1*[Apara(2:end)';zeros(L*(m-1)-m+k,1)];
    if rank(Mk)~=length(Mk)
        %disp(['The disturbance model is: [',num2str(Apara),']']);
        error(['Mk is nonsingular/badly conditioned at k=',num2str(L),', can''t find a unique solution!'])
    end
    b = Mk\a;
    W(k,:)=b(end-m+1:end)';
end

end