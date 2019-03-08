function [Mk,Mk_rdc] = Mk_prd(Apara,L,k,p)
%Mk = Mk_prd(Apara,L,k)
%create Mk matrix and the corresponding reduced matirx A,B,b1,b2.
% w = -Mk_rdc.A*f + Mk_rdc.b1
% Mk_rdc.B*f=Mk_rdc.b2
%created by Hui Xiao
%p: p+1 data prediction: W(z^-L)= w_0 + w_1*z^-L + ... + w_p*z^-pL
%10/28/2016
if Apara(1)~=1
    error('the first coefficient must be 1');
end
if ~(k>0)
    error('invalid k');
end
m=length(Apara)-1;
if (p*L+k-m)<0
    error('invalid p');
end
Mk=zeros(p*L+k,p*(L+1)+k-m+1);
for i=1:(p*L+k-m)
    for j = 1:length(Apara)
        Mk(j+i-1,i)=Apara(j);
    end
end
index=(p*L+k-m+1):(p*(L+1)+k-m+1);
I=1:p*L+k;  %initial the index for swapping rows of Mk
for i=0:p
    Mk(k+i*L,index(i+1))=1;
    I(I==k+i*L)=[]; %delete the element that is equal to k+i*L
    I=[I,k+i*L];    %then add it to the end of the index vector
end
b=zeros(p*L+k,1);
for i=1:m
    b(i)=-Apara(i+1);
end
Mk_comb=[Mk,b];
Mk_comb(1:p*L+k,:)=Mk_comb(I,:); %swapping rows of Mk
Mk_rdc.A = Mk_comb(end-p:end,1:p*L+k-m);
Mk_rdc.B = Mk_comb(1:p*(L-1)+k-1,1:p*L+k-m);
Mk_rdc.b1 = Mk_comb(end-p:end,end);
Mk_rdc.b2 = Mk_comb(1:p*(L-1)+k-1,end);
%Mk_rdc;
end