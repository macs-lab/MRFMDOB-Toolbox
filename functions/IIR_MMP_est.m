function d_est = IIR_MMP_est(dL,B,a)
%IIR_MMP_EST recover a fast sampled d_est from slow sampled dL using 
%IIR-MMP method.
%==========================================================================
% Copyright (c) 2019 Hui Xiao
%==========================================================================
% Created 1/22/2019
L = size(B,1) + 1;
d_est = zeros(1, length(dL)*L);
phi_B = zeros(size(B,2), 1);
phi_a = zeros(size(B,2), L-1);

for i = 1:length(d_est)
    k = mod(i-1,L);
    if k == 0
        measurement = dL((i-1)/L + 1);
        d_est(i) = measurement;
        phi_B = [measurement; phi_B(1:(end-1))];
    else
        d_est(i) = B(k,:)*phi_B - a*phi_a(:,k);
        phi_a(:,k) = [d_est(i); phi_a(1:(end-1),k)];
    end
end

end