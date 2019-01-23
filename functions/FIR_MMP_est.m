function d_est = FIR_MMP_est(dL, W)
%FIR_MMP_EST Recover d from sub sampled sequence dL
%==========================================================================
% Copyright (c) 2019 Hui Xiao
%==========================================================================
[L,p] = size(W);
L = L+1;
p = p+1;
d_est = zeros(1, (length(dL)-1)*L);
for i = 1:length(d_est)
    time = i-1;
    k = mod(time,L);
    if k == 0
        n = round(time/L)+1;
        d_est(i) = dL(n);
    else
        n = ceil(time/L);
        if n-p+2 < 1
            continue
        end
        d_est(i) = W(k,:)*fliplr(dL(n-p+2:n))';
    end
end

end

