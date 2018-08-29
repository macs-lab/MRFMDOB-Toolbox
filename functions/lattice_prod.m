function varargout = lattice_prod(freq_Hz,bw_Hz,Ts)
% 
% 2012-10-04
% Xu Chen
% 
n = length(freq_Hz);
if length(bw_Hz)~=n
    bw_Hz = bw_Hz*ones(n,1);
end
G = 1;
for ii = 1:n
    [G0,W,sys_A] = lattice_filter(freq_Hz(ii),bw_Hz(ii),Ts);
    G = G*G0;
end
BP = 1-G;
if nargout == 1
    figure, 
    try
        xbodeplot({G,BP})
    catch
        bodeplot(G,BP)
    end
elseif nargout == 1
    varargout{1} = BP;
elseif nargout == 2
    varargout{1} = BP;
    varargout{2} = G;
end