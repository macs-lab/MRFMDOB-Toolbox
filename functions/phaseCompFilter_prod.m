function varargout = phaseCompFilter_prod(freqresP, CenterFreqHz, Tu, varargin)
% (freqresp(Pn_dTu,MNDist.FreqInHz*2*pi),MNDist.FreqInHz,Tu)
% phi = phase(freqresP);

if nargin<4
    nx = 0;
else
    nx = length(varargin{1}); % fixed gains
end

n = length(freqresP(:));

M = zeros(2*n+nx,2*n+nx);
for ii = 1:n
    for jj = 1:(2*n+nx)
        M(ii,jj) = cos((jj-1)*CenterFreqHz(ii)*2*pi*Tu);
    end
end
for ii = n+1:(2*n)
    for jj = 1:(2*n+nx)
        M(ii,jj) = sin((jj-1)*CenterFreqHz(ii-n)*2*pi*Tu);
    end
end

RHS = [real(freqresP(:))./(abs(freqresP(:)).^2);
    imag(freqresP(:))./(abs(freqresP(:)).^2)];

if nx>=1
    for ii = (2*n+1):(2*n+nx)
        for jj = 1:(2*n+nx)
            M(ii,jj) = (varargin{1}(ii-2*n))^(jj-1);
        end
    end
    RHS = [RHS;
        zeros(nx,1)];
end

b = M\RHS;

F = tf(b',[1 zeros(1,2*n+nx-1)],Tu);
% F = F/abs(freqresp(F,CenterFreqHz*2*pi))/abs(freqresP);
if 0
    %%
    phase(freqresp(F,CenterFreqHz*2*pi))
    phi
end
if nargout == 0
    figure,
    try
        xbodeplot(F)
    catch
        bodeplot(F)
    end
elseif nargout == 1
    varargout{1} = F;
else
end
