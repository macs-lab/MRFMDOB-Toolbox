function [dc,d,dL] = SensorHarmonic(f, A, phi, b, dc_time, d_time, dL_time, noise_amp )
%MULTISIN Generate sensor measurements of harmonic signal.
%==========================================================================
% Copyright (c) 2019 Hui Xiao
%==========================================================================

dc = multisin(f,A,phi,b,dc_time);
d = multisin(f,A,phi,b,d_time);
dL = multisin(f,A,phi,b,dL_time);
dL = addNoise(dL,noise_amp);

    function y = multisin(f,A,phi,b,time_series)
        y = zeros(1,length(time_series));
        for i = 1:length(f)
            y = y + A(i)*sin(2*pi*f(i)*time_series + phi(i)) + b(i);
        end
    end
    function y = addNoise(y,amp)
        normNoise = amp*normrnd(zeros(size(y)),1);
        y = y + normNoise;
    end

end

