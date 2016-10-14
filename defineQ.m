function varargout = defineQ(freqinHz,Pdl,Ts,Bw)
% generate the Q filter, details in paper.
% freqinHz = Central frequency
% Pdl = Plant model
% Ts = Sampling time
% Bw = 3-dB bandwidth of the notch filter
% Hui Xiao

[~,stdBP,~] = lattice_filter(freqinHz,Bw,Ts);
[mag,phi] = bode(Pdl,freqinHz*2*pi);
phi = phi*pi/180;

b1 = (sin(phi)/sin(freqinHz*2*pi*Ts))/mag;
b0 = (cos(phi)-sin(phi)*cot(freqinHz*2*pi*Ts))/mag;
F = tf([b0, b1],[1 0],Ts);

varargout{1} = series(stdBP,F);

end

