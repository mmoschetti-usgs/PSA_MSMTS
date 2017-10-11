function [FREQ,AMP] = getFFT(time,series)

Ts=time(3)-time(2);
Fs=1/Ts;
N=size(series,1);

SPECTRA=fft(series)/N;
FREQ=((1:N/2-1)/N/Ts)';
SPECTRA=SPECTRA(2:size(FREQ,1)+1);
AMP=abs(SPECTRA);

end
