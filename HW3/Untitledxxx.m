load PPG % PPG: PPG signal, Fs: sampling rate in Hz

Fs = 100; % sampling rate/sampling frequency, in MHz or Msamples/sec
T = 1/Fs; % time resolution, sampling interval in time domain
t_axis = (1:T:14998*T); % time axis
%x = cos(2*pi*F0*t_axis); % sampled cosine/discrete time sinusoid ,time domain
x = PPG;
Npoint = length(x); % number of points in sampled cosine

UnitImpulse = [1 zeros(1,59)];
z = UnitImpulse;
h = [repmat(z,1,249) zeros(1,57)]; % impulse response of the noise remover
g = x + fliplr(h);
PPG_NoiseSuppressed = PPG + fliplr(g);

figure
plot((0:length(PPG)-1), PPG_NoiseSuppressed);
xlabel('Time (n)')
ylabel('x[n]');

PPG_NoiseSuppressedfft = fft(PPG_NoiseSuppressed);
f_axis_forFFT = ((0:1:(Npoint-1))); %frequency axis for fft(x), from 0 to Fs or equivalently from 0 to 2*pi for normalized anaular frequency
Mag_X = abs(PPG_NoiseSuppressedfft);

figure
plot(f_axis_forFFT, Mag_X);
xlabel('Time (\mus)');
title('Mag_X');
