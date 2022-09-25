Fs = 100; % Sampling rate in Hz
tau = 10; % Time duration of the linear FM signal in sec
B = 10; % in Hz
beta = B/tau;
A = 1;
t = 0:1/Fs:tau; % time axis
x = A*sin(pi*beta*t.^2); % the linear FM signal

%% --- Hamming windowed linear FM signal, the transmitted sonar signal from the parking sonar
x = hamming(length(x)).'.*x; 

figure
plot(t, x)
xlabel('Time (sec)');
ylabel('Amplitude')
title('Hamming windowed linear FM signal')

%% --- Perfect echo signal without noise, served as ground truth
y_woNoise = [zeros(1, 1300) x zeros(1,1000)];

% Estimate the echo time from the envelope (what is envelope, see the definition in slide 37, in Topic1_SignalsAndSystems_Part2_HandWriting0326_2020.pdf)
% You can follow the codes here to find out the echo time of noisy y
Envelope = abs(hilbert(y_woNoise)); % Envelope detection, hilbert() converts cosine or sine into its complex exponential counterpart (cos(theta) to exp(jtheta))
[PeakValue, EchoTimeIndex] = max(Envelope); % here we use the peak time as the echo time
EchoTime = EchoTimeIndex*(1/Fs); 

figure
plot( (0:(length(y_woNoise)-1))*(1/Fs), y_woNoise)
xlabel('Time (sec)');
ylabel('Amplitude')
title('Perfect echo signal with its envelope')
hold
plot( (0:(length(y_woNoise)-1))*(1/Fs), Envelope, 'r')

%% ---- Echo signal contaminated by noise
y = y_woNoise + randn(1, length(y_woNoise))*0.5; % randn() generate random numbers with zero-mean normal distribution, i.e., Gaussian distribution

figure
plot( (0:(length(y_woNoise)-1))*(1/Fs), y)
xlabel('Time (sec)');
ylabel('Amplitude')
title('Noisy echo signal, Can you find the echo time using your eyes?')

h = [zeros(1, 1300) x zeros(1,1000)]; % impulse response of the noise remover
g = y + fliplr(h);
y_NoiseSuppressed = y + fliplr(g); % noise suppressed by the noise remover, help conv(), what SHAPE do you use? 'full', 'same', or 'valid', for MATLAB 2018 or up. Give SHAPE a try and see the difference, and Is this an issue we've mentioned in our lecture? (see support and length change after convolution) 

%Envelope = abs(hilbert(y_NoiseSuppressed)); 
%[PeakValue, EchoTimeIndex] = max(Envelope); % here we use the peak time as the echo time
%EchoTime = EchoTimeIndex*(1/Fs) % should be close to the one you got for perfect echo signal