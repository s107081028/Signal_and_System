%%
%   Computer HW2: Reverberator implementation and your first experience of system design
%	sample codes (Matlab script)
%					
%
%
%                                   Edited by Meng-Lin Li, 03/28/2019
%									Modified by Meng-Lin Li, 04/09/2020

clear all; close all;

%% ----------
%% ---------- Part 1 Reverberator implementation
%% ----------
% ---------- (a) ----------

a = 0.7; % attenuation coef.
D = 5; % digital time delay
UnitImpulse = [1 zeros(1,100)]; % create DT unit impulse, starting from n=0;
x = UnitImpulse;
y = filter([-a zeros(1, D-1) 1], [1 zeros(1, D-1) -a], x); % help filter, do you know the usage of filter() now? filter([all the bk in order], [all the ak in order], x) where ak and bk are the coef. of the LCCDE in slide 21, in Topic2_LTISystems_Part2_wNote.pdf
n = 0:(length(y)-1); % check length of y[n]
figure
stem(n, y); % is it the same as that in your derivation?
xlabel('n')
title('Impulse response')

% ---------- (b) ----------
% Please comment whether this reverberator can be potentially implemented in real time and under what condition of a and D this reverberator is stable. 
%comparison of a
a = 0.5; % attenuation coef.
b = 1;
c = 1.05;
D = 5; % digital time delay
UnitImpulse = [1 zeros(1, 100)]; % create DT unit impulse, starting from n=0;
x = UnitImpulse;
y = filter([-a zeros(1, D-1) 1], [1 zeros(1, D-1) -a], x); % help filter, do you know the usage of filter() now? filter([all the bk in order], [all the ak in order], x) where ak and bk are the coef. of the LCCDE in slide 21, in Topic2_LTISystems_Part2_wNote.pdf
z = filter([-b zeros(1, D-1) 1], [1 zeros(1, D-1) -b], x); % help filter, do you know the usage of filter() now? filter([all the bk in order], [all the ak in order], x) where ak and bk are the coef. of the LCCDE in slide 21, in Topic2_LTISystems_Part2_wNote.pdf
w = filter([-c zeros(1, D-1) 1], [1 zeros(1, D-1) -c], x); % help filter, do you know the usage of filter() now? filter([all the bk in order], [all the ak in order], x) where ak and bk are the coef. of the LCCDE in slide 21, in Topic2_LTISystems_Part2_wNote.pdf

n = 0:(length(y)-1); % check length of y[n]
figure
stem(n, y); % is it the same as that in your derivation?
hold on;
stem(n, z);
hold on;
stem(n, w);
xlabel('n')
title('Impulse response')
legend("a = 0.5", "a = 1", "a = 1.05")

%comparison of d
a = 0.5; % attenuation coef.
b = 1;
c = 3;
D = 5; % digital time delay
UnitImpulse = [1 zeros(1, 100)]; % create DT unit impulse, starting from n=0;
x = UnitImpulse;
y = filter([-a zeros(1, b-1) 1], [1 zeros(1, b-1) -a], x); % help filter, do you know the usage of filter() now? filter([all the bk in order], [all the ak in order], x) where ak and bk are the coef. of the LCCDE in slide 21, in Topic2_LTISystems_Part2_wNote.pdf
z = filter([-a zeros(1, c-1) 1], [1 zeros(1, c-1) -a], x); % help filter, do you know the usage of filter() now? filter([all the bk in order], [all the ak in order], x) where ak and bk are the coef. of the LCCDE in slide 21, in Topic2_LTISystems_Part2_wNote.pdf
w = filter([-a zeros(1, D-1) 1], [1 zeros(1, D-1) -a], x); % help filter, do you know the usage of filter() now? filter([all the bk in order], [all the ak in order], x) where ak and bk are the coef. of the LCCDE in slide 21, in Topic2_LTISystems_Part2_wNote.pdf

n = 0:(length(y)-1); % check length of y[n]
figure
stem(n, y); % is it the same as that in your derivation?
hold on;
stem(n, z);
hold on;
stem(n, w);
xlabel('n')
title('Impulse response')
legend("d = 1", "d = 3", "d = 5")


% ---------- (c) ----------
% FIR system design and implementation
[x, Fs] = audioread('Halleluyah.wav'); % wavread() for the old version of MATLAB, Assume the sound is recorded at n = 0
a = 0.7;
D = 819;
%sound(x, Fs);
% FIR approximation
h = [filter([-a zeros(1, D-1) 1], [1 zeros(1, D-1) -a], [1 zeros(1, 100)])]; % impulse response of the FIR-approximated reverberator
y = conv(x, h); % help conv(), what SHAPE do you use? 'full', 'same', or 'valid', for MATLAB 2018 or up. Give SHAPE a try and see the difference, and Is this an issue we've mentioned in our lecture? (see support and length change after convolution)
             % what initial condition is assumed when you use conv()?
             
sound(y, Fs); % play the created sound with reverberation 
audiowrite('Halleluyah_FIRecho.wav', y, Fs); % wavwrite() for the old version of MATLAB

figure
plot((0:1:73112), x, '-', 'linewidth', 2);
hold on;
plot((0:1:73212), y, '-', 'linewidth', 2);
legend("sample", "FIR");

% ---------- (d) ----------
% IIR system implementation and initial condition
[x, Fs] = audioread('Halleluyah.wav'); % wavread() for the old version of MATLAB
a = 0.7;
D = 819;
sound(x, Fs);
% IIR system implementation (see slide 48, in Topic2_LTISystems_Part2_wNote.pdf), verify your implementation by MATLAB function filter()
y(1, 1) = (-1) * x(1,1);
for n = 1:(1.4*length(x))
    if(n >= length(x)) 
        y(1, n + 1) = a*y(1, n+1 - D);
    else if(n >= D)
            y(1, n + 1) = a*y(1, n+1 - D) - a*x( n + 1,1) + x(n+1 - D, 1);
        end
    end
end
sound(y, Fs); % play the created sound with reverberation 
audiowrite('Halleluyah_IIRecho.wav', y, Fs); % wavwrite() for the old version of MATLAB

figure
plot((1:1:124293), y, '-', 'linewidth', 2);
hold on;
plot((0:1:73112), x, '-', 'linewidth', 2);
legend("IIR", "sample");


% ---------- (e) ----------
% Tuning the attenuation coef. so that you're going to have an unstable system
% What does the output of an unstable reverberator shoulds like?
[x, Fs] = audioread('Halleluyah.wav'); % wavread() for the old version of MATLAB
a = 1.01;
D = 819;
sound(x, Fs);
% IIR system implementation (see slide 48, in Topic2_LTISystems_Part2_wNote.pdf), verify your implementation by MATLAB function filter()
y(1, 1) = (-1) * x(1,1);
for n = 1:(1.4*length(x))
    if(n >= length(x)) 
        y(1, n + 1) = a*y(1, n+1 - D);
    else if(n >= D)
            y(1, n + 1) = a*y(1, n+1 - D) - a*x( n + 1,1) + x(n+1 - D, 1);
        end
    end
end
sound(y, Fs); % play the created sound with reverberation 
audiowrite('Halleluyah_IIRecho.wav', y, Fs); % wavwrite() for the old version of MATLAB

figure
plot((1:1:102359), y, '-', 'linewidth', 2);
hold on;
plot((0:1:73112), x, '-', 'linewidth', 2);
legend("IIR", "sample");

% ---------- (f) ----------
% Tuning the initial condition
[x, Fs] = audioread('Halleluyah.wav'); % wavread() for the old version of MATLAB
a = 0.7;
D = 819;
sound(x, Fs);
% IIR system implementation (see slide 48, in Topic2_LTISystems_Part2_wNote.pdf), verify your implementation by MATLAB function filter()
for n = 1:(1.4*length(x))
    if(n >= length(x)) 
        y(1, n + 1) = a*y(1, n+1 - D);
        w(1, n + 1) = a*w(1, n+1 - D);
    else if(n < 819)
            y(1, n + 1) = a*1 - a*x(n + 1, 1) + 1;
            w(1, n + 1) = a*rand() - a*x(n + 1, 1) + rand();
    else if(n >= D)
            y(1, n + 1) = a*y(1, n + 1 - D) - a*x( n + 1, 1) + x(n + 1 - D, 1);
            w(1, n + 1) = a*w(1, n + 1 - D) - a*x( n + 1, 1) + x(n + 1 - D, 1);
        end
        end
    end
end
sound(y, Fs); % play the created sound with reverberation 
%audiowrite('Halleluyah_Initialecho.wav', y, Fs); % wavwrite() for the old version of MATLAB

figure
plot((1:1:102359), y, '-', 'linewidth', 2);
hold on;
plot((1:1:102359), w, '-', 'linewidth', 2);
hold on;
plot((0:1:73112), x, '-', 'linewidth', 2);
legend("= 1", "random", "sample");

% ---------- (g) ----------
% Verify the system output = zero input response + zero-state response
[x, Fs] = audioread('Halleluyah.wav'); % wavread() for the old version of MATLAB
a = 0.7;
D = 819;
sound(x, Fs);
% IIR system implementation (see slide 48, in Topic2_LTISystems_Part2_wNote.pdf), verify your implementation by MATLAB function filter()
for n = 1:(1.4*length(x))
    if(n >= length(x)) 
        y(1, n + 1) = a*y(1, n+1 - D);
        w(1, n + 1) = a*w(1, n+1 - D);
        z(1, n + 1) = a*z(1, n+1 - D);
    else if(n < 819)
            y(1, n + 1) = a*1 - a*x(n + 1, 1) + 1;
            w(1, n + 1) = y(1, n + 1);
    else if(n >= D)
            y(1, n + 1) = a*y(1, n + 1 - D) - a*x(n + 1, 1) + x(n + 1 - D, 1);
            w(1, n + 1) = a*w(1, n + 1 - D);
            z(1, n + 1) = a*z(1, n + 1 - D) - a*x(n + 1, 1) + x(n + 1 - D, 1);
        end
        end
    end
end
sound(y, Fs); % play the created sound with reverberation 
%audiowrite('Halleluyah_Initialecho.wav', y, Fs); % wavwrite() for the old version of MATLAB

figure
plot((1:1:102359), y, '-', 'linewidth', 2);
hold on;
plot((1:1:102359), z, '-', 'linewidth', 2);
hold on;
plot((1:1:102359), w, '-', 'linewidth', 2);
legend("origin", "ZIR", "ZSR");


%% ----------
%% ---------- Part 2 First experience of system design 
%% ----------
% Linear FM/Chirp signal x (if you're interested, you can google "Linear FM signal" or "Chirp signal")
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


%% ---- Design your noise reduction system here
h = [zeros(1, 1300) x zeros(1,1000)]; % impulse response of the noise remover
g = y + fliplr(h);
y_NoiseSuppressed = y + fliplr(g); % noise suppressed by the noise remover, help conv(), what SHAPE do you use? 'full', 'same', or 'valid', for MATLAB 2018 or up. Give SHAPE a try and see the difference, and Is this an issue we've mentioned in our lecture? (see support and length change after convolution) 

%Envelope = abs(hilbert(y_NoiseSuppressed)); 
%[PeakValue, EchoTimeIndex] = max(Envelope); % here we use the peak time as the echo time
%EchoTime = EchoTimeIndex*(1/Fs) % should be close to the one you got for perfect echo signal


