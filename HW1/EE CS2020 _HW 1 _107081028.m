%%
%   Computer HW1: Experiencing sinusoidal signals and sampling
%	sample codes (Matlab script)
%					
%
%
%                                   Edited by Meng-Lin Li, 02/27/2019
%									Modified by Meng-Lin Li, 03/12/2020

%% -----------------------------------------------------
%% ---------- Codes for Problems 1 and 2 ----------
%% ------------------------------------------------------
% ---------- Generate sampled cosine/discrete-time sinusoid ----------
f0 = 880; % frequency in Hz
total_time = 5; % in sec.

% !!! Sampling in time 
fsRatio = 20;
fs = f0*fsRatio; % sampling sampling rate in Hz
T = 1/fs;  % sampling interval in time
A = 1; % magnitude
phi = pi/8; % phase

t = (0:T:total_time);  % time axis
% x(t) = A*cos(2*pi*f0*t+phi) where A: magnitude (>= 0), f0: 1/(fundamental period), phi: phase
x_CT = A*cos(2*pi*f0*(0:1/(f0*100):total_time)+phi); % x(t), sampling rate is high enough so that x_CT is a good approximation to the CT signal 
%x_CT2 = A*cos(2*pi*f0*(0 + (1/f0) : 1/(f0*100) : total_time + (1/f0))+phi);
x_DT = A*cos(2*pi*f0*t+phi);  % x[n] = x(nT), sampled cosine/discrete time sinusoid
%Npoint = length(x_DT);   % number of points in sampled cosine
sound(x_DT,fs); % play the signal by the sound card, type "help sound" (without the double quote) under MATLAB console to see the usage of sound()

figure
% type "help plot" (without the double quote) under MATLAB console to see the usage of plot().
% type "help stem" to see the usage of stem().
%plot((0:1/(f0*100):total_time), x_CT,'-o', 'linewidth', 2);
plot((0:1/(f0*100):total_time), x_CT,'-', 'linewidth', 2); % CT signal
%plot((0:1/(f0*100):total_time), x_CT2,'-*', 'linewidth', 2); % CT signal
hold on
stem(t, x_DT,'r', 'linewidth', 2); % DT signal
plot(t, x_DT,'r', 'linewidth', 2); % connet the dots, i.e., connect DT x[n]
xlabel('Time (sec.)');
ylabel('x(nT)');
title('Discrete time sinusoid (time domain)');
axis([0 1/f0*5 -A A]); % only observe the signal from time 0 to time 1/f0*3.  once you remove this line, you can see the whole sampled signal (from 0 to 5 sec.)
legend('x(t)', 'x[n]', 'connected x[n]')

% ---------- Problem 1 ----------
% (a) proof x_CT is a periodic signal with a frequency of 1/f0 mathematically and via MATLAB graphic illustration
% (b) Change A from 0.25 to 0.5, to 1, and tell the changes in the signal you observe and hear
% (c) Change f0 from 220 to 440, to 880, and tell the changes in the signal you observe and hear
% (d) Change phi from 0 to pi/8, to pi/4, to pi/2, to pi, to 3*pi/2, to 2*pi and tell the changes in the signal you observe and hear, what type of transformation of the independent variable? justify your answer.


% ---------- Problem 2 ----------
% (a) change fsRatio from 20 down to 1.2 (at least try fsRatio = 20, 4, 2.5, 2.2, 2, 1.8, 1.4, and 1.2), and tell any differences among the sounds or any differnces among the DT signals
% (b) following (a), please tell, what feature of the signal, magnitude, frequency, phase
% is changed after the sampling so that you hear the incorrect sound
% (i.e.,the sound is the not same as the sound of x(t))
% Please tell once you lower down the fs, will the oscillation rate of the origial signal be kept?
% To kept the osillation rate, what should be the smallest fs
% (c) Modify the code sound(x_DT,fs) at Line 29 to sound(x_DT,2*fs) and sound(x_DT,fs/2), respectively. Tell the changes and why.

% ---------- for fun, following 2(c) ----------
% Why is Fs chosen as 44100 Hz?
[x, Fs] = audioread('sister_12sec.wav'); % wavread() for the old version of MATLAB, Assume the sound is recorded at n = 0
										 % if you want, you can read 'sister_32sec.wav'
										 % x: the DT audio signal
										 % Fs: sampling rate in Hz
sound(x,Fs);
pause
sound(x,2*Fs);
pause
sound(x,Fs/2);



%% -----------------------------------------------------
%% ---------- Codes for Problem 3 ----------
%% -----------------------------------------------------
N = 100; % total time in normalized time
n = (0:1:N-1);  % time axis, normalized time
k = 3; % integer, 1, 2, 3, 4, 5, 6, 7, 8, 9, ... 
f1 = 1/N*k; % frequency
m = 5; % integer, from 1 to 2, 3, 4, 5, 6, 7, 8, 9, ...
f2 = 1/N*m; % frequency

% DT complex sinusoids or DT complex exponential signals
x1 =  exp(sqrt(-1)*(2*pi*f1*n)); % type "help sqrt" under MATLAB console to see the definition of sqrt()
x2 =  exp(sqrt(-1)*(2*pi*f2*n));

% type "help sum" under MATLAB console to see the usage of sum().
% type "help conj" under MATLAB console to see the usage of conj().
y = sum(x1.*conj(x2));

% ---------- Problem 3 (pure DT signal) ----------
% (a) please write down the math eq. of y, what mathematical operation is performed by the eq.? i.e., Generally, what do we call this mathematical operation, the name of this operation?
% (b) vary k and m from 1 to 50, respectively, and tell the change in y you observe and do you find any rule?
%   You can try the following codes which varies k and m via looping and
%   shows you how to display y as a function of k and m

K = 50;
M = 50;
y = zeros(K,M);
for k = 1:1:K,
   for  m = 1:1:M,
       f1 = 1/N*k; % frequency 
       f2 = 1/N*m; % frequency
       x1 =  exp(sqrt(-1)*(2*pi*f1*n));
       x2 =  exp(sqrt(-1)*(2*pi*f2*n));
       y(k,m) = sum(x1.*conj(x2));
   end
end   
figure
imagesc(abs(y)) % display y as a function of k and m
colormap(gray)
colorbar
axis image
xlabel('m')
ylabel('k')
title('y')


