clear all; close all;

audioinfo('MissingYou_13sec.mp3')	% find the audio info from the mp3 header
[x, Fs] = audioread('MissingYou_13sec.mp3'); % wavread() for the old version of MATLAB, Assume the sound is recorded at n = 0										 
											% x: the DT audio signal
											% Fs: sampling rate in Hz
											
x = x(:,1);	% we only process one of the two channels since both are quite similar. 

sound(x,Fs); % play the music, i.e., reconstruct the music.
             % Note that Values in x are assumed to be in the range -1.0 <= y <= 1.0. Values outside that range will be clipped by sound().

M = length(x); % length of the input signal
t = (0:M-1)*(1/Fs);      % Time axis

%% ---------- Problem 1. Plot the signal in time domain and frequency domain (magnitude spectrum) ----------
%% Can you find the signals of drum beats? What should the signals look like in time domain and in frequency domain?
X = fft(x);
mag_X = abs(X);

T = 1/Fs;
t_axis = t;
Npoint = length(x);
dF = Fs/Npoint;
f_axis = ((1:1:Npoint)-(Npoint+1)/2)*dF;
  
figure
plot(f_axis, mag_X,'linewidth',2);
figure
plot(t_axis, x,'linewidth',2);
