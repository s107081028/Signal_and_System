clear all; close all;

audioinfo('MissingYou_13sec.mp3')	% find the audio info from the mp3 header
[x, Fs] = audioread('MissingYou_13sec.mp3'); % wavread() for the old version of MATLAB, Assume the sound is recorded at n = 0										 
											% x: the DT audio signal
											% Fs: sampling rate in Hz
											
x = x(:,1);	% we only process one of the two channels since both are quite similar. 

%sound(x,Fs); % play the music, i.e., reconstruct the music.
             % Note that Values in x are assumed to be in the range -1.0 <= y <= 1.0. Values outside that range will be clipped by sound().

M = length(x); % length of the input signal
t = (0:M-1)*(1/Fs);      % Time axis

Fc = 60*10^3;   % The frequency of this cosine wave in Hz (you CAN change this value)
Mag = 0.1;  % magnitude of the cosine single tone (DON'T CHANGE this value for your hearing safty)

y_watermarked = x + Mag*cos(2*pi*Fc*[0:1:M-1]'/Fs); % Add this cosine signal to the original music. 
sound(y_watermarked,Fs);   
audiowrite('MusicWithInaudibleWaterMark.ogg',y_watermarked,Fs);
%figure
%stem(abs(x));
%figure
%stem(abs(y_watermarked));
%figure
%freqz(y_watermarked,1);