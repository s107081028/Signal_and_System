%%
%   Computer HW4:  Discrete-time processing of continuous-time audio signals
%				  
%
%
%                                   Edited by Meng-Lin Li, 05/28/2020
%									Dept. of Electrical Engineering,
%									National Tsing Hua University, Taiwan
%



%% ---------- FIR filtering by convolution ----------
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

%% ---------- Problem 2. Filter the music ----------

Fcut = 1000; % Hz, cut off frequency (you may try 4 kHz cutoff frequency for LPF and 0.5 kHz cutoff for HPF
FilterOrder = 8; % filter order
FilterOrder1 = 16; % filter order
FilterOrder2 = 32; % filter order
FilterOrder3 = 64; % filter order
FilterOrder4 = 128;% filter order
FilterOrder5 = 256; % filter order
flags.lowpass = 1; % 1: low pass filter
% low-pass and high-pass filter design using fir1()
% perform filtering using conv()
if flags.lowpass,
	h = fir1(FilterOrder, Fcut/(Fs/2));   % help or doc fir1(), frequency normalization is done by normalization with Fs/2 instead of Fs in our lectures
else
	h = fir1(FilterOrder, Fcut/(Fs/2), 'high'); % help fir1(), frequency normalization is done by normalization with Fs/2 instead of Fs in our lectures	
end	

h = h.'; % convert to column vector, because x is a column vector
y = conv(x, h,'same'); % filtering by convolution, why 'same'? remove the group delay introduced by the FIR LPF filter (think about relationship between group delay and the support change after convolution		
sound(y, Fs);
audiowrite('FilteredMusic.ogg',y,Fs);

figure
stem(h); % Plot the impulse response.
figure
freqz(h,1); % Plot the frequency response - log magnitude response and phase response


%% FYI
%% Equivalent codes to freqz(h,1)
%% The function of freqz() can be done by the approximated CTFT codes in the computer HW3 or by fft() based on the following codes
Nfft = 512; % zero-padding h up to Nfft points (i.e., sampling the frequency axis with Nfft sample points)
H = fft(h, Nfft); 
dF = Fs/Nfft;
AbsFreqAxis_forFFT = (0:Nfft-1)*dF; %absolute frequency axis for fft(), from 0 to Fs or equivalently from 0 to 2*pi for normalized anaular frequency
NormalizedAngularFreqAxis_forFFT = AbsFreqAxis_forFFT/Fs*2*pi; % normalized angular frequency
figure
subplot(2,1,1)
plot(NormalizedAngularFreqAxis_forFFT, 20*log10(abs(H)));
xlabel('Normalized Angular Frequency')
ylabel('Magnitude (dB)')
axis([0 pi -120 20])
grid on
subplot(2,1,2)
plot(NormalizedAngularFreqAxis_forFFT, unwrap(angle(H))*180/pi); % angle(): find the "principal phase" from -pi to pi, unwrap(): find the continuous phase spectrum
xlabel('Normalized Angular Frequency')
ylabel('Phase (degrees)')
axis([0 pi -2500 0])
grid on
%% End of equivalent codes

%% ---------- Problem 3 Re-sampling the music ----------
I = 3;
D = 2;
y_resampled = resample(x, I, D); % used to verify your own resampler
sound(y_resampled, Fs*I/D);
audiowrite('HigherSamplingRateMusic.ogg',y_resampled,Fs*I/D);

I = 2;
D = 9;
y_resampled = resample(x, I, D); % used to verify your own resampler
sound(y_resampled, Fs*I/D);
audiowrite('LowerSamplingRateMusic.ogg',y_resampled,Fs);

%% ---------- Problem 4 Noise the music to test your hearing sensitivity ----------
% Add noise (Gaussian white noise)
SNR = 20; % in dB
noise_std = 10^(-SNR/20); % standard deviation of the noise
y_noise = x + noise_std*randn(M,1);  % add noise 
sound(y_noise, Fs);
audiowrite('NoisyMusic.ogg',y_noise,Fs)

figure
freqz(y_noise,1);
figure
freqz(x,1);
%% ---------- Problem 5 Watermark the music with a single tone ----------
Fc = 65*10^3;   % The frequency of this cosine wave in Hz (you CAN change this value)
Mag = 0.1;  % magnitude of the cosine single tone (DON'T CHANGE this value for your hearing safty)

y_watermarked = x + Mag*cos(2*pi*Fc*[0:1:M-1]'/Fs); % Add this cosine signal to the original music. 
sound(y_watermarked,Fs);   
audiowrite('MusicWithInaudibleWaterMark.ogg',y_watermarked,Fs);

figure
stem(abs(x));
figure
stem(abs(y_watermarked));
figure
freqz(y_watermarked,1);