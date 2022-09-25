%
%   Computer HW3:  Experiencing your first Fourier analysis using a computer
%	sample codes (Matlab script - Example of CTFT and ICTFT implementation)
%					
%                                   Edited by Meng-Lin Li, 05/02/2019
%									Revised by Meng-Lin Li, 05/07/2020 
%									Dept. of Electrical Engineering,
%									National Tsing Hua University, Taiwan
%

%% ---------- Part 1 ----------
%% ---------- Generate sampled cosine/discrete-time sinusoid ----------
F0 = 5; % in MHz
F05_1 = 0;
F05_2 = 25;
F05_3 = 50;
F05_4 = 75;
F05_5 = 95;
F05_6 = 100;
Fs = 100; % sampling rate/sampling frequency, in MHz or Msamples/sec
T = 1/Fs;  % time resolution, i.e., sampling interval in time domain
total_time = 1; % in us

% !!! Sampling in time 
t_axis = (0:T:total_time);  % time axis
x = cos(2*pi*F0*t_axis);  % sampled cosine/discrete time sinusoid ,time domain
Npoint = length(x);   % number of points in sampled cosine

figure
plot(t_axis, x,'linewidth',2);
hold
stem(t_axis, x,'r','linewidth',2);
xlabel('Time (\mus)');
ylabel('x(nT)');
title('Discrete time sinusoid (time domain)');

figure
stem(0:1:Npoint-1, x,'b', 'linewidth', 2);
xlabel('Time (n)')
ylabel('x[n]');
title('Discrete time sinusoid (time domain)')

%% ---------- Fourier transform - Analysis ----------
% !!! Sampling in frequency
dF = Fs/Npoint; % frequency resolution, i.e., sampling interval in frequency domain
%f_axis = (0:1:(Npoint-1))*dF;   % frequency axis (from 0 to Fs or equivalently from 0 to 2*pi for normalized angular frequency)
f_axis = ((1:1:Npoint)-(Npoint+1)/2)*dF; % frequency axis (from -Fs/2 to Fs/2 or equivalently from -pi to +pi for normalized angular frequency)
f_axis4_1 = ((1:1:Npoint*2)-(Npoint+1))*dF;
f_axis4_2 = ((1:1:Npoint*4)-(Npoint+1)*2)*dF;
X = zeros(1,length(f_axis)); % spectrum

% implementatoin of X(Fk) = summation x(nT)*exp(-j*2*pi*Fk*(nT))*T 
for iFreq = 1:length(f_axis),
    iFreq
   for iTime = 1:length(t_axis),
       X(iFreq) = X(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*T; % what if t_axis is not starting from time = 0?
   end
   X(iFreq)
end    

t_axis3 = (0:T*0.1:total_time);  % time axis
x3 = 0.1*cos(2*pi*F0*t_axis3);  % sampled cosine/discrete time sinusoid ,time domain
Npoint3 = length(x3);   % number of points in sampled cosine
dF3 = Fs/Npoint3; % frequency resolution, i.e., sampling interval in frequency domain
f_axis3 = ((1:1:Npoint3)-(Npoint3+1)/2)*dF3; % frequency axis (from -Fs/2 to Fs/2 or equivalently from -pi to +pi for normalized angular frequency)
X3 = zeros(1,length(f_axis3)); % spectrum
for iFreq3 = 1:length(f_axis3),
   for iTime3 = 1:length(t_axis3),
       X3(iFreq3) = X3(iFreq3) + x3(iTime3)*exp(-sqrt(-1)*2*pi*f_axis3(iFreq3)*t_axis3(iTime3))*T; % what if t_axis is not starting from time = 0?
   end
end  
mag_X3 = abs(X3);   % magnitude
figure
plot(f_axis3, mag_X3,'linewidth',2);
hold
stem(f_axis3, mag_X3, 'r', 'linewidth',1)
xlabel('Frequency (MHz)');
ylabel('abs(X3(F))')
title('Magnitude spectrum')
%% Part 1.7
% !!! You can compare the result with that from MATLAB fft()
X = fft(x); % spectrum of sampled cosine, frequency domain, complex; Any difference from the spectrum X obtained by our own CTFT codes?
f_axis_forFFT = Fs*(0:(length(x)/2))/length(x);; %frequency axis for fft(x), from 0 to Fs or equivalently from 0 to 2*pi for normalized anaular frequency
 
figure

title('Magnitude spectrum')

mag_X = abs(X);   % magnitude
pha_X = angle(X); % phase  

figure
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
xlabel('Frequency (MHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum')

figure
plot(f_axis, pha_X,'linewidth',2);
hold
stem(f_axis, pha_X, 'r', 'linewidth',1)
xlabel('Frequency (MHz)');
ylabel('phase(X(F))')
title('Phase spectrum')

 figure
 subplot(2,1,1)
 plot(t_axis, x,'linewidth',2);
 hold
 stem(t_axis, x,'r','linewidth',2);
 axis([0.2 0.8 -1 1]);   % to zoom in
 set(gca,'fontsize',14);
 set(gca,'linewidth',2);
 set(gca,'box','off');
 xlabel('Time (\mus)');
 title('Sampled cosine (time domain): x(nT)');
 legend('Original cosine','Sampled cosine','0');	
 legend('boxoff')
 
 subplot(2,1,2)
 plot(f_axis, mag_X,'linewidth',2);
 set(gca,'fontsize',14);
 set(gca,'linewidth',2);
 set(gca,'box','off');
 xlabel('Frequency (MHz)');
 title('Magnitude spectrum (frequency domain)')
 shg
 set(gca,'Xtick',[-50 -40 -30 -20 -10 -5 0 5 10 20 30 40 50]);   % set(gca,'Ytick', [ ]); set(gca,'XtickLabel',[...]);

%% ---------- Inverse Fourier transform - Synthesis (For your information if you're interested) ----------
xx = zeros(1, length(t_axis));
for iTime = 1:length(t_axis),
	iTime
	for iFreq = 1:length(f_axis),
		xx(iTime) = xx(iTime) + X(iFreq)*exp(sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*dF;
	end
end
% --- Mean Squared Errors ---
MSE = sum(abs(xx-x).^2)


figure
subplot(2,1,1)
plot(t_axis, x, 'r', 'linewidth',2);
hold
plot(t_axis, real(xx), 'k-.', 'linewidth',2); % imag(xx) should be close to 0
axis([0.2 0.8 -1 1]);   % to zoom in
set(gca,'fontsize',14);
set(gca,'linewidth',2);
set(gca,'box','off');
xlabel('Time (\mus)');
legend('Original cosine','Synthesized cosine','0');	% new !!!!!
legend('boxoff')

subplot(2,1,2)
plot(t_axis, abs(xx-x), 'linewidth',2);
xlabel('Time (\mus)');
title('Trashogram');
set(gca,'fontsize',14);
set(gca,'linewidth',2);

%% ---------- Part 2 ---------- 
%% 1 
load PPG % PPG: PPG signal, Fs: sampling rate in Hz

Fs = 100; % sampling rate/sampling frequency, in MHz or Msamples/sec
T = 1/Fs; % time resolution, sampling interval in time domain
t_axis = (1:T:1500*T); % time axis
%x = cos(2*pi*F0*t_axis); % sampled cosine/discrete time sinusoid ,time domain
x= PPG(1:1500);
Npoint = length(x); % number of points in sampled cosine

figure
plot((0:length(x)-1)/Fs, x, 'linewidth', 2);
hold
stem((0:length(x)-1)/Fs, x,'b', 'linewidth', 2);
xlabel('Time (n)')
ylabel('x[n]');
figure

plot((0:length(PPG)-1)/Fs, PPG); 
xlabel('Time (in sec.)')
ylabel('Amplitude (in mV)');
title('PPG signal')

dF = Fs/Npoint; % frequency resolution, i.e., sampling interval in frequency domain
%f_axis = (0:1:(Npoint-1))*dF;   % frequency axis (from 0 to Fs or equivalently from 0 to 2*pi for normalized angular frequency)
f_axis = ((1:1:Npoint)-(Npoint+1)/2)*dF; % frequency axis (from -Fs/2 to Fs/2 or equivalently from -pi to +pi for normalized angular frequency)
X = zeros(1,length(f_axis)); % spectrum

% implementatoin of X(Fk) = summation x(nT)*exp(-j*2*pi*Fk*(nT))*T 
for iFreq = 1:length(f_axis),
    iFreq
   for iTime = 1:length(t_axis),
       X(iFreq) = X(iFreq) + x(iTime)*exp(-sqrt(-1)*2*pi*f_axis(iFreq)*t_axis(iTime))*T; % what if t_axis is not starting from time = 0?
   end
   X(iFreq)
end    

mag_X = abs(X);

figure
plot(f_axis, mag_X,'linewidth',2);
hold
stem(f_axis, mag_X, 'r', 'linewidth',1)
xlabel('Frequency (MHz)');
ylabel('abs(X(F))')
title('Magnitude spectrum')

Fs = 100; % sampling rate/sampling frequency, in MHz or Msamples/sec
T = 1/Fs; % time resolution, sampling interval in time domain
t_axis = (1:T:14998*T); % time axis
%x = cos(2*pi*F0*t_axis); % sampled cosine/discrete time sinusoid ,time domain
x= PPG;
Npoint = length(x); % number of points in sampled cosine

figure
plot((0:length(PPG)-1)/Fs, x, 'linewidth', 2);
hold
stem((0:length(PPG)-1)/Fs, x,'b', 'linewidth', 2);
xlabel('Time (n)')
ylabel('x[n]');
figure

plot((0:length(PPG)-1)/Fs, PPG); 
xlabel('Time (in sec.)')
ylabel('Amplitude (in mV)');
title('PPG signal')

df = Fs/Npoint; % frequency resolution
f_axis = (0:1:(Npoint-1))*df;   % frequency axis
X = fft(x); % spectrum of sampled cosine, freqeuncy domain, complex
mag_X = abs(X);   % magnitude

figure
plot(f_axis, mag_X);
hold
stem(f_axis, mag_X,'r');
xlabel('Time (\mus)');
title('mag_X(time domain)');

%% 2
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

%% 3
% magnitude response of the comb reverberator in our lectures

% 
%   Demof of Echo generation and reverberation
%                                               Edited by Meng-Lin Li,05/02/2019
%                                               Dept. of Electrical Engineering, 
%												National Tsing Hua University, Taiwan

load handel
whos

%sound (y, Fs);
%pause

% ----- Echo -----
n = 1:length(y);
a = 0.7; % stable

sprintf('Echo/reverberation delay = 100 ms, attenuation factor a = 0.7')
tau = 100e-3;
D = floor(tau*Fs);
ye = filter(1, [1 zeros(1, D-1) -a], y); % comb reverberator
sound (ye, Fs);


ye2 = filter([-a zeros(1,D-1) 1], [1 zeros(1,D-1) -a], y); % all pass reverberator
%sound(ye2,Fs);

a = 1.1; % unstable
ye = filter(1, [1 zeros(1, D-1) -a], y);
%sound(ye, Fs);

Npoint = length(ye);
Fs = 100;
df = Fs/Npoint; % frequency resolution
f_axis = (0:1:(Npoint-1))*df; 

yy = fft(ye);
mag_yy = abs(yy);

figure
plot(f_axis, mag_yy);
xlabel('Frequency (MHz)');
ylabel('abs(yy)')
title('Magnitude spectrum')







