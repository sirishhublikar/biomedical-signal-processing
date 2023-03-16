clc; clear; close all;
%% Load the ECG Signal and Plot it
ecg = load('ecg_lfn.dat');
fs = 1000;
l = length(ecg);
t = [1:l]/fs;
plot(t,ecg)
xlabel('Time'); ylabel('Amplitude'); title('Noisy ECG Signal');

%% Design Derivative based filter with a zero and no poles to remove low frequency artifacts (baseline wandering)
% zeros location = cos(w) (+/-) j*sin(w); 
w = (2*pi*0)/fs;
z = [exp(1i*w) exp(-1i*w)]';
p = 0;
[b,a] = zp2tf(z,p,1); 
t1 = tf(b,a,(1/fs));
pzplot(t1) % Pole-Zero Plot

[m,f] = freqz(b,a,l,fs);
subplot(211); plot(f,abs(m)); % Magnitude Spectrum
xlabel('Freq(Hz)'); ylabel('Amplitude'); title('Magnitude Spectrum');
[p,f] = phasez(b,a,l,fs);
subplot(212); plot(f,p); % Phase Spectrum
xlabel('Phase'); ylabel('Amplitude'); title('Phase Spectrum');

%% Apply the filter to signal 
ecgfilt = filter(b,a,ecg); % Apply filter to the signal
subplot(211); plot(t,ecg);
xlabel('Time'); ylabel('Amplitude'); title('Noisy ECG Signal');
subplot(212); plot(t,ecgfilt); 
xlabel('Time'); ylabel('Amplitude'); title('Filtered ECG Signal');

%% To get a narrow bandwidth for the notch filter, add poles near the zeros
z = [exp(1i*w) exp(-1i*w)]';
r = 0.99; % Radius of pole  
p = [r*exp(1i*w) r*exp(-1i*w)];
[b,a] = zp2tf(z,p,1); 
t1 = tf(b,a,(1/fs));
pzplot(t1) % Pole-Zero Plot

[m,f] = freqz(b,a,l,fs);
subplot(211); plot(f,abs(m)); % Magnitude Spectrum
xlabel('Freq(Hz)'); ylabel('Amplitude'); title('Magnitude Spectrum');
[p,f] = phasez(b,a,l,fs);
subplot(212); plot(f,p); % Phase Spectrum
xlabel('Phase'); ylabel('Amplitude'); title('Phase Spectrum');
%% Apply the filter to signal 
ecgfilt = filter(b,a,ecg); % Apply filter to the signal
subplot(211); plot(t,ecg);
xlabel('Time'); ylabel('Amplitude'); title('Noisy ECG Signal');
subplot(212); plot(t,ecgfilt); 
xlabel('Time'); ylabel('Amplitude'); title('Filtered ECG Signal');
