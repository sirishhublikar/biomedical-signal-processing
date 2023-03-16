clc; clear; close all;
%% Load and Plot the signal
ecg = load('ecg_lfn.dat');
fs = 1000;
l = length(ecg);
t = [1:l]/fs;
plot(t,ecg);
xlabel('Time'); ylabel('Amplitude'); title('Noisy ECG Signal');

%% Butterworth High-Pass Filter
n = 4; fc = 0.5;
[b,a] = butter(n,fc/(fs/2),'high');
[m,f] = freqz(b,a,l,fs);
subplot(211); plot(f,abs(m)); % Magnitude Spectrum
xlabel('Freq(Hz)'); ylabel('Amplitude'); title('Magnitude Spectrum');
[p,f] = phasez(b,a,l,fs);
subplot(212); plot(f,p); % Phase Spectrum
xlabel('Phase'); ylabel('Amplitude'); title('Phase Spectrum');

%% Apply filter to signal
ecgfilt = filter(b,a,ecg); % Apply filter to the signal
subplot(211); plot(t,ecg);
xlabel('Time'); ylabel('Amplitude'); title('Noisy ECG Signal');
subplot(212); plot(t,ecgfilt); 
xlabel('Time'); ylabel('Amplitude'); title('Filtered ECG Signal');
