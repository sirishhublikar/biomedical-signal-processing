clc; clear; close all;
%% Load the signal 
eeg = load('eeg1-c3.dat');
fs = 100;
l = length(eeg);
t = (1:l)/fs;
figure(1); plot(t,eeg)
xlabel('Time'); ylabel('Amplitude'); title('EEG signal with Alpha Rhythm');
%% Choose the template of alpha wave
temp = eeg(470:570);
tempL = length(temp);
figure(2); plot((470:570)/fs,temp)
xlabel('Time'); ylabel('Amplitude'); title('Template of Alpha Rhythm');
%% Cross-Correlation 
eeg1 = [zeros(50,1); eeg; zeros(50,1)];
% Normalized Cross correlation
for j = 1:length(eeg)
    prod = temp.*eeg1(j:j+tempL-1);
    cc(j) = sum(prod) / (norm(eeg1(j:j+tempL-1))*norm(temp));
end
% plot(1:length(cc),cc);

% Cross-Correlation thresholding
th = 0.5*max(cc);
bin_cc = cc > th;
cc_th = bin_cc .* cc;
figure(3)
subplot(211); plot(1:length(cc),cc);
title('Cross-Correlation of input eeg and template');
subplot(212); plot(1:length(cc_th),cc_th);
title('Thresholded cross-correlation');
%% Peak detection
n=1;
for j = 2:length(cc_th)-1
    if cc_th(j) > cc_th(j-1) && cc_th(j)> cc_th(j+1)
        peak_loc(n) = j;
        n = n+1;
    end
end
% Delay between peaks
peak_loc_delay = diff(peak_loc);

% Histogram of delay between peaks
[counts, ~] = hist(peak_loc_delay,1:50);
figure; stem(counts);
xlabel('Index'); ylabel('Counts');
title('Histogram of delay between peaks');
[~, ind] = max(counts); % index of max count of historgram
freq = 1/(ind/fs); % dominant frequency
disp(['Corresponding frequency for input is ' num2str(freq) ' Hz']);