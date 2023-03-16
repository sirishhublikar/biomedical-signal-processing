clc; clear; close all;
%% Load and Plot the signal
ecg = load('ecg_hfn.dat');
fs = 1000;
l = length(ecg);
t = [1:l]/fs;
plot(t,ecg);
xlabel('Time'); ylabel('Amplitude'); title('Noisy ECG Signal');

%% Cropping the template
start = 750; % Manually Chosen by visual inspection
tempL = 801; % Manually Chosen by visual inspection
temp = ecg(start:(start+tempL-1));
t1 = (1:tempL)/fs;
plot(t1,temp)
title('Template');
%% Cross-correlation between Template and input
n = floor(tempL/2);
ecg1 = [zeros(n,1); ecg; zeros(n,1)];
for i = 1:length(ecg)
    prod = temp.*ecg1(i:i+tempL-1);
    % Normalized cross-correlation
    cc(i) = sum(prod)/(norm(temp)*norm(ecg1(i:i+tempL-1)));
end
subplot(211); plot(t,ecg);
xlabel('Time'); ylabel('Amplitude'); title('Noisy ECG Signal');
subplot(212); plot(t,cc);
xlabel('Time'); ylabel('Cross-correlation Coeff.'); title('Cross-correlation');

%% Synchronized Averaging
th = 0.8;
idx = find(cc > th);
for i = 1:length(idx)-1
    if ((idx(i)-n)>=1) && ((idx(i)+n)<=length(ecg))
        sync_avg(i,:) = ecg(idx(i)-n:idx(i)+n);
    end
end
sync_avg = mean(sync_avg);
close; plot(sync_avg)
title('Synchronized Averaging Output');
