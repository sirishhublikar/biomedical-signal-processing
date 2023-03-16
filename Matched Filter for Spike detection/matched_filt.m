clc; clear; close all;
%% Load the EEG files
fs = 100;
u = ["c3", "c4", "f3", "f4", "o1", "o2", "p3", "p4", "t3", "t4"];

for i = 1:length(u)
    data = "eeg2-" + u(i) + ".dat";
    EEG(i,:) = load(data);
    EEG(i,:) = EEG(i,:)/max(abs(EEG(i,:)));
end

l = length(EEG);
t = (1:l)/fs;
plot(t,EEG(3,:)) % Plot the EEG signal

%% Extract the template of Spike-and-wave
temp = EEG(3,60:82);

%% Create matched filter
mfilt = wrev(temp); % reverses the template
num = mfilt; den = [1 0];
mf1 = dfilt.df2(num,den); % Digital filter (Direct Form II)

%% Apply the filter to input
c = filter(mf1,EEG(3,:));
subplot(211); plot(t,EEG(3,:));
title('EEG Signal');
subplot(212); plot(t,c);
title('Output of matched filter');
