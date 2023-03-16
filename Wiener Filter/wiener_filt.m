clc; clear; close all;
%% Load and Plot the signal
ecg = load('ecg_hfn.dat');
fs = 1000;
l = length(ecg);
t = [1:l]/fs;
plot(t,ecg);
xlabel('Time'); ylabel('Amplitude'); title('Noisy ECG Signal');

%% Choose the points to approximate a linear model
ecg_samp = ecg(1:700);
figure;
hold on
plot(ecg_samp)
[x,y] = ginput(14);
plot(x,y,'*');

y = max(ecg_samp)*y;

%% Generate piecewise linear ECG cycle
linModel = [];
for i = 1:length(x)-1
    if isequal(y(i+1),y(i))
        a = y(i)*ones(1,floor(x(i+1) - x(i)));
    else
        a = y(i):(y(i+1) - y(i))/(x(i+1) - x(i)):y(i+1);
        a = a(1:end-1);
    end
    linModel = [linModel,a];
end
t1 = (1:length(linModel))/fs;
figure; plot(t1,linModel)

%% PSD of desired signal
nfft = max(256,2^nextpow2(length(linModel)));
[Pxx,F] = periodogram(linModel,[],nfft,fs);
figure;
plot(F,10*log10(Pxx));

%% PSD of noise
necg = ecg(2776:2948);
necg = necg - mean(necg);
[Pxx1, F] = periodogram(necg,[],nfft,fs);
% Taking Avg PSD of noise
necg2 = ecg(4205:4381);
necg2 = necg2 - mean(necg2);
[Pxx2, F] = periodogram(necg2,[],nfft,fs);

necg3 = ecg(2776:2948);
necg3 = necg3 - mean(necg3);
[Pxx3, F] = periodogram(necg3,[],nfft,fs);

Pxx_avg = (Pxx1 + Pxx2 + Pxx3)/3;
plot(F,10*log10(Pxx_avg))

%% Transfer function of Wiener filter
W = zeros(1,length(F));
for i = 1:length(F)
    W(i) = 1/(1+(Pxx_avg(i))/(Pxx(i)));
end
% W in time domain
Y = ifftshift(abs(ifft(W,200)));
ecgfilt = conv(ecg,Y);
t2 = (1:length(ecgfilt))/fs;
plot(t2,ecgfilt)


