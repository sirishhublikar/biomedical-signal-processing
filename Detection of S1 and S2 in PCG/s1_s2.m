clc; clear; close all;
%% Load the signal
pec = load('pec1.dat');

% Pre-processing of signal to remove baseline drift in ECG
z1 = 1; p1 = 0.995;
num = [1 -z1]; den = [1 -p1];
ft = dfilt.df2(num,den); % Digital filter (Direct Form II)
ecg = filter(ft,pec(:,2));

% Plot
PCG = pec(1050:6800,1);
crot = pec(1050:6800,3);
ecg = ecg(1050:6800);
ecg = ecg/max(abs(ecg));
PCG = PCG/max(abs(PCG));
fs = 1000;
l = length(ecg);
t = [1:l]/fs;
subplot(311); plot(t,PCG);
title('PCG Signal');
subplot(312); plot(t,ecg);
title('ECG Signal');
subplot(313); plot(t,crot);
title('Carotid Signal');

%% Pan-Tompkins Algorithm to detect QRS peaks
% Low-Pass Filter
b1 = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
a1 = [1 -2 1] * 32;

lp = dfilt.df2(b1,a1); % Create digital filter (Direct Form II)
% Apply filter to input
ecg_filt1 = filter(lp,ecg);
ecg_filt1 = ecg_filt1 - mean(ecg_filt1);
ecg_filt1 = ecg_filt1/max(abs(ecg_filt1));

% High-Pass Filter
b2 = [-1,zeros(1,15),32,-32,zeros(1,14),1];
a2 = [1 -1]*32;
hp = dfilt.df2(b2,a2);
% Apply filter to input
ecg_filt2 = filter(b2,a2,ecg_filt1);
ecg_filt2 = [ecg_filt2(1:40)*0.25; ecg_filt2(41:end)];
ecg_filt2 = ecg_filt2/max(abs(ecg_filt2));

% Plot the Result of BPF
% figure(2); plot(t,ecg_filt2);
% xlabel('Time'); ylabel('Amplitude'); title('Band-Pass Filtered Signal');
% Derivative Operator
b3 = [2 1 0 -1 -2];
a3 = [1 ] * 8;
df = dfilt.df2(b3,a3);
% Apply filter to input
ecg_filt3 = filter(df,ecg_filt2);
ecg_filt3 = ecg_filt3/max(abs(ecg_filt3));
figure(3); plot(t,ecg_filt3);
xlabel('Time'); ylabel('Amplitude'); title('Output of Derivative Operator');

% Squaring
ecg_filt4 = ecg_filt3.^2;
ecg_filt4 = ecg_filt4/max(abs(ecg_filt4));
figure(4); plot(t,ecg_filt4);
xlabel('Time'); ylabel('Amplitude'); title('Output of Squaring Operator');

% Moving window Integration
ecg_filt4Pad = [zeros(1,29) ecg_filt4' zeros(1,29)];
for i = 30:length(ecg_filt4Pad)-29
    ecg5 (i-29) = sum(ecg_filt4Pad(i-29:i))/30;
end
ecg5 = ecg5';
ecg5 = ecg5/max(abs(ecg5));
figure(5); plot(t,ecg5);
xlabel('Time'); ylabel('Amplitude'); title('Output of MA Integration');

% Thresholding
th = mean(ecg5);
ecg6 = zeros(l,1);
w = (ecg5>th);
ecg6(w) = 1;
x = find(diff([0 w']) == 1);
y = find(diff([w' 0]) == -1);
x = x - (6 + 16); % Cancel delay bcz of LPF and HPF
y = y - (6 + 16);

% Detect QRS peaks
for i = 1:length(x)
    [R_val(i), R_loc(i)] = max(ecg(x(i):y(i)));
    R_loc(i) = R_loc(i) - 1 + x(i);
    [Q_val(i), Q_loc(i)] = min(ecg(R_loc(i):-1:R_loc(i)-8));
    Q_loc(i) = R_loc(i) - Q_loc(i)+1;
    [S_val(i), S_loc(i)] = min(ecg(R_loc(i):R_loc(i)+10));
    S_loc(i) = R_loc(i) + S_loc(i)-1;
end

%% Lehner and Rangayyan method to detect dicrotic notch
% Low-Pass Filter
n = 8; fc = 40; 
h = fdesign.lowpass('N,F3dB',n,fc,fs);
ft_40 = design(h,'butter');
crot1 = filter(ft_40,crot);
crot1 = [crot1(34:end); crot1(l-32:l)]; % remove some unwanted samples
plot(t,crot1);

% Differentiator
c = [crot1(2) crot1(1) crot1' crot1(l) crot1(l-1)];

% Lehner Rangayyan Differentiator
for n = 3:l+2
    p(n-2) = 2*c(n-2) - c(n-1) - 2*c(n) - c(n+1) + 2*c(n+2);
end
p = p.^2;

% Moving Average
m = 164; % window size
p = [wrev(p(1,1:m)),p,wrev(p(1,l-m+1:l))];
s1 = zeros(1,l);
for i = m:l+m-1
    for j=1:m
        s1(i-m+1) = s1(i-m+1) + (p(i-j+1)*(m-j+1));
    end
end
s1(1:200) = 0; % Initial samples
plot(t,s1);

% Finding peaks in MA output
[pks, locs] = findpeaks(s1);
[pks1, locs1] = findpeaks(pks);
s3 = zeros(1,l);
s3(locs(locs1)) = pks1;
s4 = s3>0.5;
[~,a2] = find(s4==1);
k=1;
for i=1:length(a2)-1
    val = a2(i+1) - a2(i);
    if val > 180 && val < 340
        dpks(k) = a2(i+1);
        k=k+1;
    end
end

% Find D-notch within 20ms of second peak
for i=1:length(dpks)
    for j=dpks(i):-1:dpks(i)-50
        if (crot1(j+7) > crot1(j)) && (crot1(j-7) > crot1(j))
            a3(i) = dpks(i) - j;
            break;
        end
    end
end
dpks1 = dpks - round(mean(a3));

%% Plot
v = ones(1,length(dpks1));
figure; 
subplot(311);
plot(t,PCG,'g',t(R_loc),v,'r^',t(dpks),v,'ro');
legend('PCG','S1', 'S2'); 
subplot(312);
plot(t,ecg/max(ecg),'k',t(R_loc),R_val,'r^',t(S_loc),S_val,'r*', t(Q_loc), Q_val, 'ro');
legend('ECG','R', 'S', 'Q');
subplot(313);
plot(t,crot1,'g',t(dpks1),crot1(dpks1),'ro');
legend('Carotid Pulse', 'Dicrotic Notch');