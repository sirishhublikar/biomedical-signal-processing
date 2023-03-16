function [dur, HR] = pantomp(ecg)
% PANTOMP Calculate QRS duration and Heart Rate of ECG Signal using
% Pan-Tompkins Algorithm with simple thresholding
%
% dur = pantomp(ecg) gives the duration of QRS peaks
%
% [dur, hr] = pantomp(ecg) gives duration and heart rate of ECG signal

fs = 200;
ecg = ecg - mean(ecg); % Normalizing
ecg = ecg/(max(abs(ecg))); % Normalizing
l = length(ecg);
t = [1:l]/fs;

% Band-Pass Filter

% Low-Pass Filter
b1 = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
a1 = [1 -2 1] * 32;

% Create digital filter (Direct Form II)
lp = dfilt.df2(b1,a1);

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


% Derivative Operator
b3 = [2 1 0 -1 -2];
a3 = [1 ] * 8;
df = dfilt.df2(b3,a3);

% Apply filter to input
ecg_filt3 = filter(df,ecg_filt2);
ecg_filt3 = ecg_filt3/max(abs(ecg_filt3));

% Squaring
ecg_filt4 = ecg_filt3.^2;
ecg_filt4 = ecg_filt4/max(abs(ecg_filt4));

% Moving window Integration
ecg_filt4Pad = [zeros(1,29) ecg_filt4' zeros(1,29)];
for i = 30:length(ecg_filt4Pad)-29
    ecg5 (i-29) = sum(ecg_filt4Pad(i-29:i))/30;
end
ecg5 = ecg5';
ecg5 = ecg5/max(abs(ecg5));

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


% QRS duration and Heart Rate
c = 0;
for i = 1:l-1
    if(ecg6(i) == 0 && ecg6(i+1) >=1 && ecg6(i+2) >= 1)
        c = c + 1;
    end
end
ecg7 = diff([ecg6; 0]);
x = find(ecg7>0);
y = find(ecg7<0);
z = y-x;

dur = mean(z)*(1/fs);
HR = c * 3;

end