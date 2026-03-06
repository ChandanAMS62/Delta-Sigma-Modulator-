# Delta-Sigma-Modulator-
# MATLAB IMPLEMENTATION OF DSM
# 1.) With inbuilt MATLAB functions

clear all;
clc;
osr = 128;
amp = [-100:2:0];
OBG = 1.5;
order = 3;
ntf = synthesizeNTF(order,osr,0,OBG,0);

%% parameters
Fs = 8192;
f1 = 11;
N = 8192;
n = (0:N-1)/Fs;
bits = 6;
nlevel = 2^bits;

%% calculate sutff

[a,g,b,c] = realizeNTF(ntf,'CIFF');
ABCD = stuffABCD(a,g,b,c,'CIFF');

%% simulateSNR

sqnr= simulateSNR(ABCD,osr,amp,0,nlevel,1/(2*osr),log2(N),0);
plot(amp, sqnr, '-b', 'LineWidth', 1.5);
hold on;
plot(amp, sqnr, '*', 'MarkerSize', 6)  
grid on;
hold off;
max(sqnr)

% input signal
sqnr_max = max(sqnr);
index = 0;
for i =1:51
    if sqnr(i) == sqnr_max
        index = i;
    end
end
msa = amp(index);
u = (10^(msa/10))*(nlevel - 1)*sin(2*pi*f1*n);

%% simulateDSM
x0 = zeros(order,1);
[v,xn,xmax,y] = simulateDSM(u,ABCD,nlevel,x0);
f = (0:N/2-1);
fB = ceil(Fs/(2*osr));
%% time domain plot
figure;
t = 0:1000; 
stairs(t, u(t+1),'g'); hold on; 
stairs(t,v(t+1),'b');  
ylabel('u, v'); 
title 'time domain input(green) and output(blue) of DSM'

%% spectrum analysis
w_hann = sum(ds_hann(N).^2)/N;
Xh = fft(v.*ds_hann(N))/(sqrt(N*sum(ds_hann(N).^2)));  % power correction
 
% % Xh = fft(v.*ds_hann(N))/(sum(ds_hann(N))/4); % amplitude correction

figure;
plot(10*log10(f),10*log(abs(Xh(1:N/2))));
grid on;
% manual power calculation
signal_power = 2*sum(abs(Xh(f1:f1+2)).^2);
total_power = 2*sum(abs(Xh(1:fB)).^2);
noise_power = total_power - signal_power;
snr_manual_calcualted = 10*log10(signal_power/noise_power)

%snr calculation using function
snr = calculateSNR(Xh(1:fB), f1)

text(0.7*(Fs/2), -300, sprintf('SNR = %4.1f dB', snr));

text(0.7*(Fs/2), -350, sprintf('NBW = %7.5f', 1.5/N));


%% filtering 
filter_order = 300;
wn = 3/(osr);
h = fir1(filter_order,wn,hann(filter_order+1));
h_norm  = h / sum(h);
v_filtered_delayed = filter(h_norm,1,v);
delay = filter_order/2;
v_filtered = v_filtered_delayed(delay+1:end);

%% plottting of input and filter output

figure;
plot(t,u(t+1),'b');
hold on;
plot(t,v_filtered(t+1),'r');
title 'time domain';
grid on;
legend('input signal (u)', ' filter output (v filtered)');

% frequency domain
u_ffth = fft(u.*ds_hann(N)) / (N/4);
v_filtered_ffth = fft(v_filtered_delayed.*ds_hann(N))/(N/4);
figure;
plot(10*log10(f),10*log10(abs(u_ffth(1:N/2))));
hold on;
plot (10*log10(f),10*log10(abs(v_filtered_ffth(1:N/2))));
grid on;
legend ('input signal', 'filtered ouput after DSM');

% manual power after filtering calculation
filtered_signal_power = 2*sum(abs(v_filtered_ffth(f1-2:f1+4)).^2);
filtered_total_power = 2*sum(abs(v_filtered_ffth(1:end/2)).^2);
filtered_noise_power = filtered_total_power - filtered_signal_power;
snr_after_filtering = 10*log10(filtered_signal_power/filtered_noise_power)

%% decimation
figure;
v_decimated = decimate(v_filtered,osr/4,'fir');
t1 = (0:length(v_decimated)-1)*(osr/4);
stairs(t, u(t+1),'g'); hold on; 
plot(t1,v_decimated(1:length(v_decimated)));
xlim ([0 length(t)-1]);
title 'time domain comparision of input and output after decimation'
figure;
v_decimated_fft = fft(v_decimated,length(v_decimated));
f_axis_new = (0:(length(v_decimated)/2)-1);
plot(10*log10(f_axis_new),10*log10(abs(u_ffth(1:(length(v_decimated)+1)/2))));
plot(10*log10(f_axis_new),10*log10(abs(v_decimated_fft(1:length(v_decimated)/2))));
title 'frequency domain comparision of input and output after decimation'

% manual power after filtering calculation
decimated_signal_power = 2*sum(abs(v_decimated_fft(f1:f1+2)).^2);
decimated_total_power = 2*sum(abs(v_decimated_fft(1:fB)).^2);
decimated_noise_power = decimated_total_power - decimated_signal_power;
snr_after_decimation = 10*log10(decimated_signal_power/decimated_noise_power)
