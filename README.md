# Delta-Sigma-Modulator-
# MATLAB IMPLEMENTATION OF DSM
# 1.) Without DAC capacitor mismatch MATLAB functions

clear all;
clc;

%% parameters
Fs = 16384;
f1 = 13;
N = 16384;
n = (0:2*N-1)/Fs;
bits = 1;
nlevel = 2^bits;
f = (0:N/2-1);
t = 0:16384; 
osr = 128;
fB = Fs/(2*osr);
amp = [-100:2:0];
OBG = 1.5;
order = 3;
ntf = synthesizeNTF(order,osr,0,OBG,0);

v = zeros(1,2*N);
v_dac = zeros(1,2*N);


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

y1 = zeros(1,2*N);
y2 = zeros(1,2*N);
y3 = zeros(1,2*N);
x = 0;
y = zeros(1,2*N);
y_min = min(u);
y_max = max(u);
%% DSM

for i = 2:2*N
    % integration
    y1(i) = y1(i-1) + c(1)*x + b(2); 
    y2(i) = y2(i-1) + c(2)*y1(i-1) + b(3);
    y3(i) = y3(i-1) + c(3)*y2(i-1) ;

    y(i) = b(4)*u(i) + y1(i)*a(1) + y2(i)*a(2) + y3(i)* a(3);

    % quantizer
    
    delta =(y_max - y_min)/nlevel ;
    quantized_y = round( (y(i) - y_min)/delta )*delta + y_min ;
    
    v(i) = quantized_y;
    v_dac(i) =  quantized_y;
    x = u(i) - v_dac(i);

end

%% time domain plot
figure;
t = 1:2048; 
stairs(t, u(t+1),'g'); hold on; 
stairs(t,v(t+1),'b');  
ylabel('u, v'); 
title 'time domain input(green) and output(blue) of DSM'


%% spectrum analysis
% w_hann = sum(ds_hann(N).^2)/N;
Xh = fft(v(1:N).*ds_hann(N))/(sqrt(N*sum(ds_hann(N).^2)));  % power correction
mag_Xh = abs(Xh);
% % Xh = fft(v(1:N).*ds_hann(N))/(sum(ds_hann(N))/4); % amplitude correction

figure;
plot(10*log10(f),20*log(abs(Xh(1:N/2))));
grid on;
% manual power calculation
signal_power = 2*sum(mag_Xh((f1:f1+2)).^2);
total_power = 2*sum(abs(Xh(1:fB)).^2);
noise_power = total_power - signal_power;
snr_manual_calcualted = 10*log10(signal_power/noise_power)

%snr calculation using function
snr = calculateSNR(Xh(1:fB), f1)

text(0.7*(Fs/2), -300, sprintf('SNR = %4.1f dB', snr));

text(0.7*(Fs/2), -350, sprintf('NBW = %7.5f', 1.5/N));

%% filtering 
filter_order = 800;
wn = (Fs/(2*osr))/(Fs/2);
h = fir1(filter_order,wn,hann(filter_order+1));
h_norm  = h / sum(h);
% % h_quant_norm =fi(h_norm, 1, 16, 15);   % signed, 16-bit word, 15-bit fraction
% % v_filtered_trunc_delayed = filter(h_quant_norm,1,fi(v,1,16,15));

B = 8;                     % number of coefficient bits
scale = 2^(B);  
h_quant = round(h.*scale)/scale;   % truncated / quantized coefficients
h_quant_norm = h_quant / sum(h_quant);
v_filtered_trunc_delayed = filter(h_quant_norm,1,v);
v_filtered_trunc_delayed_trunc_scale = (v_filtered_trunc_delayed.*scale)/scale;

v_filtered_trunc_delayed_trunc = v_filtered_trunc_delayed_trunc_scale;
for i=1:length(v_filtered_trunc_delayed_trunc_scale)
    if v_filtered_trunc_delayed_trunc_scale >= (scale-1)
       v_filtered_trunc_delayed_trunc = scale-1;
    elseif v_filtered_trunc_delayed_trunc_scale <= -scale
        v_filtered_trunc_delayed_trunc = -scale;
    else
        v_filtered_trunc_delayed_trunc = v_filtered_trunc_delayed_trunc_scale;
    end
end

v_filtered_delayed = filter(h_norm,1,v);
delay = filter_order/2;
v_filtered = v_filtered_delayed(delay+1:delay+N);
v_filtered_trunc = v_filtered_trunc_delayed(delay+1:delay+N);
v_filtered_trunc_trunc = v_filtered_trunc_delayed_trunc(delay+1:delay+N);

figure;
stairs(t, v_filtered_delayed(t+1),'b'); hold on;
stairs(t, v_filtered_trunc_delayed(t+1),'r');
stairs(t, v_filtered_trunc_delayed_trunc(t+1),'g');

legend('Ideal filter','Truncated coeff filter','trucated cofficient & output both');
title('Effect of coefficient truncation');
grid on;

%% plottting of input and filter output

figure;
plot(t,u(t+1),'b');
hold on;
plot(t,v_filtered(t+1),'r');
plot(t,v_filtered_trunc_trunc(t+1),'g');
title 'time domain';
grid on;
legend('input signal (u)', ' filter output (without truncation)','filter output (with truncation)');
% N_changed = pow2(floor(log2(length(v_filtered))));
% frequency domain
u_ffth = fft(u(1:N).*ds_hann(N)) / (sqrt(N*sum(ds_hann(N).^2)));
v_filtered_ffth = fft(v_filtered(1:N).*ds_hann(N))/ (sqrt(N*sum(ds_hann(N).^2)));
v_filtered_trunc_ffth = fft(double(v_filtered_trunc(1:N).*ds_hann(N)))/ (sqrt(N*sum(ds_hann(N).^2)));
v_filtered_trunc_trunc_ffth = fft(double(v_filtered_trunc_trunc(1:N).*ds_hann(N)))/ (sqrt(N*sum(ds_hann(N).^2)));

figure;
plot(10*log10(f),10*log10(abs(u_ffth(f+1))),'r');
hold on;
plot(10*log10(f),20*log(mag_Xh(f+1)),'g');
plot (10*log10(f),10*log10(abs(v_filtered_ffth(f+1))),'b');
plot (10*log10(f),10*log10(abs(v_filtered_trunc_ffth(f+1))));
plot (10*log10(f),10*log10(abs(v_filtered_trunc_trunc_ffth(f+1))));

grid on;
legend ('input signal','DSM output without filtering' ,'filtered ouput after DSM','filtered output with trucated coefficient','truncated filter output');

% manual power after filtering calculation
filtered_signal_power = 2*sum(abs(v_filtered_ffth(f1:f1+2).^2));
filtered_total_power = 2*sum(abs(v_filtered_ffth(1:N/2)).^2);
filtered_noise_power = filtered_total_power - filtered_signal_power;
snr_after_filtering = 10*log10(filtered_signal_power/filtered_noise_power)

% power calculation for filter output with trucated filter coefficiet
filtered_signal_power_trunc = 2*sum(abs(v_filtered_trunc_ffth(f1:f1+2).^2));
filtered_total_power_trunc = 2*sum(abs(v_filtered_trunc_ffth(1:N/2)).^2);
filtered_noise_power_trunc = filtered_total_power_trunc - filtered_signal_power_trunc;
snr_after_filtering_trunc = 10*log10(filtered_signal_power_trunc/filtered_noise_power_trunc)

% power calculation for filter output with trucated filter coefficiet
filtered_signal_power_trunc_trunc = 2*sum(abs(v_filtered_trunc_trunc_ffth(f1:f1+2).^2));
filtered_total_power_trunc_trunc = 2*sum(abs(v_filtered_trunc_trunc_ffth(1:N/2)).^2);
filtered_noise_power_trunc_trunc = filtered_total_power_trunc_trunc - filtered_signal_power_trunc_trunc;
snr_after_filtering_trunc_trunc = 10*log10(filtered_signal_power_trunc_trunc/filtered_noise_power_trunc_trunc)

%% decimation
figure;
v_decimated = decimate(v_filtered(1:N),osr,'fir');
v_trunc_trunc_decimated = decimate(v_filtered_trunc_trunc(1:N),osr,'fir');

t1 = (0:osr:N-1);
stairs(t, u(t+1),'g'); hold on; 
stairs(t1,v_decimated(1:length(t1)),'r');
stairs(t1,v_trunc_trunc_decimated(1:length(t1)),'b');
legend ('input signal','decimated signal without truncation','decimated signal after truncation')

xlim ([0 length(t)-1]);
title 'time domain comparision of input and output after decimation'
figure;
v_decimated_hann = v_decimated.*ds_hann(length(v_decimated));
v_trunc_trunc_decimated_hann = v_trunc_trunc_decimated.*ds_hann(length(v_trunc_trunc_decimated));

v_decimated_fft = abs(fft(v_decimated_hann,length(v_decimated)))/(sqrt(length(v_decimated)*sum(ds_hann(length(v_decimated)).^2)));
v_trunc_trunc_decimated_fft = abs(fft(v_trunc_trunc_decimated_hann,length(v_trunc_trunc_decimated)))/(sqrt(length(v_trunc_trunc_decimated)*sum(ds_hann(length(v_trunc_trunc_decimated)).^2)));

f_axis_new = (0:(length(v_decimated_fft)/2)-1);
plot(10*log10(f_axis_new),10*log10(abs(v_decimated_fft(1:length(v_decimated_fft)/2)))); hold on;
plot(10*log10(f_axis_new),10*log10(abs(v_trunc_trunc_decimated_fft(1:length(v_decimated_fft)/2))));
legend ('decimated signal without truncation','decimated signal after truncation')
title 'frequency domain comparision of input and output after decimation'

% manual power after decimation calculation
decimated_signal_power = 2*sum(abs(v_decimated_fft(f1:f1+2)).^2);
decimated_total_power = 2*sum(abs(v_decimated_fft(1:end/2)).^2);
decimated_noise_power = decimated_total_power - decimated_signal_power;
snr_after_decimation = 10*log10(decimated_signal_power/decimated_noise_power)

% manual power after decimation calculation of truncated signal
trunc_trunc_decimated_signal_power = 2*sum(abs(v_trunc_trunc_decimated_fft(f1:f1+2)).^2);
trunc_trunc_decimated_total_power = 2*sum(abs(v_trunc_trunc_decimated_fft(1:end/2)).^2);
trunc_trunc_decimated_noise_power = trunc_trunc_decimated_total_power - trunc_trunc_decimated_signal_power;
trunc_trunc_snr_after_decimation = 10*log10(trunc_trunc_decimated_signal_power/trunc_trunc_decimated_noise_power)

# 2.) With DAC capacitor mismatch and Data Weighted Averaging
clear all;
clc;

%% parameters
Fs = 16384;
f1 = 13;
N = 16384;
n = (0:2*N-1)/Fs;
bits = 1;
nlevel = 2^bits;
f = (0:N/2-1);
t = 0:16384; 
osr = 128;
fB = Fs/(2*osr);
amp = [-100:2:0];
OBG = 1.5;
order = 4;
ntf = synthesizeNTF(order,osr,0,OBG,0);

v = zeros(1,2*N);
% v(1) = 1; %for NTF 
v_dac = zeros(1,2*N);
% v_dac(1) = 1;
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
%clock jitter
sigma_jitter = 0;
t_jittered = n*(1+ (sigma_jitter*randn)/Fs);
u = (nlevel-1)+(10^(msa/10))*(nlevel-1)*sin(2*pi*f1*t_jittered);

% u= zeros(1,length(n)); %for NTF
y1 = zeros(1,2*N);
y2 = zeros(1,2*N);
y3 = zeros(1,2*N);
x = zeros(1,2*N);
y = zeros(1,2*N);
levels=zeros(1,nlevel);
for i=1:nlevel
    levels(i)=i-1;
end   
% mismatch details
sigma=10/100;
cap_array = 2.*(1 + sigma * randn(1,nlevel-1)); 


%% DSM
pointer = 0;
for i = 2:2*N
    sum1=0;
    % integration
    y1(i) = y1(i-1) + b(1)*u(i)-c(1)*v_dac(i-1); 
    y2(i) = y2(i-1) + c(2)*y1(i-1) + b(2)*u(i);
    y3(i) = y3(i-1) + c(3)*y2(i-1) +b(3)*u(i) ;

    y(i) = b(4)*u(i) + y1(i)*a(1) + y2(i)*a(2) + y3(i)* a(3);
    
    % quantizer
    [~,idx] = min(abs(y(i) - levels));
    if y(i)-levels(idx) <0
        idx = idx-1;
    end

    for j = 1:idx
        elem_index = mod(pointer + j - 1, nlevel-1) + 1;   % circular index
    sum1 = sum1 + cap_array(elem_index);
    end
    pointer = mod(pointer + idx, nlevel-1);

    v_dac(i)=sum1;
    v(i) = idx;

end


%% time domain plot
figure;
t = 1:2048; 
stairs(t, u(t+1),'g'); hold on; 
stairs(t,v(t+1),'b');  
ylabel('u, v'); 
legend ('input signal','DSM output (v)');
title 'time domain input and output of DSM'


%% spectrum analysis
% w_hann = sum(ds_hann(N).^2)/N;
Xh = fft(v(1:N).*ds_hann(N))/(sqrt(N*sum(ds_hann(N).^2)));  % power correction
mag_Xh = abs(Xh);
% % Xh = fft(v(1:N).*ds_hann(N))/(sum(ds_hann(N))/4); % amplitude correction

figure;
semilogx(f,10*log10(abs(Xh(1:N/2))));
grid on;
% manual power calculation
signal_power = 2*sum(mag_Xh((f1:f1+2)).^2);
total_power = 2*sum(abs(Xh(4:fB)).^2);
noise_power = total_power - signal_power;
snr_manual_calcualted = 10*log10(signal_power/noise_power)

% %snr calculation using function
% snr = calculateSNR(Xh(4:fB), f1)
% 
% text(0.7*(Fs/2), -300, sprintf('SNR = %4.1f dB', snr));
% 
% text(0.7*(Fs/2), -350, sprintf('NBW = %7.5f', 1.5/N));

%% filtering 
filter_order = 800;
wn = (Fs/(2*osr))/(Fs/2);
h = fir1(filter_order,wn,hann(filter_order+1));
h_norm  = h / sum(h);
% % h_quant_norm =fi(h_norm, 1, 16, 15);   % signed, 16-bit word, 15-bit fraction
% % v_filtered_trunc_delayed = filter(h_quant_norm,1,fi(v,1,16,15));

B = 16;                     % number of coefficient bits
scale = 2^(B);  
h_quant = round(h.*scale)/scale;   % truncated / quantized coefficients
h_quant_norm = h_quant / sum(h_quant);
v_filtered_trunc_delayed = filter(h_quant_norm,1,v);
v_filtered_trunc_delayed_trunc_scale = (v_filtered_trunc_delayed.*scale)/scale;

v_filtered_trunc_delayed_trunc = v_filtered_trunc_delayed_trunc_scale;
for i=1:length(v_filtered_trunc_delayed_trunc_scale)
    if v_filtered_trunc_delayed_trunc_scale >= (scale-1)
       v_filtered_trunc_delayed_trunc = scale-1;
    elseif v_filtered_trunc_delayed_trunc_scale <= -scale
        v_filtered_trunc_delayed_trunc = -scale;
    else
        v_filtered_trunc_delayed_trunc = v_filtered_trunc_delayed_trunc_scale;
    end
end

v_filtered_delayed = filter(h_norm,1,v);
delay = filter_order/2;
v_filtered = v_filtered_delayed(delay+1:delay+N);
v_filtered_trunc = v_filtered_trunc_delayed(delay+1:delay+N);
v_filtered_trunc_trunc = v_filtered_trunc_delayed_trunc(delay+1:delay+N);

figure;
stairs(t, v_filtered_delayed(t+1),'b'); hold on;
stairs(t, v_filtered_trunc_delayed(t+1),'r');
stairs(t, v_filtered_trunc_delayed_trunc(t+1),'g');

legend('Ideal filter','Truncated coeff filter','trucated cofficient & output both');
title('Effect of coefficient truncation');
grid on;

%% plottting of input and filter output

figure;
plot(t,u(t+1),'b');
hold on;
plot(t,v_filtered(t+1),'r');
plot(t,v_filtered_trunc_trunc(t+1),'g');
title 'time domain';
grid on;
legend('input signal (u)', ' filter output (without truncation)','filter output (with truncation)');
% N_changed = pow2(floor(log2(length(v_filtered))));
% frequency domain
u_ffth = fft(u(1:N).*ds_hann(N)) / (sqrt(N*sum(ds_hann(N).^2)));
v_filtered_ffth = fft(v_filtered(1:N).*ds_hann(N))/ (sqrt(N*sum(ds_hann(N).^2)));
v_filtered_trunc_ffth = fft(double(v_filtered_trunc(1:N).*ds_hann(N)))/ (sqrt(N*sum(ds_hann(N).^2)));
v_filtered_trunc_trunc_ffth = fft(double(v_filtered_trunc_trunc(1:N).*ds_hann(N)))/ (sqrt(N*sum(ds_hann(N).^2)));

figure;
semilogx(f,10*log(abs(u_ffth(f+1))),'r');
hold on;
semilogx(f,10*log(mag_Xh(f+1)),'g');
semilogx(f,10*log(abs(v_filtered_ffth(f+1))),'b');
semilogx(f,10*log(abs(v_filtered_trunc_ffth(f+1))));
semilogx(f,10*log(abs(v_filtered_trunc_trunc_ffth(f+1))));

grid on;
legend ('input signal','DSM output without filtering' ,'filtered ouput after DSM','filtered output with trucated coefficient','truncated filter output');

% manual power after filtering calculation
filtered_signal_power = 2*sum(abs(v_filtered_ffth(f1:f1+2).^2));
filtered_total_power = 2*sum(abs(v_filtered_ffth(4:N/2)).^2);
filtered_noise_power = filtered_total_power - filtered_signal_power;
snr_after_filtering = 10*log10(filtered_signal_power/filtered_noise_power)

% power calculation for filter output with trucated filter coefficiet
filtered_signal_power_trunc = 2*sum(abs(v_filtered_trunc_ffth(f1:f1+2).^2));
filtered_total_power_trunc = 2*sum(abs(v_filtered_trunc_ffth(4:N/2)).^2);
filtered_noise_power_trunc = filtered_total_power_trunc - filtered_signal_power_trunc;
snr_after_filtering_trunc = 10*log10(filtered_signal_power_trunc/filtered_noise_power_trunc)

% power calculation for filter output with trucated filter coefficiet
filtered_signal_power_trunc_trunc = 2*sum(abs(v_filtered_trunc_trunc_ffth(f1:f1+2).^2));
filtered_total_power_trunc_trunc = 2*sum(abs(v_filtered_trunc_trunc_ffth(4:N/2)).^2);
filtered_noise_power_trunc_trunc = filtered_total_power_trunc_trunc - filtered_signal_power_trunc_trunc;
snr_after_filtering_trunc_trunc = 10*log10(filtered_signal_power_trunc_trunc/filtered_noise_power_trunc_trunc)

%% decimation
figure;
v_decimated = decimate(v_filtered(1:N),osr,'fir');
v_trunc_trunc_decimated = decimate(v_filtered_trunc_trunc(1:N),osr,'fir');

t1 = (0:osr:N-1);
stairs(t, u(t+1),'g'); hold on; 
stairs(t1,v_decimated(1:length(t1)),'r');
stairs(t1,v_trunc_trunc_decimated(1:length(t1)),'b');
legend ('input signal','decimated signal without truncation','decimated signal after truncation')

xlim ([0 length(t)-1]);
title 'time domain comparision of input and output after decimation'
figure;
v_decimated_hann = v_decimated.*ds_hann(length(v_decimated));
v_trunc_trunc_decimated_hann = v_trunc_trunc_decimated.*ds_hann(length(v_trunc_trunc_decimated));

v_decimated_fft = abs(fft(v_decimated_hann,length(v_decimated)))/(sqrt(length(v_decimated)*sum(ds_hann(length(v_decimated)).^2)));
v_trunc_trunc_decimated_fft = abs(fft(v_trunc_trunc_decimated_hann,length(v_trunc_trunc_decimated)))/(sqrt(length(v_trunc_trunc_decimated)*sum(ds_hann(length(v_trunc_trunc_decimated)).^2)));

f_axis_new = (0:(length(v_decimated_fft)/2)-1);
semilogx(f_axis_new,10*log(abs(v_decimated_fft(1:length(v_decimated_fft)/2)))); hold on;
semilogx(f_axis_new,10*log(abs(v_trunc_trunc_decimated_fft(1:length(v_decimated_fft)/2))));
legend ('decimated signal without truncation','decimated signal after truncation')
title 'frequency domain comparision of input and output after decimation'

% manual power after decimation calculation
decimated_signal_power = 2*sum(abs(v_decimated_fft(f1:f1+2)).^2);
decimated_total_power = 2*sum(abs(v_decimated_fft(4:end/2)).^2);
decimated_noise_power = decimated_total_power - decimated_signal_power;
snr_after_decimation = 10*log10(decimated_signal_power/decimated_noise_power)

% manual power after decimation calculation of truncated signal
trunc_trunc_decimated_signal_power = 2*sum(abs(v_trunc_trunc_decimated_fft(f1:f1+2)).^2);
trunc_trunc_decimated_total_power = 2*sum(abs(v_trunc_trunc_decimated_fft(4:end/2)).^2);
trunc_trunc_decimated_noise_power = trunc_trunc_decimated_total_power - trunc_trunc_decimated_signal_power;
trunc_trunc_snr_after_decimation = 10*log10(trunc_trunc_decimated_signal_power/trunc_trunc_decimated_noise_power)
