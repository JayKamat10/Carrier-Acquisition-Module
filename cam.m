
%******Carrier Acquisition Module QPSK*****

clc
clear all
fprintf('*****Carrier Acquisition Module***** \n');
%Generating Input
type = input('Enter 1 for BPSK, or 2 for QPSK: \n');

bit_rate = input('Enter Bit Rate:');
fs = 1024*bit_rate;
delta_fc = input('Enter tolerable acquisition error: \n');
D = 500;
fc_max = 5*fs/(2*D);

fprintf('Centre Frequency of carrier must be lesser than %d corresponding to a Decimation Factor of 500\n',fc_max);

fc = input('Enter Centre Frequency of Carrier: \n');     %Carrier Frequency

bit_dt = 1/bit_rate;
dt = 1/fs;
duration = input('Enter duration of signal: \n');    %Time duration of signal

t = 0:dt:duration-dt;


N = duration * bit_rate; %Number of bits to be generated
n = randi([0 1], 1, N); 

%*****Mapping of bits in NRZ form*****
for index = 1:N
    if n(index) == 0
        nn(index) = -1
    else 
        nn(index) = 1
    end
end

%*****Conversion to Continuous time BPSK NRZ signal*****
fsb = fs/bit_rate ; %Samples per bit
i = 1; %Index of input bits
st = 0 : 1/fsb : N-1/fsb; %Time
for j = 1 : length(st) %Index of time array
    if st(j) <= i
        signal(j) = nn(i);
    else
        signal(j) = nn(i);
        i = i + 1;
    end
end


%*****Making the Root Raised Cosine filter (Default value is sqrt for the given function)***** 
rrc_filt = comm.RaisedCosineTransmitFilter("FilterSpanInSymbols",4,"RolloffFactor",0.35,"OutputSamplesPerSymbol",1024);
rrc_sig = rrc_filt(nn');
rrc_sig = rrc_sig';


carrier0 = cos(2*pi*fc*t);  %Carrier Wave of 900 Hz
carrier1 = cos(2*pi*fc*t +pi/4)
carrier2 = cos(2*pi*fc*t + 3*pi/4);
carrier3 = cos(2*pi*fc*t + 5*pi/4);
carrier4 = cos(2*pi*fc*t + 7*pi/4);

%*****BPSK Modulation of rrc_signal*****
x_bpsk = rrc_sig .* carrier0;

%*****QPSK Modulation of rrc_signal*****
for a = 1:fsb:length(st)-fsb
    if(signal(a)== -1 || signal(a+fsb) == -1) 
        x_qpsk(a:a+2*fsb-1) = rrc_sig(a:a+2*fsb-1) .* carrier1(a:a+2*fsb-1);
        
    elseif(signal(a)== -1 || signal(a+fsb) == 1) 
        x_qpsk(a:a+2*fsb-1) = rrc_sig(a:a+2*fsb-1) .* carrier2(a:a+2*fsb-1);
        
    elseif(signal(a)== 1 || signal(a+fsb) == -1) 
            x_qpsk(a:a+2*fsb-1) = rrc_sig(a:a+2*fsb-1) .* carrier3(a:a+2*fsb-1);
   
    else 
        x_qpsk(a:a+2*fsb-1) = rrc_sig(a:a+2*fsb-1) .* carrier4(a:a+2*fsb-1);
       
    end
end

if type == 1
    x = x_bpsk;
    fprintf('*****Showing results of BPSK Modulated Signal***** \n');
else
    x = x_qpsk;
    fprintf('*****Showing results of QPSK Modulated Signal***** \n');
end

x_sq = x .* x;         %Squared Signal BPSK

x_sq_mixed = x_sq .* cos(2*pi*0.9*2*fc*t);

x_sq_filt = bandpass(x_sq_mixed, [2*0.1*(fc-delta_fc) 2*0.1*(fc+delta_fc)], fs);  %Bandpass filtered output


l = length(t);  %Length of time axis

nfft = 1024;
f = (fs/nfft)*(0:nfft/2-1); %Frequency axis



decimated_op = decimate(x_sq_filt,D);

%*****Fourier transforms*****


%Input signal
INPUT = abs(fft(signal, nfft));
INPUT = INPUT(1:nfft/2);
INPUT_RRC = abs(fft(rrc_sig, nfft));
INPUT_RRC = INPUT_RRC(1:nfft/2);


%Squared Waveform
SQUARED = abs(fft(x_sq, nfft));
SQUARED = SQUARED(1:nfft/2);

%Filtered Waveform
FILTERED = abs(fft(decimated_op, nfft));
FILTERED = FILTERED(1:nfft/2);

%*****PLOTS*****

figure;
title('Original Signal in Time Domain');
subplot 211; plot(t,signal); xlabel('time'); ylabel('Amplitude'); title('NRZ Signal');
subplot 212; plot(t,rrc_sig); xlabel('time'); ylabel('Amplitude'); title('Pulse Shaped Signal');

figure;
title('Original Signal in Frequency Domain')
subplot 211; plot(f,INPUT); xlabel('frequency'); ylabel('Amplitude'); title('NRZ Signal');
subplot 212; plot(f,INPUT_RRC); xlabel('frequency'); ylabel('Amplitude'); title('Pulse Shaped Signal');

figure;
title('Squared Signal')
subplot 211; plot(t,x_sq); xlabel('time'); ylabel('Amplitude'); title('Squared Signal');
subplot 212; plot(f,SQUARED); xlabel('frequency'); ylabel('Amplitude'); title('Squared Signal');

figure;
title('Filtered Signal');
subplot 211; plot(t,x_sq_filt); xlabel('time'); ylabel('Amplitude'); title('Filtered Signal');
subplot 212; plot(f,FILTERED); xlabel('frequency'); ylabel('Amplitude'); title('Filtered Signal');


%*****RESULTS*****
[Max, k] = max(FILTERED); %finding bin of interest

fprintf('Actual carrier frequency is %f \n',fc);
fprintf('Estimate (without any interpolation) = %f \n',(10*(k-1)*fs)/(2*D*nfft));

fprintf('\n');
fprintf('\n');

fprintf('Note: The way I am performing decimation is as follows: \n');
fprintf('Firstly I am mixing 2fc (squared signal) with a 1.8fc signal, so that I obtain 0.2fc after filtering \n');
fprintf(' I decimate this mixed signal by a factor of 500, and obtain fc');
