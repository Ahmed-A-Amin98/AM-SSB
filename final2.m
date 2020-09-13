%%%%%%%%%%%%%%%%%%%%%%%%%%%    1      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read the attached audio file , plot the spectrum of this signal. 

[audio, Fs] = audioread('eric.wav');
s = length(audio) / Fs;
t = linspace(0, s, s*Fs + 1);
fs = linspace(-Fs/2, Fs/2, s*Fs + 1);

figure;
subplot(2, 1, 1)
plot(t, audio);
ylim([-0.5, 0.5]);
title('original signal Before Filter in time domaim');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(audio))));
xlim([-1e4, 1e4]);
ylim([-500, 500]);
title('original signal Before Filter in freq. domaim');
%%%%%%%%%%%%%%%%%%%%%%%%%%%    2      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use an ideal filter to remove all frequencies greater than 4KH
d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', 4e3, 'SampleRate', Fs);
filteredSig = filter(d, audio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%    3      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(2, 1, 1);
plot(t, filteredSig);
ylim([-0.5, 0.5]);
title('original signal After Filter in time domaim');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(filteredSig))));
xlim([-0.5e4, 0.5e4]);
ylim([-500, 500]);
title('original signal after Filter in freq. domaim');

%sound( real(double(filteredSig)), Fs);



%%%%%%%%%%%%%%%%%%%%%%%%%%%    4      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DSB_SC

Fc = 100e3;
message = resample(filteredSig, 5 * Fc, Fs); % Upsampling signal to 5 * Fc
Fs = 5*Fc;
s = length(message)/Fs;
t = linspace(0, s, s*Fs);
fs = linspace(-Fs/2, Fs/2, s*Fs);

%DSBSC
carrier = cos(2*pi*Fc*t);
modulatedSig = message .* transpose(carrier);

figure; 
subplot(2,1,1);
plot(t, modulatedSig);
title('DSB-SC Time Domain')
subplot(2,1,2); 
plot(fs, real(fftshift(fft(modulatedSig))));
xlim([-1.2e5, 1.2e5]);
ylim([-2.5e3, 2.5e3])
title('DSB-SC Frequency Domain')


%%%%%%%%%%%%%%%%%%%%%%%%%%%    5      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SSB_SC using ideal filter
%SSBSC
d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', Fc, 'SampleRate', Fs);
modulatedSCSig  = filter(d, modulatedSig);


figure; 
subplot(2,1,1);
plot(t, modulatedSCSig);
title('SSB-SC Time Domain')
subplot(2,1,2); 
plot(fs, real(fftshift(fft(modulatedSCSig))));
xlim([-1.2e5, 1.2e5]);
ylim([-2.5e3, 2.5e3])
title('SSB-SC Frequency Domain')


%%%%%%%%%%%%%%%%%%%%%%%%%%%    6      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

carrier = cos(2*pi*Fc*t);
demodSignal = modulatedSCSig .* transpose(carrier);
d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', 4e3, 'SampleRate', 5 * Fc);
coherentSC  = filter(d, demodSignal);

figure;
subplot(2, 1, 1);
plot(t, coherentSC);
ylim([-0.05, 0.05]);
title('demodulated signal using coherent in time domain');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(coherentSC))));
xlim([-0.5e4, 0.5e4]);
ylim([-1e3, 1e3]);
title('demodulated signal using coherent in freq. domain');

%sound( real(double(resample(coherentSC, Fs, 5 * Fc))), Fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    7   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SSB_SC using butter filter
%SSBSC
Fs = 5*Fc;
Fn=Fs/2;
Fl=(Fc-4e3)/Fn;
Fh=Fc/Fn;
Fcut=[Fl,Fh]; 

[b, a] = butter(4,Fcut);
SSBSC_practical = filter(b, a, modulatedSig);

figure; 
subplot(2,1,1);
plot(t, SSBSC_practical);
title('SSB-SC practical Time Domain')
subplot(2,1,2); 
plot(fs, real(fftshift(fft(SSBSC_practical))));
xlim([-1.2e5, 1.2e5]);
ylim([-2.5e3, 2.5e3])
title('SSB-SC practical Frequency Domain')


carrier = cos(2*pi*Fc*t);
demodSignal_practical = SSBSC_practical .* transpose(carrier);
d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', 4e3, 'SampleRate', 5 * Fc);
coherentSC_practical  = filter(d, demodSignal_practical);

figure;
subplot(2, 1, 1);
plot(t, coherentSC_practical);
ylim([-0.05, 0.05]);
title('demodulated signal using coherent in time domain');
subplot(2, 1, 2);
plot(fs, real(fftshift(fft(coherentSC_practical))));
xlim([-0.5e4, 0.5e4]);
ylim([-1e3, 1e3]);
title('demodulated signal using coherent in freq. domain');

%sound( real(double(resample(coherentSC_practical, Fs, 5 * Fc))), Fs);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    8   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generating noise signals
Signal_with_noise_1 = awgn(modulatedSCSig,0,'measured'); %0 db SNR
Signal_with_noise_2 = awgn(modulatedSCSig,10,'measured'); %10 db SNR
Signal_with_noise_3 = awgn(modulatedSCSig,30,'measured'); %30 db SNR

%demodulating
noisy_signal_demodulated_1 = Signal_with_noise_1.*transpose(carrier);
noisy_signal_demodulated_1_resampled = resample(noisy_signal_demodulated_1, Fs, 5 * Fc);
%sound( real(double(noisy_signal_demodulated_1_resampled)), Fs);

noisy_signal_demodulated_2 = Signal_with_noise_2.*transpose(carrier);
noisy_signal_demodulated_2_resampled = resample(noisy_signal_demodulated_2, Fs, 5 * Fc);
%sound( real(double(noisy_signal_demodulated_2_resampled)), Fs);


noisy_signal_demodulated_3 = Signal_with_noise_3.*transpose(carrier);
noisy_signal_demodulated_3_resampled = resample(noisy_signal_demodulated_3, Fs, 5 * Fc);
%sound( real(double(noisy_signal_demodulated_3_resampled)), Fs);



%plot them in time and frequency domain

figure; 
subplot(2,1,1);
plot(t, noisy_signal_demodulated_1);
title('recived SSB-SC with noise 0 in time domain');
subplot(2,1,2); 
plot(fs, real(fftshift(fft(noisy_signal_demodulated_1))));
xlim([-1.2e5, 1.2e5]);
ylim([-2.5e3, 2.5e3])
title('recived SSB-SC with noise 0 in freq. domain');

figure; 
subplot(2,1,1);
plot(t, noisy_signal_demodulated_2);
title('recived SSB-SC with noise 10 in time domain');
subplot(2,1,2); 
plot(fs, real(fftshift(fft(noisy_signal_demodulated_2))));
xlim([-1.2e5, 1.2e5]);
ylim([-2.5e3, 2.5e3])
title('recived SSB-SC with noise 10 in freq. domain');

figure; 
subplot(2,1,1);
plot(t, noisy_signal_demodulated_3);
title('recived SSB-SC with noise 30 in time domain');
subplot(2,1,2); 
plot(fs, real(fftshift(fft(noisy_signal_demodulated_3))));
xlim([-1.2e5, 1.2e5]);
ylim([-2.5e3, 2.5e3])
title('recived SSB-SC with noise 30 in freq. domain');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  9  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ferror = 100100;
carrier_Signal_Error = cos(2*pi*Ferror*t + degtorad(10));

noisy_signal_demodulated_3_with_phase_freq_err = Signal_with_noise_3.*transpose(carrier_Signal_Error);
noisy_signal_demodulated_3_with_phase_freq_err_resampled = resample(noisy_signal_demodulated_3_with_phase_freq_err, Fs, 5 * Fc);


figure;
plot(noisy_signal_demodulated_3_with_phase_freq_err);
title('recived SSB-SC with noise 30 and phase err and freq err in time domain');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  10  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generating of DSBTC
m=.5;
Mc = max(message);
normalized_filtered_signal_t=(message)/Mc;
Ac=2*Mc;
carrierSignal2 = Ac*cos(2*pi*Fc*t);
signal_resampled2 = resample(normalized_filtered_signal_t,5 * Fc, Fs);
signal_to_modulate2 = signal_resampled2(1:length(carrierSignal2));
DSBTC =  (1+m*signal_to_modulate2).* transpose(carrierSignal2);

%generating of SSBTC using ideal filter

d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', Fc, 'SampleRate', Fs);
modulatedSCSig2  = filter(d, DSBTC);

received_envelope_SSBTC = abs(hilbert(modulatedSCSig2));

figure;
plot(received_envelope_SSBTC);
title('SSB-TC demodulated signal');

%sound( real(double(resample(received_envelope_SSBTC, Fs, 5 * Fc))), Fs);




