clear
close all

%%Parameters
samplesPerBit = 256;
num_bits = 256;
Tb = 1/6600; %symbol period
Ts = 2*Tb;
fb = 1/Tb; %symbol rate
kT = 5;
N = samplesPerBit*num_bits; %number of samples
fs = samplesPerBit/Tb;
r = .6;
fc = 10000;
num_trials = 10;

min_EbN0 = -40;
max_EbN0 = 10;

num_carriers = 2; %number of QPSK carriers

%%
%Code starts here

freqRange = linspace(-fs/2, fs/2, samplesPerBit*num_bits);
T = (1/fs)*samplesPerBit*num_bits;

[pulse, pulse_time] = root_rcro(kT, Ts, samplesPerBit, r);

figure();
plot(pulse_time, pulse)
xlabel('Time (sec)')
ylabel('Amplitude (V)')
title(strcat('Root RCRO pulse with r =', num2str(r)))

Eb = (pulse_time(2) - pulse_time(1))*sum(pulse);


neg_val = -1;
pos_val = 1;

all_error_trials = zeros(num_trials, max_EbN0 - min_EbN0 + 1);

avg_psd = zeros(1, length(num_bits*samplesPerBit));

for j = 1:num_trials
    j
    all_error_rates = zeros(1, max_EbN0 - min_EbN0 + 1);

    for z = min_EbN0:max_EbN0
        EbN0 = z;
        noise_amp = kT*10^(-EbN0*(1/20)); %To be automatically calculated
        bits = bit_gen(num_bits, neg_val, pos_val);

        reshaped_bits = reshape(bits, (num_bits/(2*num_carriers)), 2*num_carriers);

        for k = 1:num_carriers
            reshaped_bits(:,2*k) = 1i*reshaped_bits(:,2*k);
        end

        %Add the real and imaginary numbers. There are four columns because of two
        %carriers, each with two parallel values to compute the IFFT
        qpsk_bits = [reshaped_bits(:,1)+reshaped_bits(:,2) reshaped_bits(:,3)+reshaped_bits(:,4)];

        tx_ifft = ifft(transpose(qpsk_bits));

        tx_serial = reshape(tx_ifft, 1, 2*length(tx_ifft));

        tx_i = real(tx_serial);
        tx_q = imag(tx_serial);

        i_rcro = generate_baseband(tx_i, pulse, 2*kT, Tb, samplesPerBit);
        q_rcro = generate_baseband(tx_q, pulse, 2*kT, Tb, samplesPerBit);

        carriers = [cos(2*pi*fc*(1:length(i_rcro))*(1/fs)); sin(2*pi*fc*(1:length(i_rcro))*(1/fs))];

        ofdm_signal = i_rcro.*carriers(1,:) - q_rcro.*carriers(2,:);

        msg_spectrum = fft(ofdm_signal, samplesPerBit*num_bits);
        msg_spectrum = fftshift(msg_spectrum);
        msg_spectrum = (1/fs)^2 .* (abs(msg_spectrum)).^2;
        avg_psd = avg_psd +  msg_spectrum;

        ofdm_noisy = ofdm_signal + noise_amp*randn(1, length(ofdm_signal));

        %%
        %RECEIVER CODE STARTS HERE

        rx_demod_i = ofdm_noisy .* carriers(1,:);
        rx_demod_q = ofdm_noisy .* carriers(2,:);

        rx_matched_i = conv(rx_demod_i, pulse);
        rx_matched_q = conv(rx_demod_q, pulse);

        sampled_i = zeros(1, num_bits/2);
        sampled_q = zeros(1, num_bits/2);

        index = length(pulse);
        for i = 1:num_bits/2
            sampled_i(i) = rx_matched_i(index);
            sampled_q(i) = rx_matched_q(index);
            index = index + samplesPerBit;
        end

        rx_serial = sampled_i - 1i*sampled_q;
        rx_parallel = reshape(rx_serial, 2, length(rx_serial)/2);

        rx_fft = fft(rx_parallel);
        rx_fft = transpose(rx_fft);

        rx_reshaped_bits = [real(rx_fft(1:end,1)) imag(rx_fft(1:end,1)) real(rx_fft(1:end, 2)) imag(rx_fft(1:end, 2))];
        rx_bits = reshape(rx_reshaped_bits, 1, num_bits);

        for h = 1:length(rx_bits)
            if rx_bits(h) >=0
                rx_bits(h) = 1;
            else
                rx_bits(h) = -1;
            end
        end

        error_rate = sum(abs(rx_bits - bits)/2)/num_bits;
        all_error_rates(z-min_EbN0+1) = error_rate;
    end %End of a single trial

    all_error_trials(j, 1:end) = all_error_rates;

end

avg_error_rate = mean(all_error_trials, 1);

figure();
plot((1/fs)*(1:length(i_rcro)),i_rcro)
xlabel('Time (Sec)')
ylabel('Amplitude (V)')
title('I Channel RCRO Message')

figure();
plot((1/fs)*(1:length(carriers(1, 1:end))), carriers(1, 1:end))
xlabel('Time (Sec)')
ylabel('Amplitude (V)')
title('I Channel Cosine Carrier')

figure();
plot((1/fs)*(1:length(ofdm_signal)), ofdm_signal)
xlabel('Time (Sec)')
ylabel('Amplitude (V)')
title('Noise Free OFDM Signal')

figure();
plot((1/fs)*(1:length(ofdm_noisy)), ofdm_noisy)
xlabel('Time (Sec)')
ylabel('Amplitude (V)')
title('Noisy OFDM Signal')

figure();
plot((1/fs)*(1:length(rx_demod_i)), rx_demod_i)
xlabel('Time (Sec)')
ylabel('Amplitude (V)')
title('Demodulated I Channel')

figure();
plot((1/fs)*(1:length(rx_matched_i)), rx_matched_i)
xlabel('Time (Sec)')
ylabel('Amplitude (V)')
title('Match Filtered I Channel')

figure();
plot(-40:1:10, avg_error_rate)
hold on
EbN0 = 10.^((-40:1:10)/10);
plot(-40:1:10, qfunc(sqrt(2*(EbN0))))
xlabel('Eb/No (dB)')
ylabel('Probability of Error')
title('BER for OFDM w/ QPSK')
legend('Measured', 'Ideal')


% Generate theoretical PSD
% Bandpass signal PSD has shape of RCRO
% Two carriers at frequencies fc1 & fc2
Ts = 2* Tb;
T_ofdm = 2*Ts;
fc1 = fc + 1/T_ofdm;
fc2 = fc - 1/T_ofdm;
theoretical_1 = zeros(1,length(freqRange));
theoretical_2 = zeros(1,length(freqRange));
f0 = 1/(2*Ts);
fd = r*f0;
f1 = f0 - fd;
B = fd +f0;
theoretical_1(abs(freqRange - fc1) < f1) = 1;
theoretical_1((abs(freqRange - fc1) > f1) & ...
  (abs(freqRange - fc1) < B)) = 1 * 1/2 * (1+cos( ...
  pi/(2*fd)*(abs(freqRange((abs(freqRange - fc1) > f1) & ...
  (abs(freqRange - fc1) < B))- fc1) - f1)));

theoretical_2(abs(freqRange - fc2) < f1) = 1;
theoretical_2((abs(freqRange - fc2) > f1) & ...
  (abs(freqRange - fc2) < B)) = 1 * 1/2 * (1+cos( ...
  pi/(2*fd)*(abs(freqRange((abs(freqRange - fc2) > f1) & ...
  (abs(freqRange - fc2) < B))- fc2) - f1)));

% Add 1/4 factor because of upconversion.
% PSD is proportional to FFT^2 so (1/2)^2 = 1/4
% Ts is twice Tb because 2 carriers means each channel
% has a symbol time of 2*Tb.
% T_ofdm is even longer because we split the data stream between
% two different carriers.
theoretical_1 = 1/4 .* theoretical_1/ T_ofdm;
theoretical_2 = 1/4 .* theoretical_2/ T_ofdm;
theoretical = theoretical_1 + theoretical_2;
% Display PSD Graph
figure();
% avg_psd = avg_psd ./ (num_trials* (max_EbN0 - min_EbN0 + 1) *T);
avg_psd = avg_psd ./ (num_trials *T);

scaleFactor = max(abs(avg_psd))/max(theoretical);
theoretical = scaleFactor*theoretical;
plot(freqRange, abs(avg_psd), freqRange, theoretical);
% xlim([fc-fb, fc+fb])
xlim([5000, 15000])
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Frequency Spectrum for OFDM w/ QPSK')
legend('Measured', 'Ideal')

figure();
plot(freqRange, 10*log(abs(avg_psd)))
xlim([0, 20000])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Frequency Spectrum for OFDM w/ QPSK in dBW')