%% Prepare Customized Signals for Polytec Signal Generator
% Created on 09/09/2018
% -------------------------------------------------------------------------
Fs = 10000; % Sampling Frequency (Hz)
dt = 1/Fs; % Sampling duration (secs)
SignalDuration = 1; % Signal duration

pauseSig = zeros(0.1*Fs,1); % Pause 100 ms between trials

t = dt:dt:SignalDuration; % Time stamps

sine_freq = 4000;

sinewave = sin(2*pi*sine_freq.*t);

[ spect, freq ] = powerSpectrum(sinewave, Fs);

figure
subplot(2,1,1);
plot(t,sinewave);
subplot(2,1,2);
plot(freq,spect);

% Output signal as .csv file
file_id = fopen(sprintf('sine%dhz.txt',sine_freq), 'w');
fprintf(file_id, 'ArbWave = \n');
fprintf(file_id, '%f\n', sinewave);
fclose(file_id);


%% Prepare signal: Bandpassed random noise
f_bprn = [20, 40, 80, 160, 320, 640];
RN_sig = randn(SignalDuration*Fs,1);

BPRN_Sig = [];
for i = 1:length(f_bprn)
    band_width = 0.1*f_bprn(i);
    temp = bandpass(RN_sig,[f_bprn(i)-band_width,f_bprn(i)+band_width],Fs);
    BPRN_Sig = [BPRN_Sig;temp;pauseSig];
end

%% Local function: FFT - Spectrum
function [ abs_Y_fft, f ] = powerSpectrum( input_data, Fs )
L=length(input_data); % Sampling data number
NFFT = 2^nextpow2(L); % Next power of 2 from length of y

Y_fft = fft(input_data, NFFT)/L; % FFT 
f = Fs/2*linspace(0, 1, NFFT/2+1); % Frequency scale

abs_Y_fft = 2*abs(Y_fft(1:NFFT/2+1)); % Power spectrum 
end