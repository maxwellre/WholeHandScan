%% Prepare Customized Signals for Polytec Signal Generator
% Created on 09/09/2018
% -------------------------------------------------------------------------
Fs = 10000; % Sampling Frequency (Hz)
dt = 1/Fs; % Sampling duration (secs)
SignalDuration = 1; % Signal duration (secs)

pauseSig = zeros(0.1*Fs,1); % Pause 100 ms between trials

t = (dt:dt:SignalDuration)'; % Time stamps

signalType = 4;

switch signalType
    case 0
        sine_freq = 4000;
        outSig = sin(2*pi*sine_freq.*t);
        
        saveName = sprintf('Sine%dhz',sine_freq);
        
    case 1
        f_bprn = [20, 40, 80, 160, 320, 640];
        cycle_num = 20; 

        outSig = [];
        for i = 1:length(f_bprn)
            RN_sig = randn(round(cycle_num*Fs/f_bprn(i)),1);
            band_width = 0.1*f_bprn(i);
            
            bpFilter = designfilt('bandpassfir', ...       % Response type
               'FilterOrder',100, ...            % Filter order
               'StopbandFrequency1',f_bprn(i)-band_width-5, ... % Frequency constraints
               'PassbandFrequency1',f_bprn(i)-band_width, ...
               'PassbandFrequency2',f_bprn(i)+band_width, ...
               'StopbandFrequency2',f_bprn(i)+band_width+5, ...
               'DesignMethod','ls', ...         % Design method
               'SampleRate',Fs);               % Sample rate
           
            temp = filtfilt(bpFilter,RN_sig);
            
%             [spect,freq]=powerSpectrum(temp, Fs);plot(freq,spect);xlim([0 700])
            
            outSig = [outSig;temp;pauseSig];
        end
        
        saveName = 'BandpassRandomNoise';
        
    case 2
        impulse_width = 0.0005; % (secs)
        impulse_period = 0.1; % (secs)
        impulse_num = SignalDuration/impulse_period;

        temp = zeros(impulse_period*Fs,1);
        temp(1:impulse_width*Fs) = 1;

        outSig = repmat(temp,[impulse_num,1]);
        
        saveName = 'ImpulseTrain';
        
    case 3
        impulse_width = 0.01; % (secs)
        impulse_period = 0.1; % (secs)
        impulse_num = SignalDuration/impulse_period;

        temp = zeros(impulse_period*Fs,1);
        temp(1:impulse_width*Fs) = 1;

        outSig = repmat(temp,[impulse_num,1]);
        
        saveName = 'TapPulseTrain';
        
    case 4
        sigma2 = 1e-5;
        pulse_length = 0.01; % (secs)
        pulse_period = 0.1; % (secs)
        cosine_freq = 250; % (Hz)
        pulse_num = SignalDuration/pulse_period;
        
        t_temp = (dt:dt:pulse_length)';
        temp = cos(2*pi*cosine_freq.*(t_temp - (0.5*pulse_length)));
        temp = temp.*exp(-(t_temp - (0.5*pulse_length)).^2./sigma2);
        temp = [temp;zeros(round((pulse_period-pulse_length)*Fs),1)];
        
        outSig = repmat(temp,[pulse_num,1]);
        
        saveName = 'GaborPulseTrain';
        
    case 5
        t_delay = (0:0.2:10)/1000; % (secs)
        impulse_width = 0.0005; % (secs)
        impulse_period = 0.1; % (secs)
        
        impulse_ind = 1:round(impulse_width*Fs);
        
        outSig = [];
        for i = 1:length(t_delay)
            temp1 = zeros(impulse_period*Fs,1);
            temp1(impulse_ind) = 1;
            temp2 = zeros(impulse_period*Fs,1);
            temp2(round(t_delay(i)*Fs)+impulse_ind) = 1;
            
            outSig = [outSig;(temp1+temp2)];
        end
        
        saveName = 'ImpulsePair';
end
RF = Fs/length(outSig); % Repeating Frequency (Hz)
saveName = [saveName,sprintf('_RF%.6f.txt',RF)];
disp(['Prepare: ',saveName]);

%% Output signal as .txt file
file_id = fopen(saveName, 'w');
fprintf(file_id, 'ArbWave = \n');
fprintf(file_id, '%f\n', outSig);
fclose(file_id);

%% Local function: FFT - Spectrum
function [ abs_Y_fft, f ] = powerSpectrum( input_data, Fs )
L=length(input_data); % Sampling data number
NFFT = 2^nextpow2(L); % Next power of 2 from length of y

Y_fft = fft(input_data, NFFT)/L; % FFT 
f = Fs/2*linspace(0, 1, NFFT/2+1); % Frequency scale

abs_Y_fft = 2*abs(Y_fft(1:NFFT/2+1)); % Power spectrum 
end