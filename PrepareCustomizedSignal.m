%% Prepare Customized Signals for Polytec Signal Generator
% Created on 09/09/2018
% -------------------------------------------------------------------------
clear all
% -------------------------------------------------------------------------
Fs = 20000; % Sampling Frequency (Hz)
dt = 1/Fs; % Sampling duration (secs)
SignalDuration = 1; % Signal duration (secs)
ExtraMeasureLen = 0.2*Fs; % 200 ms extra measurement time

beginSig = zeros(0.02*Fs,1); % Pause 20 ms at beginning trials
pauseSig = zeros(0.1*Fs,1); % Pause 100 ms between trials

t = (dt:dt:SignalDuration)'; % Time stamps

% 0: Sinusoid
% 1: Bandpass Filtered Random Noise
% 2: Impulse Train
% 3: Tapping Impulse Train
% 4: ImpulsePair
signalType = 4; 

outSig = [];
saveName = [];
switch signalType
    case 0
        sine_freq = [20,40,80,120,160,200,240,280,320,360,400,440,480,520,...
            560,600,640];
        cycle_num = 10; 

        outSig = beginSig;
        for i = 1:length(sine_freq)
            t = (dt:dt:cycle_num/sine_freq(i))'; 
            temp = sin(2*pi*sine_freq(i).*t);
            
            hann_win = hann(length(temp));
            temp = temp.*hann_win;
            
            spectAnalysis(temp, Fs);
            
            outSig = [outSig;temp;pauseSig];
        end

        saveName = sprintf('Sinusoid');
        
    case 1
        f_bprn = [20,40,80,120,160,200,240,280,320,360,400,440,480,520,...
            560,600,640];
        cycle_num = 10; 

        outSig = beginSig;
        for i = 1:length(f_bprn)
            RN_sig = randn(round(cycle_num*Fs/f_bprn(i)),1);
            band_width = (0.1*f_bprn(i));
            
%             bpFilter = designfilt('bandpassfir', ... % Response type
%                'FilterOrder',100, ...            % Filter order
%                'CutoffFrequency1',f_bprn(i)-band_width,...
%                'CutoffFrequency2',f_bprn(i)+band_width,...
%                'SampleRate',Fs);               % Sample rate
%             temp = filtfilt(bpFilter,RN_sig);
            
            temp = bandpass(RN_sig,[f_bprn(i)-band_width,...
                f_bprn(i)+band_width],Fs,...
                'ImpulseResponse','iir','Steepness',0.95);
            
            hann_win = hann(length(temp));
            temp = temp.*hann_win;
            
            spectAnalysis(temp, Fs);
            
            outSig = [outSig;temp;pauseSig];
        end
        
        saveName = 'BandpassRandomNoise';
        
    case 2
        impulse_width = 0.0005; % (secs)
        impulse_period = 0.1; % (secs)
        impulse_num = SignalDuration/impulse_period;

        temp = zeros(impulse_period*Fs,1);
        temp(1:impulse_width*Fs) = 1;

        outSig = [beginSig;repmat(temp,[impulse_num,1])];
        
        saveName = 'ImpulseTrain';
        
    case 3
        impulse_width = 0.01; % (secs)
        impulse_period = 0.1; % (secs)
        impulse_num = SignalDuration/impulse_period;

        temp = zeros(impulse_period*Fs,1);
        temp(1:impulse_width*Fs) = 1;

        outSig = [beginSig;repmat(temp,[impulse_num,1])];
        
        saveName = 'TapPulseTrain';
        
    case 4
        t_delay = (0:0.2:10)/1000; % (secs)
        impulse_width = 0.0005; % (secs)
        impulse_period = 0.1; % (secs)
        
        impulse_ind = 1:round(impulse_width*Fs);
        
        outSig = beginSig;
        for i = 1:length(t_delay)
            temp1 = zeros(impulse_period*Fs,1);
            temp1(impulse_ind) = 1;
            temp2 = zeros(impulse_period*Fs,1);
            temp2(round(t_delay(i)*Fs)+impulse_ind) = 1;
            
            outSig = [outSig;(temp1+temp2)];
        end
        
        saveName = 'ImpulsePair';
        
    case 5
        error('Gabor has been removed!')
%         sigma2 = 1e-5;
%         pulse_length = 0.01; % (secs)
%         pulse_period = 0.1; % (secs)
%         cosine_freq = 250; % (Hz)
%         pulse_num = SignalDuration/pulse_period;
%         
%         t_temp = (dt:dt:pulse_length)';
%         temp = cos(2*pi*cosine_freq.*(t_temp - (0.5*pulse_length)));
%         temp = temp.*exp(-(t_temp - (0.5*pulse_length)).^2./sigma2);
%         temp = [temp;zeros(round((pulse_period-pulse_length)*Fs),1)];
%         
%         outSig = [beginSig;repmat(temp,[pulse_num,1])];
%         
%         saveName = 'GaborPulseTrain';
        
end
SN = length(outSig) +ExtraMeasureLen; % Sample Number
RF = Fs/length(outSig); % Repeating Frequency (Hz)
saveName = [saveName,sprintf('_SN-%d_RF-%.6f.txt',SN,RF)];
disp(['Prepare: ',saveName]);

%% Output signal as .txt file
file_id = fopen(saveName, 'w');
fprintf(file_id, 'ArbWave = \n');
fprintf(file_id, '%f\n', outSig);
fclose(file_id);

%% Local function: Spectrum Analysis
function spectAnalysis(data_in, Fs)
[spect,freq]=powerSpectrum(data_in, Fs);
avg_freq = sum(freq'.*(spect/sum(spect)));

[peakSpect,peakInd] = max(spect);

plot(freq,spect);
hold on
plot(freq(peakInd),peakSpect,'.r');
text(freq(peakInd)+10,peakSpect,sprintf('%.0f Hz',freq(peakInd)));
text(freq(peakInd)+10,0.01,sprintf('Avg. Freq = %.0f Hz',avg_freq));
hold off
xlim([0 800]);
xlabel('Frequency (Hz)');
pause(1);
end
%% Local function: FFT - Spectrum
function [ abs_Y_fft, f ] = powerSpectrum( input_data, Fs )
L=length(input_data); % Sampling data number
NFFT = 2^nextpow2(L); % Next power of 2 from length of y

Y_fft = fft(input_data, NFFT)/L; % FFT
f = Fs/2*linspace(0, 1, NFFT/2+1); % Frequency scale

abs_Y_fft = 2*abs(Y_fft(1:NFFT/2+1)); % Power spectrum
end