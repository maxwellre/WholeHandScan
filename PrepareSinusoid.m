%% Prepare Sinusoid
% Created on 10/01/2018
% -------------------------------------------------------------------------
% Configuration
Fs = 125000; % Sampling Frequency (Hz)

SineFs = 1000; % Sinusoid Frequency (Hz)

% -------------------------------------------------------------------------
SignalDuration = 1/SineFs; % Duration of one Sinusoid cycle (secs)
dt = 1/Fs; % Sampling duration (secs)

t = (0:dt:SignalDuration)'; % Time stamps for one cycle of sinusoid 
sineSig = sin(2*pi*SineFs.*t); % One cycle of sinusoid

beginSig = zeros(0.1*Fs,1); % Pause 100 ms at the beginning
pauseSig = zeros(0.2*Fs,1); % Pause 200 ms between trials

% Output signal: group of four sinusoids with 200ms interval between them
outSig = [beginSig;...
    sineSig;pauseSig;...
    sineSig;pauseSig;...
    sineSig;pauseSig;...
    sineSig;pauseSig]; % Signal to be used

t = (1:length(outSig))/Fs; % Corresponding time stamps

% -------------------------------------------------------------------------
% Write the output signal into a .txt file
SN = length(outSig); % Sample Number
RF = Fs/SN; % Repeating Frequency (Hz)

% Name of the saved .txt file
saveName = sprintf('Sine%dHz_SN-%d_RF-%.6f.txt',SineFs,SN,RF);
disp(['Saved to: ',saveName]);

% Write the signal
file_id = fopen(saveName, 'w');
fprintf(file_id, 'ArbWave = \n');
fprintf(file_id, '%f\n', outSig);
fclose(file_id);