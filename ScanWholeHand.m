%% Scanning the whole hand
% Created on 06/26/2018 
% -------------------------------------------------------------------------
clearvars -except s32AO
close all
warning ('off','all');
% -------------------------------------------------------------------------
Fs = 40000; % Sampling frequency = 40 kS/sec
P2PAmp = 4.0; % Peak to peak ammplitude (Volt)
trial_T = 0.5; % Duration of each trial 
win_len = 0.25*Fs; % Each stimuli last 0.25 seconds
% -------------------------------------------------------------------------
% Initialize the NI DAQ interface to drive the vibration actuators and
% acquiring data from sensors
if ~exist('sNI','var')
    sNI = daq.createSession('ni');

    % Analog output channel (actuators)
    ao0 = sNI.addAnalogOutputChannel('Dev3','ao0','Voltage');
    
% %     % Analog input channel (LDV)
% %     ch0 = sNI.addAnalogInputChannel('Dev3','ai0','Voltage');
% %     ch0.TerminalConfig = 'SingleEndedNonReferenced';
% %     ch0.Range = [-5,5]; % Input voltage range from -5 to 5V
end
sNI.Rate = Fs;

%--------------------------------------------------------------------------
% Prepare stimuli: Chirp
Ts = 1/Fs;
t = 0:Ts:6;
chirp_f0 = 50;
chirp_f1 = 2000;
Chirp_sig = chirp(t,chirp_f0,t(end),chirp_f1,'logarithmic'); 
Chirp_sig = diff(Chirp_sig,2); % Second order derivative
Chirp_sig = P2PAmp.*Chirp_sig'./max(abs(Chirp_sig));
ind = find(Chirp_sig>=0,1);
Chirp_sig(1:ind) = 0;

% Prepare stimuli: Sinewaves (1/3 Octave)
fn = [40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,2000];
f_num = length(fn);
DSW_sig = discreteSineSweep(fn,Fs,P2PAmp,trial_T,win_len);

% Prepare stimuli: Sawtooth
f_saw = [100, 200, 250, 300];
[ST_sig, ~] = sawtoothImpulse(f_saw,Fs,P2PAmp,trial_T,win_len);

% Prepare stimuli: Impulse
f_imp = [10, 50, 100];
[~, Imp_sig] = sawtoothImpulse(f_imp,Fs,P2PAmp,2*trial_T,2*win_len);

% Prepare stimuli: Bandpassed random noise
f_bprn = [50, 100, 200, 300, 500, 1000];
RN_sig = randn(trial_T*Fs,1);

BPRN_Sig = [];
for i = 1:length(f_bprn)
    temp = bandpass(RN_sig,[f_bprn(i)-20,f_bprn(i)+20],Fs);
    temp((win_len+1):end) = 0;
    BPRN_Sig = [BPRN_Sig;temp];
end
BPRN_Sig = P2PAmp.*BPRN_Sig./max(abs(BPRN_Sig));

% Grouping all signals
outQueue = [Chirp_sig;DSW_sig;ST_sig;Imp_sig;BPRN_Sig];

out_len = size(outQueue,1);
fprintf('Total output time = %.2f secs\n',out_len/Fs);

%% Test
% test_sig = RN_sig;
% sig_len = length(test_sig);
% bin_num = sig_len/2+1;
% f = Fs/2*linspace(0,1,bin_num); % Frequency scale
% Y_fft = fft(test_sig,sig_len)/sig_len; % FFT     
% abs_Y_fft = 2*abs(Y_fft(1:bin_num)); % Spectrum     
% plot(f,abs_Y_fft)
% % xlim([0 1001])
% xlabel('Frequency (Hz)')
% ylabel('Amplitude (Volt)')
%--------------------------------------------------------------------------
%% Run measurement
if ~isempty(find(abs(outQueue)>4, 1))
    error('Output must be within 4 V!');
else    
    queueOutputData(sNI,outQueue);

    sNI.startForeground; 
end
%--------------------------------------------------------------------------
disp('        [ End of Measurement ]        ')

