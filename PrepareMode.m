%% Touch Amp Signal Preparation
Path_Source = './TactileMode/';
Path_Target = './TactileMode/';
SignalDuration = 2; % Signal duration (secs)
Fs = 20000; % Sampling Frequency (Hz)
dt = 1/Fs; % Sampling duration (secs)

begin_len = 0.1*Fs;
win_len = 0.1*Fs;
pad_len = 0.2*Fs;

is_disp = 1;

path_info = dir([Path_Source,'*.mat']);
path_len = length(path_info);

slct_signal_ind = [2, 11, 19, 20, 23, 25, 34];

for i = 1:path_len
    a_name = path_info(i).name;
    disp([a_name,' converted']);

    importSig = load([Path_Source,a_name]);
    importSig = importSig.TimeDomAct;
    t = (0:(length(importSig)-1))*dt; % Time stamps
    
    outSig = importSig(1:win_len);
    
    outSig = [outSig,zeros(1,pad_len)];
    
    outSig = repmat(outSig',[3,1]);
    
    outSig = [zeros(begin_len,1);outSig];

    SN = length(outSig); % Sample Number
    RF = Fs/length(outSig); % Repeating Frequency (Hz)

    % % % % Output signal as .txt file
    file_id = fopen(sprintf('%s%s_SN-%d_RF-%.6f.txt',...
        Path_Target,a_name,SN,RF),'w');
    fprintf(file_id, 'ArbWave = \n');
    fprintf(file_id, '%f\n', outSig);
    fclose(file_id);
end

