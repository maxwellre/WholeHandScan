%% Touch Amp Signal Preparation
Path_Source = './out2/';
Path_Target = './TouchAmpSignal/';
SignalDuration = 2; % Signal duration (secs)
Fs = 44100; % Sampling Frequency (Hz)
dt = 1/Fs; % Sampling duration (secs)

outFs = 50000; % Output Sampling Frequency (Hz)

begin_len = 0.1*outFs;
pad_len = 0.2*outFs;

is_disp = 1;

path_info = dir([Path_Source,'*.csv']);
path_len = length(path_info);

% t_ind = 1:SignalDuration*Fs;

% slct_signal_ind = [2, 11, 19, 20, 23, 25, 34];
% slct_signal_ind = 24; t_ind = 16601:22000;

slct_signal_ind = 12; t_ind = 16801:36600;
% slct_signal_ind = 25; t_ind = 16801:36600;

allSig = zeros(begin_len,1);
for i = 1:path_len
    a_name = path_info(i).name;
    num_ind = 1:(strfind(a_name,'.csv')-1);
    signal_i = str2double(a_name(num_ind));
    if ~isempty(find(slct_signal_ind == signal_i,1))
        disp([a_name,' converted']);
        
        importSig = csvread([Path_Source,a_name]);
        importSig = importSig(t_ind);
        
        t = (0:(length(importSig)-1))*dt; % Time stamps
        
        outSig = resample(importSig,t,outFs);
        t_out = (0:(length(outSig)-1))/outFs;
        
        if is_disp
            plot(t,importSig);
            hold on
            plot(t_out,outSig,'--');
            hold off
        end
        
%         allSig = [allSig;outSig];
        outSig = [zeros(begin_len,1);outSig;zeros(pad_len,1)];
        
        SN = length(outSig); % Sample Number
        RF = outFs/length(outSig); % Repeating Frequency (Hz)
        
        % % % % % Output signal as .txt file
        file_id = fopen(sprintf('%sSig%02d_SN-%d_RF-%.6f.txt',...
            Path_Target,signal_i,SN,RF),'w');
        fprintf(file_id, 'ArbWave = \n');
        fprintf(file_id, '%f\n', outSig);
        fclose(file_id);
    end
end

% SN = length(allSig); % Sample Number
% RF = outFs/length(allSig); % Repeating Frequency (Hz)
% file_id = fopen(sprintf('%sAllSig_SN-%d_RF-%.6f.txt',...
%             Path_Target,SN,RF),'w');
% fprintf(file_id, 'ArbWave = \n');
% fprintf(file_id, '%f\n', allSig);
% fclose(file_id);