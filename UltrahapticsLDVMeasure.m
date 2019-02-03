%% Sampling NI Terminal Analog Input in Real Time
% Created on 02/02/2016 based on 'MeasureProgram'
% -------------------------------------------------------------------------
warning ('off','all');
% -------------------------------------------------------------------------
Fs = 44000; % Sampling frequency (44 kHz)
DispDuration = 20; % Display time duration of * seconds
fprintf('Sampling Frequency = %d Hz\nDisplay time duration = %d seconds\n',...
    Fs,DispDuration);

ScaleFactor = 25; % 25(mm/s)/V for 100mm/s measurement scale

LDV_info.LaserPointSize = 16^2; % (px)
LDV_info.LDVHeight = 79.4; % (mm)

AIData = struct('LDV',[],'ref',[],'Fs',Fs,...
    'ScaleFactor',ScaleFactor,'cosTheta',[]);

% -------------------------------------------------------------------------
% Initialize the NI DAQ interface to acquire analog data
if ~exist('sNI','var')  
    % Analog input session ------------------------------------------------
    sNI(1) = daq.createSession('ni');
    sNI(1).IsContinuous = 1;
    sNI(1).Rate = Fs;

    % Analog input channel (LDV)
    chNI0 = sNI(1).addAnalogInputChannel('Dev3','ai0','Voltage');
    chNI0.TerminalConfig = 'SingleEndedNonReferenced';
    chNI0.Range = [-5,5]; % Input voltage range from -5 to 5V
  
    % Analog input channel (reference signal)
    chNI1 = sNI(1).addAnalogInputChannel('Dev3','ai1','Voltage');
    chNI1.TerminalConfig = 'SingleEndedNonReferenced';
    chNI1.Range = [-5,5]; % Input voltage range from -5 to 5V
end

% -------------------------------------------------------------------------
% Display Initialization
figSize = [20,100,1880,800];
fig_h = figure('Name','NI Data Acquisition','Position',figSize,...
    'Color',[0.98,0.98,0.98]);

% Measurement signal display ----------------------------------------------
sigDisp_h(1) = subplot('Position',[0.25 0.56 0.72 0.4]);
plt_h(1) = plot(0,0);
ylim([-4.1 4.1]); % Dsiplay voltage range from -4.1 to 4.1V
ylabel('LDV Amplitude (V)')
set(plt_h(1), 'UserData', AIData);

sigDisp_h(2) = subplot('Position',[0.25 0.06 0.72 0.4]);
hold on
plt_h(2) = plot(0,0);
ylim([0.0 0.2]); 
xlim([0 DispDuration])
xlabel('Time (secs)');
ylabel('Ref Amplitude (V)')

% % % plt_h(3) = ; % Reserved for more AI channels
% % % plt_h(4) = ; % Reserved for more AI channels
% % % plt_h(5) = ; % Reserved for more AI channels

% Textbox for message
plt_h(6) = uicontrol('Style','text','BackgroundColor','w',...
    'Position',[36 figSize(4)-420 380 30],...
    'String',[],'FontSize',14);

% Start button
uicontrol('Style', 'pushbutton', 'String', 'START SAMPLING',...
    'Position', [30 figSize(4)-100 100 40],...
    'Callback', {@startSamp, plt_h, sNI, sigDisp_h, DispDuration},...
    'BackgroundColor',[0.8,1,0.8]);

% Stop button
uicontrol('Style', 'pushbutton', 'String', 'STOP SAMPLING',...
    'Position', [30 figSize(4)-160 100 40],...
    'Callback', {@stopSamp, plt_h, sNI},...   
    'BackgroundColor',[1,0.8,0.8]);

% Analysis button
uicontrol('Style', 'pushbutton', 'String', 'LDV ANALYSIS',...
    'Position', [30 figSize(4)-220 100 40],...
    'Callback', {@analyzeData, plt_h, sNI, Fs},...
    'BackgroundColor',[0.8,1,1]);

% Save button
uicontrol('Style', 'pushbutton', 'String', 'SAVE DATA',...
    'Position', [30 figSize(4)-280 100 40],...
    'Callback', {@saveData, plt_h, sNI},...
    'BackgroundColor',[1,1,0.8]);

% Save CSV button
uicontrol('Style', 'pushbutton', 'String', 'SAVE CSV',...
    'Position', [30 figSize(4)-340 100 40],...
    'Callback', {@saveCSV, plt_h, sNI},...
    'BackgroundColor',[1,1,0.75]);

% Event listener
lh(1) = addlistener(sNI(1),'DataAvailable',...
    @(src,event) streamingAnalogData(src,event, plt_h));

% Figure close button (end the program)
set(fig_h, 'CloseRequestFcn',{@closeReq, sNI, fig_h, lh});

% -------------------------------------------------------------------------
% Main loop: clean display and save data 
disp_length = DispDuration*Fs;
while ishandle(fig_h)
    curr_time = get(plt_h(1),'XData'); % Time stamps of current display

    if length(curr_time) > disp_length
        curr_data = get(plt_h(1),'YData'); % Data of current display 1

        curr_data2 = get(plt_h(2),'YData'); % Data of current display 2 
    
        AIData = get(plt_h(1), 'UserData');
        AIData.LDV = [AIData.LDV, [curr_time;curr_data]];
        AIData.ref = [AIData.ref, [curr_time;curr_data2]];
        set(plt_h(1), 'UserData', AIData);
        
        set(plt_h(1), 'XData',[], 'YData',[]);
        set(plt_h(2), 'XData',[], 'YData',[]);
        xlim(sigDisp_h(1),[curr_time(end) curr_time(end)+DispDuration]);       
        xlim(sigDisp_h(2),[curr_time(end) curr_time(end)+DispDuration]);
    end  
    
    pause(0.001)
end

% =========================================================================
%% Even Listener function: Plot the streaming data
function streamingAnalogData(~, event, plt_h)
    t_stamps = get(plt_h(1),'XData');
    e_data = get(plt_h(1),'YData'); % LDV
    t_stamps = [t_stamps,event.TimeStamps'];
    e_data = [e_data,event.Data(:,1)'];
    set(plt_h(1), 'XData',t_stamps, 'YData',e_data);
    
    t_stamps2 = get(plt_h(2),'XData');
    e_data2 = get(plt_h(2),'YData'); % Ref
    t_stamps2 = [t_stamps2,event.TimeStamps'];
    e_data2 = [e_data2,event.Data(:,2)'];
    set(plt_h(2), 'XData',t_stamps2, 'YData',e_data2);
end

%% Local function: Start sampling button
function startSamp(~, ~, p_h, s_h, sigDisp_h, DispDuration)
    fprintf('Sampling...\n');

    if ~s_h(1).IsLogging
        p_h(6).String = 'Starting A New Measurement';

        AIData = get(p_h(1),'UserData');
        AIData.cosTheta = 1;
        AIData.LDV = [];
        AIData.ref = [];
        set(p_h(1), 'XData',[], 'YData',[],'UserData',AIData);
        set(p_h(2), 'XData',[], 'YData',[]);
        xlim(sigDisp_h(1),[0 DispDuration]);
        xlim(sigDisp_h(2),[0 DispDuration]);
        
        s_h(1).startBackground();

        tic
    end
end
    
%% Local function: Stop sampling button
function stopSamp(~, ~, p_h, s_h)
    s_h(1).stop;
    p_h(6).String = sprintf('%.1f seconds of Measurement Completed\n', toc); 
    disp(p_h(6).String);

    pause(0.01);
    AIData = get(p_h(1), 'UserData');
    curr_t = get(p_h(1),'XData'); % Time stamps of current display
    curr_d = get(p_h(1),'YData'); % Data of current display 1
    curr_d2 = get(p_h(2),'YData'); % Data of current display 2
        
    AIData.LDV = [AIData.LDV, [curr_t;curr_d]];
    AIData.ref = [AIData.ref, [curr_t;curr_d2]];
    set(p_h(1), 'UserData', AIData);
end

%% Local function: Analyze LDV data
function analyzeData(~, ~, p_h, s_h, samp_freq)
    if  s_h(1).IsLogging
        msgbox('Stop sampling before analyzing the data','Failed','warn')
    else
        s_data = get(p_h(1), 'UserData');
        AIData = s_data.LDV';
        
        disp('Analyze data: DC removed by average subtraction!')
        AIData(:,2) = AIData(:,2) -mean(AIData(:,2));
        
        figure('Position',[40,40,1800,960])
        subplot(2,1,1)
        plot(AIData(:,1),AIData(:,2));
        title('Temporal Domain')
        xlabel('Time (secs)');
        ylabel('Amplitude (V)');
        
        [abs_Y_fft,freq]=powerSpectrum(AIData(:,2), samp_freq); 
        mean_Y_fft = conv(abs_Y_fft, ones(3,1)/3, 'same');
        [~,descend_ind] = sort(abs_Y_fft,'descend');
        peak_freq = mean(freq(descend_ind(1:10)));
        abs_Y_fft = 20*log10(abs_Y_fft);
        mean_Y_fft = 20*log10(mean_Y_fft);
        subplot(2,1,2)
        plot(freq,abs_Y_fft)
        hold on
        plot(freq,mean_Y_fft)
        plot([peak_freq,peak_freq],ylim,'g');
        text(peak_freq,-10,sprintf('%.3f Hz',peak_freq));
        hold off
        title('Frequency Domain')
        xlabel('Frequency (Hz)');
        ylabel('Amplitude (dB)');
    end
end

%% Local function: Save data as .mat file
function saveData(~, ~, p_h, s_h)
    if  s_h(1).IsLogging
        msgbox('Stop sampling before saving','Failed','warn')
    else
        AIData = get(p_h(1), 'UserData');
        LDVvel = AIData.ScaleFactor*single(AIData.LDV(2,:)); % Converted to single precision
        refSig = single(AIData.ref(2,:)); % Converted to single precision

        figure('Position',[40,80,1600,900]);
        subplot(2,1,1)
        plot(LDVvel); ylabel('LDV (mm/s)'); box off;
        subplot(2,1,2)
        plot(refSig); ylabel('Ref (V)'); box off;
        xlabel('Sample')
        
        Fs = AIData.Fs;
        cosTheta = AIData.cosTheta;
        file_name = 'TouchAmp_measure';
        save(file_name,'LDVvel','refSig','Fs','cosTheta');
        p_h(6).String = [file_name,' Saved'];
        disp(p_h(6).String);
    end
end

%% Local function: Output data to .csv file
function saveCSV(~, ~, p_h, s_h)
    if  s_h(1).IsLogging
        msgbox('Stop sampling before output','Failed','warn')
    else
        file_name = inputdlg('File name:','Save as .csv file',...
            [1 60],{'AIData'});
        if ~isempty(file_name)
            AIData = get(p_h(1), 'UserData');
            
            dataTable = [AIData.LDV(1,:)',...
                AIData.ScaleFactor*AIData.LDV(2,:)',...
                AIData.ref(2,:)'];
            dataTable = array2table(dataTable,'VariableNames',...
                {'t','LDV','ref'});
            
            writetable(dataTable,sprintf('%s.csv',file_name{1}),...
                'Delimiter',',','QuoteStrings',true)

            p_h(6).String = 'CSV Data Saved';
            disp(p_h(6).String);
        else
            disp('Failed to save to .csv file!');
        end
    end
end

%% Local function: Figure Close button
function closeReq(~, ~, s_h, f_h, l_h)
    if  s_h(1).IsLogging
        disp('Sampling ceased! Data lost!')
    end
    delete(l_h);
    pause(0.01)
    delete(f_h);
    s_h(1).stop;
    disp('---------- Program forced shutdown ----------')
end

%% Local function: FFT - Spectrum
function [ abs_Y_fft, f ] = powerSpectrum( input_data, Fs )
% FFT - Spectrum
% frequenc range (optional) must be a 1-by-2 vector 
% Author: Yitian Shao
% Adopted as local function on 11/07/2016
%==========================================================================
L=length(input_data); % Sampling data number
NFFT = 2^nextpow2(L); % Next power of 2 from length of y

Y_fft = fft(input_data, NFFT)/L; % FFT 
f = Fs/2*linspace(0, 1, NFFT/2+1); % Frequency scale

abs_Y_fft = 2*abs(Y_fft(1:NFFT/2+1)); % Power spectrum 
end