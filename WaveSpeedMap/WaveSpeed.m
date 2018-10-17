%% Visualize wave speed distribution
% Created on 10/16/2018 based on 'CompareDorsalWithVolar2.m'
% -------------------------------------------------------------------------
clear all
% -------------------------------------------------------------------------
sbj_label = 'SBJ1';
% sbj_label = 'SBJ2';

% scan_side = 'Dorsal';
scan_side = 'Volar';

% -------------------------------------------------------------------------
DataPath = sprintf('../Data_TapDigitII/TapDigitII_%s/',sbj_label);

dataName = [DataPath,sprintf('TapDigitII_%s_%s.svd',...
    scan_side,sbj_label)];

if ~exist('data_info','var')
    [t,y,data_info] = GetPointData(dataName, 'Time', 'Vib', 'Velocity',...
        'Samples', 0, 0);
    
    dt = mean(diff(t));
    Fs = 1/dt;
end

if ~exist('XYZ','var')
    XYZ = GetXYZCoordinates(dataName, 0);
end

locator_num = size(y,1);
for i = 1:locator_num
    try
        pt_i = GetIndexOfPoint(dataName,i);
    catch
        warning(sprintf('Error occurred reading point %d',i));
    end
    
    if (pt_i ~= i)
        warning(sprintf('Measurement point index mismatch! %d ~= %d',...
            pt_i,i));
    end
end

% y_tavg = rms(y,2);
% [~,source_ind] = max(y_tavg);

%% Measurement location
% maskImg = imread([DataPath,'Dorsal_Mask.jpg']);
% maskImg = rgb2gray(maskImg);
pointImg = imread([DataPath,sprintf('%s_MP.jpg',scan_side)]);
pointImg = rgb2gray(pointImg); 

% Low pass filter the chosen grid and computable area
lpF = ones(3)/9; % 3-by-3 mean filter
ptImg_PFilter=uint8(filter2(lpF,pointImg));
MP_Posi = findMP(ptImg_PFilter, locator_num, 128);
MP_Posi = round(MP_Posi); % Round the position to get index
xDiff = diff(MP_Posi(:,2));
xDiff(abs(xDiff)<5) = 0;
MP_Score = cumsum([0; xDiff]).*10000 + MP_Posi(:,1);
[~,ind] = sort(MP_Score);
MP_Posi = MP_Posi(ind,:);

if 0 %---------------------------------------------------------------------
figure('Position',[60,20,820,940],'Color','w')
imshow(imread([DataPath,sprintf('%s_MP.jpg',scan_side)]))
for i = 1:locator_num
    text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','r')
end
end %----------------------------------------------------------------------

%% Estimate wave speed
MP_radius = 40; % (pixel)

exclude_ind = []; % Volar
% exclude_ind = 337; % SBJ1 - Dorsal
exclude_ind = 373; % SBJ2 - Dorsal


slct_ind = setdiff(1:locator_num,exclude_ind);
slct_len = length(slct_ind);
MP_slct = MP_Posi(slct_ind,:);
y_slct = y(slct_ind,:)';
XYZ_slct = XYZ(slct_ind,:);


BPF_freq = 75:50:725;
BPF_num = length(BPF_freq) -1;

avg_speed = NaN(BPF_num,1);
for f_i = 1:BPF_num
%     y_BPF = (bandpass(y_slct,[BPF_freq(f_i),BPF_freq(f_i+1)],Fs))';
    bpFilter = designfilt('bandpassfir', ... % Response type
       'FilterOrder',400, ...            % Filter order
       'CutoffFrequency1',BPF_freq(f_i),...
       'CutoffFrequency2',BPF_freq(f_i+1),...
       'SampleRate',Fs);               % Sample rate
    y_BPF = (filtfilt(bpFilter,y_slct))';    
    
    estimate_num = zeros(slct_len,1);
    wavespeed = NaN(slct_len,7);
    
    for i = 1:slct_len
        ind = find((MP_slct(:,1) < (MP_slct(i,1)+MP_radius)) &...
            (MP_slct(:,1) > (MP_slct(i,1)-MP_radius)) &...
            (MP_slct(:,2) < (MP_slct(i,2)+MP_radius)) &...
            (MP_slct(:,2) > (MP_slct(i,2)-MP_radius)));
        ind = ind(ind~=i);
        
        if ~isempty(ind)
            for j = 1:length(ind)
                [r,lags] = xcorr(y_BPF(i,:),y_BPF(ind(j),:));
                [~,max_i] = max(r);
                delay_t = lags(max_i)*dt; % (sec)
                dist_3D = sum((XYZ_slct(i,:) - XYZ_slct(ind(j),:)).^2).^0.5; % (m)
                
                if delay_t ~= 0
                    estimate_num(ind(j)) = estimate_num(ind(j)) +1;
                    wavespeed(ind(j),estimate_num(ind(j))) = dist_3D/delay_t; % (m/sec)
                end
            end
        end
    end
    avg_speed(f_i) = mean(max(abs(wavespeed),[],2),'omitnan');
    fprintf('Band %d to %d Hz - average speed = %.1f m/s\n',...
        BPF_freq(f_i),BPF_freq(f_i+1),avg_speed(f_i));
end

figure('Position',[160,80,720,600],'Color','w')
plot(BPF_freq(2:end),avg_speed,'.-');
xlabel('Frequency (Hz)'); ylabel('Speed (m/s)');
box off; set(gca,'FontSize',16); 
title(sprintf('%s - %s Scan',sbj_label,scan_side));
% ylim([0 30])

%%
% figure('Position',[60,20,820,940],'Color','w')
% for i = 1:slct_len
%     ind = find((MP_slct(:,1) < (MP_slct(i,1)+MP_radius)) &...
%         (MP_slct(:,1) > (MP_slct(i,1)-MP_radius)) &...
%         (MP_slct(:,2) < (MP_slct(i,2)+MP_radius)) &...
%         (MP_slct(:,2) > (MP_slct(i,2)-MP_radius)));
%     ind = ind(ind~=i);
%     imshow(imread([DataPath,sprintf('%s_MP.jpg',scan_side)]))
%     hold on
%     if ~isempty(ind)
%         scatter(MP_slct(i,2),MP_slct(i,1),20,'r','filled')
%         plot(MP_slct(ind,2),MP_slct(ind,1),'or')
%     else
%         scatter(MP_slct(i,2),MP_slct(i,1),20,'g','filled')
%     end
%     pause(0.01)
%     hold off
%     fprintf('MP %d [%d]\n',i,length(ind));
% end




