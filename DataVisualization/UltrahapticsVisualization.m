%% Visualize PDV data 
% Created on 01/14/2019 based on 'FrequencyBPNoiseCompare.m'
% -------------------------------------------------------------------------
% clear all
% -------------------------------------------------------------------------
% dataName = 'Greg_MovingSpot_1ms_Dir1_1';
dataName = 'Greg_MovingSpot_11ms_Dir1_1';
% -------------------------------------------------------------------------
MapPath = './UltrahapticsMeasurement/';
DataPath = sprintf('../Data_Ultrahaptics/%s.svd',dataName);

cleanDataPath = sprintf('../Data_Ultrahaptics/%s.mat',dataName(1:end-2));
load(cleanDataPath);

Alpha = 25.5; % (mm)
% Alpha = 32;
C = 0.4;
px2mm = 0.289;  % (mm/pixel)
interp_radius = 100;
maskThreshold = 100;

Fs = 125000;
% -------------------------------------------------------------------------
if ~exist('data_info','var')
    [t,y,data_info] = GetPointData(DataPath, 'Time', 'Vib', 'Velocity',...
        'Samples', 0, 0);
    
    [t_ref,ref,~] = GetPointData(DataPath, 'Time', 'Ref1', 'Voltage',...
        'Samples', 0, 0);
end

locator_num = size(y,1);
for i = 1:locator_num
    try
        pt_i = GetIndexOfPoint(DataPath,i);
    catch
        warning(sprintf('Error occurred reading point %d',i));
    end
    
    if (pt_i ~= i)
        warning(sprintf('Measurement point index mismatch! %d ~= %d',...
            pt_i,i));
    end
end

% maskImg = imread([MapPath,sprintf('%s_Mask.jpg',dataName)]);
maskImg = imread([MapPath,'Greg_MovingSpot_1ms_Dir1_1_Mask.jpg']);
maskImg = rgb2gray(maskImg);
% pointImg = imread([MapPath,sprintf('%s_MP.jpg',dataName)]);
pointImg = imread([MapPath,'Greg_MovingSpot_1ms_Dir1_1_MP.jpg']);
pointImg = rgb2gray(pointImg); 

% Low pass filter the chosen grid and computable area
lpF = ones(3)/9; % 3-by-3 mean filter
ptImg_PFilter = uint8(filter2(lpF,pointImg));
MP_Posi = findMP(ptImg_PFilter, locator_num, 128);
MP_Posi = round(MP_Posi); % Round the position to get index
xDiff = diff(MP_Posi(:,2));
xDiff(abs(xDiff)<5) = 0;
MP_Score = cumsum([0; xDiff]).*10000 + MP_Posi(:,1);
[~,ind] = sort(MP_Score);
MP_Posi = MP_Posi(ind,:);

% figure('Position',[60,60,1840,880],'Color','w');
% imshow(imread([MapPath,'Greg_MovingSpot_1ms_Dir1_1_MP.jpg']))
% for i = 1:locator_num
%     text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','r')
% end

%%
allXYPosi = allPositions;
remain_ind = zeros(locator_num,1);
for i = 1:locator_num
    posi_diff = bsxfun(@minus,positions,allXYPosi(i,:)); 
    remain_ind(i) = (sum(posi_diff(:) == 0) == 3);
end
remain_ind = logical(remain_ind);

%% Pre-filtering and Synchronization
% y_filt = y;
% 
% start_i = 900;
% truncLen = 27400 - start_i;
% 
% truncData = NaN(locator_num,truncLen);
% truncData(end,:) = y_filt(end,start_i+(1:truncLen));
% % plot(truncData(end,:)); hold on;
% for i = 1:locator_num-1
%     [acor,lags] = xcorr(ref(end,:),ref(i,:));
%     [~,max_i] = max(abs(acor));
%     lagDiff = lags(max_i);
%     truncData(i,:) = y_filt(i,start_i-lagDiff+(1:truncLen));
% %     plot(truncData(i,:));
% end

%%
truncData = 0.5*y_vib_clean1 + 0.5*y_vib_clean2;

avgWinLen = 20;

y_avg = [];
% for i = 0:88
for i = 0:450
% for i = 0:348
% for i = 0:1750
    slct_ind = i*avgWinLen+(1:(2*avgWinLen));
    y_avg = [y_avg;mean(abs(truncData(:,slct_ind)),2)'];
end
% USE WINDOW FUNCTION
% CONSIDER RMS?

% slow_factor = 25;
% slow_factor = 125;
slow_factor = 250;
frame_rate = Fs/avgWinLen/slow_factor;

v_h = VideoWriter(sprintf('%s_slow%dx_%.ffps.avi',dataName,slow_factor,...
    frame_rate));
v_h.FrameRate = frame_rate;
open(v_h);

colorRange = [min(y_avg(:)),max(y_avg(:))*0.8*17/Alpha];
% colorRange = [min(y_avg(:)),max(y_avg(:))*17/Alpha];

curr_fig = figure('Position',[60,60,1840,880],'Color','w');
colormap(jet(1000));

for i = 1:size(y_avg,1)
interpImg = interpMP(maskImg, MP_Posi(remain_ind,:), y_avg(i,:),...
        maskThreshold, interp_radius, Alpha, C, px2mm);      
surf(flipud(interpImg),'EdgeColor','none');
caxis(colorRange);
view(2)
axis equal; axis off;
title(sprintf('Time = %.0f ms',i*avgWinLen*1000/Fs))
drawnow;
writeVideo(v_h,getframe(curr_fig));
end
close(v_h);
