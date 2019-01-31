%% Visualize PDV data 
% Created on 01/29/2019 based on 'UltrahapticsVisualization.m'
% -------------------------------------------------------------------------
% clear all
% -------------------------------------------------------------------------
% dataName = 'Greg_MovingSpot_1ms_Dir1';
dataName = 'Greg_MovingSpot_11ms_Dir1';
% -------------------------------------------------------------------------
MapPath = './UltrahapticsMeasurement/';
DataPath = sprintf('../Data_Ultrahaptics/%s_1.svd',dataName);
OutputPath = './UltrahapticsVideo/';

cleanDataPath1 = sprintf('../Data_Ultrahaptics/%s_1.mat',dataName);
load(cleanDataPath1);
% cleanDataPath2 = sprintf('../Data_Ultrahaptics/%s_2.mat',dataName);
% load(cleanDataPath2);

Alpha = 25.5; % (mm)
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

maskImg = imread([MapPath,'Greg_MovingSpot_1ms_Dir1_1_Mask.jpg']);
maskImg = rgb2gray(maskImg);
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

%%
% slow_factor = 25;
slow_factor = 125;
% slow_factor = 250;

% Sequence 20,40,80,120,160,200,240,280,320,360,400,440,480,520,560,600,640 Hz
freqBand = [20,40,80,160,240,320,400,480,560,640];
freqBandNum = length(freqBand)-1;

% truncData = 0.5*y_vib_clean1 + 0.5*y_vib_clean2;

dataEn = rms(y_vib_sync');
sorted_En = sort(dataEn);
threshold = sorted_En(ceil(0.92*length(sorted_En)));
% threshold = sorted_En(ceil(0.88*length(sorted_En)));
remain_ind = (dataEn < threshold);

temp_fig = figure('Position',[60,60,1840,880],'Color','w');
imshow(imread([MapPath,'Greg_MovingSpot_1ms_Dir1_1_MP.jpg']))
hold on;
for i = find(remain_ind > 0)
    text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','c')
end
for i = find(remain_ind == 0)
    scatter(MP_Posi(i,2),MP_Posi(i,1),18,'r','filled');
end
title(sprintf('%d Point Removed',sum(remain_ind == 0)));
hold off;
close(temp_fig);

avgWinLen = 20;
hammWin = repmat(hamming(2*avgWinLen),[1,sum(remain_ind)]);

for f_i = 1:freqBandNum
filteredData = bandpass(y_vib_sync(remain_ind,:)',...
    [freqBand(f_i),freqBand(f_i+1)],Fs);

y_avg = [];
for i = 0:100
% for i = 0:1310
    slct_ind = i*avgWinLen+(1:(2*avgWinLen));
    y_avg = [y_avg;rms(hammWin.*filteredData(slct_ind,:))];
end

frame_rate = Fs/avgWinLen/slow_factor;

v_h = VideoWriter(sprintf('%s%s_slow%dx_%.ffps_%d-%dHz.avi',OutputPath,...
    dataName,slow_factor,frame_rate,freqBand(f_i),freqBand(f_i+1)));
v_h.FrameRate = frame_rate;
open(v_h);

colorRange = [min(y_avg(:)),max(y_avg(:))*0.7*17/Alpha];

curr_fig = figure('Position',[60,60,1840,880],'Color','w');
colormap(jet(1000));

for i = 1:size(y_avg,1)
interpImg = interpMP(maskImg, MP_Posi(remain_ind,:), y_avg(i,:),...
        maskThreshold, interp_radius, Alpha, C, px2mm);      
surf(flipud(interpImg),'EdgeColor','none');
caxis(colorRange);
view(2)
axis equal; axis off;
title(sprintf('[%d - %d Hz] Time = %.0f ms',freqBand(f_i),freqBand(f_i+1),...
    i*avgWinLen*1000/Fs))
drawnow;
writeVideo(v_h,getframe(curr_fig));
end
close(v_h);
pause(0.5);
close(curr_fig);
end