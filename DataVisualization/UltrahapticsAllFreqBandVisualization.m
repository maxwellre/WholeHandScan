%% Visualize PDV data 
% Created on 01/29/2019 based on 'UltrahapticsVisualization.m'
% -------------------------------------------------------------------------
% clear all
% -------------------------------------------------------------------------
slctData = 3; 

switch slctData
    case 1
        dataName = 'Greg_MovingSpot_1ms_Dir1';
    case 2
        dataName = 'Greg_MovingSpot_1ms_Dir2';
    case 3
        dataName = 'Greg_MovingSpot_11ms_Dir1';
    case 4
        dataName = 'Greg_MovingSpot_11ms_Dir2';
end

% -------------------------------------------------------------------------
MapPath = './UltrahapticsMeasurement/';
DataPath = sprintf('../Data_Ultrahaptics/%s_1.svd',dataName);
OutputPath = './UltrahapticsVideo/';

cleanDataPath1 = sprintf('../Data_Ultrahaptics/%s_1.mat',dataName);
% cleanDataPath1 = sprintf('../Data_Ultrahaptics/%s_2.mat',dataName);
load(cleanDataPath1);

Alpha = 25.5; % (mm)
% Alpha = 32; % (mm)

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
    
     XYZ = GetXYZCoordinates(DataPath, 0);
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

%% Configuration
% slow_factor = 25;
slow_factor = 125;
% slow_factor = 250;

% -------------------------------------------------------------------------
% % % Sequence: 20,40,80,120,160,200,240,280,320,360,400,440,480,520,560,600,640 Hz
% -------------------------------------------------------------------------
% freqBand = [20,40,80,160,240,320,400,480,560,640];
freqBand = [40,240];

freqBandNum = length(freqBand)-1;

% truncData = 0.5*y_vib_clean1 + 0.5*y_vib_clean2;

%% Discard ill measurement points
% dataEn = rms(y_vib_sync');
% sorted_En = sort(dataEn);
% threshold = sorted_En(ceil(0.92*length(sorted_En)));
% % threshold = sorted_En(ceil(0.88*length(sorted_En)));
% % threshold = sorted_En(ceil(0.80*length(sorted_En)));
% remain_ind = (dataEn < threshold);

discard_ind = [19, 20, 38, 39, 48, 49, 58, 68, 69, 79, 80, 176, 190, 193, 205, 216,...
    220, 222, 225, 228, 233, 247, 283, 289];

if slctData
    discard_ind = [discard_ind, 2,26,37,67,74, 94, 106,144, 219,223, 242, 281];
else
    discard_ind = [discard_ind, 8, 99, 120, 138, 144, 158, 195, 200, 277];
end

remain_ind = true(1,locator_num);
remain_ind(discard_ind) = false;

remainMP_Posi = MP_Posi(remain_ind,:);

%% Display selected measurement points and provide distance measure
temp_fig = figure('Position',[60,60,1840,880],'Color','w');
imshow(imread([MapPath,'Greg_MovingSpot_1ms_Dir1_1_MP.jpg']))
hold on;
for i = find(remain_ind > 0)
    text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','c')
end
for i = find(remain_ind == 0)
    scatter(MP_Posi(i,2),MP_Posi(i,1),18,'r','filled');
    text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','y')
end
title(sprintf('%d Point Removed',sum(remain_ind == 0)));

MP_Dist = ((MP_Posi(2:end,1)-MP_Posi(1,1)).^2+...
    (MP_Posi(2:end,2)-MP_Posi(1,2)).^2).^0.5;
MP_realDist = ((XYZ(2:end,1)-XYZ(1,1)).^2 +...
    (XYZ(2:end,2)-XYZ(1,2)).^2 +...
    (XYZ(2:end,3)-XYZ(1,3)).^2).^0.5; % Unit:m
realToMap = median(MP_Dist./MP_realDist);

fprintf('Real distance to image distance ratio = %f pixel/m\n',realToMap);

% scatter(MP_Posi(remain_ind,2),MP_Posi(remain_ind,1),20,y_avg(1,:),'filled');
hold off;

close(temp_fig);

%% Moving Average Windowing
avgWinLen = 20;
hammWin = repmat(hamming(2*avgWinLen),[1,sum(remain_ind)]);

for f_i = 1:freqBandNum
filteredData = bandpass(y_vib_sync(remain_ind,:)',...
    [freqBand(f_i),freqBand(f_i+1)],Fs);

% filteredData = y_vib_sync(remain_ind,:)';
% filteredData = highpass(y_vib_sync(remain_ind,:)',40,Fs);
% filteredData = lowpass(filteredData,240,Fs);

% plot(filteredData(1:28000,:));

y_rect = [];
y_avg = [];

if slctData
    frame_num = 1200;
else
    frame_num = 120;
end

for i = 0:frame_num
% for i = 801:900
    slct_ind = i*avgWinLen+(1:(2*avgWinLen));
    y_rect = [y_rect;rms(hammWin.*filteredData(slct_ind,:))]; % RMS Rectify
    y_avg = [y_avg;mean(hammWin.*filteredData(slct_ind,:))]; % (Bipolar)
end

% y_avg = y_rect;

%% Estimate focus point path
% Image coordinates
[meshX,meshY] = meshgrid(1:size(maskImg,2),1:size(maskImg,1));
switch slctData
    case 1
    case 2
    case 3
        estim_frame = 32:114;
        focus_speed = 11*realToMap; % (pixel/sec);
    case 4
end
frame_num = length(estim_frame);
focus_Posi = NaN(frame_num,2);
wb_h = waitbar(0, 'O', 'Name','Estimating focus point path ...');
for i = 1:frame_num
    waitbar(i/frame_num,wb_h,sprintf('Frame %d',i));
    % Track the focus point
    trackImg = interpMP(maskImg, MP_Posi(remain_ind,:),...
        y_rect(estim_frame(i),:), maskThreshold, interp_radius, Alpha,...
        C, px2mm); 
    valid_ind = ~isnan(trackImg);
    Phi = exp(200000*trackImg(valid_ind))';
    
    switch slctData
        case 1
            Phi(Phi < 2e7) = 0;
        case 3
            Phi(Phi < 2e7) = 0;
    end

    if sum(Phi > 0)
        Phi = Phi./sum(Phi);
        focus_Posi(i,2) = Phi*meshX(valid_ind);
        focus_Posi(i,1) = Phi*meshY(valid_ind);
    end
end
close(wb_h);
valid_ind = ~isnan(focus_Posi(:,1));
focus_Posi = focus_Posi(valid_ind,:);
pf_h = polyfit(focus_Posi(:,1),focus_Posi(:,2),1);
% focus_endPoints = [min(focus_Posi(:,1)),max(focus_Posi(:,1))];
focus_endPoints = [focus_Posi(1,1),focus_Posi(end,1)];
focus_fit2 = polyval(pf_h,focus_endPoints); 
focus_endPoints = [focus_endPoints;focus_fit2];
imshow(pointImg);
hold on;
plot(focus_endPoints(2,:),focus_endPoints(1,:),'-or')

pathCosine = abs(focus_endPoints(1,2) - focus_endPoints(1,1))./...
    (((focus_endPoints(1,2) - focus_endPoints(1,1)).^2+...
    (focus_endPoints(2,2) - focus_endPoints(2,1)).^2).^0.5);

%% Plot selected frames
if 1 %---------------------------------------------------------------switch
row_num = 3;

switch slctData
    case 1
        slct_frame = (0:100:990)+50+2; % (1m/s included in the paper)
        % slct_frame = (0:100:900)+10+2; % (not good)
        % slct_frame = 1:5:size(y_avg,1); % (test only);
    case 3
        frame_interval = 9;
%         slct_frame = -1+(24:9:106);
        slct_frame = 32:frame_interval:114;
end


frame_num = length(slct_frame);

% focus_points = linspace(focus_endPoints(1,1),focus_endPoints(1,end),...
%     frame_num);

% (pixel/sec)*(cos[Theta])*(sample/frame)/(1/sec)*(frame indcies)
focus_points = focus_endPoints(1,1) +...
    (focus_speed*pathCosine*frame_interval*avgWinLen/Fs).*(0:frame_num-1);

focus_fit2 = polyval(pf_h,focus_points); 
focus_points = [focus_points;focus_fit2]';

% plot(focus_points(:,2),focus_points(:,1),'*b')

y_slct = y_avg(slct_frame,:);
% colorRange = [min(y_slct(:)),max(y_slct(:))*0.7*17/Alpha];
% colorRange = [-0.0002, 0.0002];
% colorRange = [0, 0.00016];
% colorRange = [-0.00016, 0.00016];
colorRange = [-0.0001, 0.0001];

curr_fig = figure('Position',[60,60,1840,580],'Color','w','Name',...
    sprintf('Bandpass [%d - %d Hz]',freqBand(f_i),freqBand(f_i+1)));

colormap(jet(1000));

switch slctData
    case 1
        subplot_width = 0.098;
    case 3
        subplot_width = 0.098;
end

for i = 1:frame_num
%     subplot(row_num,ceil(frame_num/row_num),i)
    subplot('Position',[0.005+(i-1)*subplot_width,0.2,subplot_width,0.6]);  
    
    interpImg = interpMP(maskImg, MP_Posi(remain_ind,:),...
      y_slct(i,:),maskThreshold, interp_radius, Alpha, C, px2mm);      
    sc_h = imagesc(interpImg); 
    set(sc_h,'AlphaData',~isnan(interpImg));
    caxis(colorRange);
    view(2)
    axis equal; axis off;
    
    hold on
    scatter(focus_points(i,2),focus_points(i,1),30,'w','Linewidth',1);
    hold off
    
%     title(sprintf('(%d) Time = %.0f ms',slct_frame(i),...
%         slct_frame(i)*avgWinLen*1000/Fs))
    title(sprintf('%.2f ms',...
        (slct_frame(i)-slct_frame(1))*avgWinLen*1000/Fs));
    drawnow;
    
    set(gca,'FontSize',16)   
    if i == round(0.5*frame_num)
        c_h = colorbar('Location','south','Ticks',[]);
        c_h.Label.String = sprintf('RMS Velocity %.1f - %.1f (mm/s)\n',...
            1000*colorRange);
        disp(c_h.Label.String{1});
    end
end

end %------------------------------------------------------------switch end 

%% Produce Videos
if 0 %---------------------------------------------------------------switch
    
% Image coordinates
[meshX,meshY] = meshgrid(1:size(maskImg,2),1:size(maskImg,1));
    
frame_rate = Fs/avgWinLen/slow_factor;

v_h = VideoWriter(sprintf('%s%s_slow%dx_%.ffps_%d-%dHz.avi',OutputPath,...
    dataName,slow_factor,frame_rate,freqBand(f_i),freqBand(f_i+1)));
v_h.FrameRate = frame_rate;
open(v_h);

% colorRange = [min(y_avg(:)),max(y_avg(:))*0.8*17/Alpha];

if is1ms
    % colorRange = [0, 0.00065];
    colorRange = [-0.0002, 0.0002];
else
%     colorRange = [0, 0.00022];
    colorRange = [-0.0002, 0.0002];
end

curr_fig = figure('Position',[60,60,1840,880],'Color','w');
colormap(jet(1000));

focus_Posi = NaN(size(y_avg,1),2);
for i = 1:size(y_avg,1)
interpImg = interpMP(maskImg, MP_Posi(remain_ind,:), y_avg(i,:),...
        maskThreshold, interp_radius, Alpha, C, px2mm);          
sc_h = imagesc(interpImg); 
set(sc_h,'AlphaData',~isnan(interpImg));
caxis(colorRange);
view(2)
axis equal; axis off;
title(sprintf('(%d) [%d - %d Hz] Time = %.1f ms',i,...
    freqBand(f_i),freqBand(f_i+1),i*avgWinLen*1000/Fs))

% Track the focus point
trackImg = interpMP(maskImg, MP_Posi(remain_ind,:), y_rect(i,:),...
        maskThreshold, interp_radius, Alpha, C, px2mm); 
valid_ind = ~isnan(trackImg);
Phi = exp(200000*trackImg(valid_ind))';
if is1ms
    Phi(Phi < 2e7) = 0;
else
    Phi(Phi < 2e7) = 0;
end

if sum(Phi > 0)
    Phi = Phi./sum(Phi);
    focus_Posi(i,2) = Phi*meshX(valid_ind);
    focus_Posi(i,1) = Phi*meshY(valid_ind);
else
    focus_Posi(i,:) = [-100 -100]; 
end

hold on
scatter(focus_Posi(i,2), focus_Posi(i,1),1200,'k','Linewidth',1);
hold off

drawnow;
writeVideo(v_h,getframe(curr_fig));
end
close(v_h);
pause(0.5);
close(curr_fig);
end %------------------------------------------------------------switch end

end 