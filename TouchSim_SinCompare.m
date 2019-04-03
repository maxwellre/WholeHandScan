%% TouchSim Simulation of the SLDV data 
% Created on 04/02/2019 based on 'FrequencySinCompare.m'
% -------------------------------------------------------------------------
MapPath = './DataVisualization/DigitIISinusoid/';
DataPath = './Data_DigitIISinusoid/SubjectBha_Session04.svd';
TouchSimPath = '../SurfaceWaveStudyVI/Lib_TouchSim';
% -------------------------------------------------------------------------
addpath('./DataVisualization/');
addpath([TouchSimPath '/base/model']);      % base model code
addpath([TouchSimPath '/base/internal']);   % internal functions
addpath([TouchSimPath '/base/GUI']);        % GUI
addpath([TouchSimPath '/base/skinmech/']);  % Skin mechanics code
addpath([TouchSimPath '/docs']);            % Toolbox documentation
% -------------------------------------------------------------------------
Fs = 20000; % (Sampling frequency = 20 kHz)

% % Create a PC afferent
% a = Afferent('PC','location',[0 0],'idx',1);
% t = 0:(1/Fs):2; trace = 5*(1+sin(2*pi*240.*t))';
% s = Stimulus(trace,[0 0],Fs,0.05); % (~0.1 mm diameter "pin")
% % Calculate response
% r = a.response(s); r.rate

% -------------------------------------------------------------------------
if ~exist('data_info','var')
    [t,y,data_info] = GetPointData(DataPath, 'Time', 'Vib', 'Velocity',...
        'Samples', 0, 0);
    
    [t_ref,ref,~] = GetPointData(DataPath, 'Time', 'Ref1', 'Voltage',...
        'Samples', 0, 0);
end
fprintf('y: %s (Unit: %s)\n',data_info.YName,data_info.YUnit);

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

maskImg = imread([MapPath,'DigitIISinusoid_Mask.jpg']);
maskImg = rgb2gray(maskImg);
pointImg = imread([MapPath,'DigitIISinusoid_MP.jpg']);
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

y_filt = y;

% -------------------------------------------------------------------------
sinNum = 17;
tavgData = cell(sinNum,2);

zero_ind = (abs(ref(1,:)) < 5e-4);
onsite_ind = movmedian(double(zero_ind),19);

onsite_ind = [0,diff(onsite_ind)];
onsite_ind = (onsite_ind == -1);
temp_fig = figure('Position',[60,60,1840,880],'Color','w');
subplot(3,1,1)
plot(onsite_ind);
hold on
plot(ref(1,:)./max(ref(1,:)));
plot(y_filt(end,:)./max(abs(y_filt(end,:))));
ylabel('Reference Signal')

seg_ind = find(onsite_ind == 1);
seg_ind = [seg_ind, (2*seg_ind(end)-seg_ind(end-1))];
hold off

subplot(3,1,2)
plot(t,y_filt(end,:));
xlabel('Time (secs)'); ylabel('Velocity (m/s)')

subplot(3,1,3)
hold on
for i = 1:sinNum
    slct_ind = seg_ind(i):seg_ind(i+1);
    [ seg_spect, f ] = spectr(ref(1,slct_ind), Fs, [0 800]);
    plot(f,seg_spect); xlabel('Frequency (Hz)'); 
    tavgData{i,2} =...
        sprintf('%.0fHz', 10*round(0.1*f*(seg_spect./sum(seg_spect))));
    tavgData{i,1} = rms(y_filt(:,slct_ind),2)';
    
    fprintf('%s maximum amplitude = %f (mm/s)\n',tavgData{i,2},...
        1000*max(tavgData{i,1}));
end
legend(tavgData(:,2));
hold off

close(temp_fig);

%% 
a = Afferent('PC','location',[0 0],'idx',1);

PC_spikes = cell(sinNum,locator_num);
PC_resp = nan(sinNum,locator_num);
PC_duration = nan(sinNum,locator_num);

wb_h = waitbar(0, 'O', 'Name','Estimating afferent response ...');
for i = 1:sinNum
    slct_ind = seg_ind(i):seg_ind(i+1);
    skinVel = y_filt(:,slct_ind)';
    
    skinDisp = cumsum(skinVel)./Fs;
    skinDisp = detrend(skinDisp);
    skinDisp = skinDisp.*1000; % Convert unit from m to mm
    skinDisp = bsxfun(@minus, skinDisp,min(skinDisp));
    
    for j = 1:locator_num
        waitbar(j/locator_num,wb_h,sprintf('(%d) %s',i,tavgData{i,2}));
        s = Stimulus(skinDisp(:,j),[0 0],Fs,0.05); % (~0.1 mm diameter "pin")
        r = a.response(s); 
        
        PC_spikes{i,j} = r.spikes;
        PC_resp(i,j) = r.rate;
        PC_duration(i,j) = r.duration;
    end
end
close(wb_h);

%% Raster plot
ref_posi = [1106 438]; % Fingertip of digitII

dist = ((MP_Posi(:,2) - ref_posi(1)).^2 +...
    (MP_Posi(:,1) - ref_posi(2)).^2).^0.5;

a1 = -0.1379;
b1 = 562.9572;
a2 = -0.1562;
b2 = 657.3926;
slct_ind = (MP_Posi(:,1) > a1*MP_Posi(:,2)+b1) &...
            (MP_Posi(:,1) < a2*MP_Posi(:,2)+b2); 
        
remain_ind = 1:locator_num; remain_ind = remain_ind(slct_ind);
remain_num = length(remain_ind);
        
[~,reorder_ind] = sort(dist(slct_ind),'ascend');

for i = 1:sinNum
fig_h = figure('Position',[60 120 1800 860],'Color','w');
colormap(jet(1000));
subplot('Position',[0.54 0.05 0.42 0.9]);
imshow(imread([MapPath,'DigitIISinusoid_HandContour.jpg']));
hold on;
scatter(MP_Posi(slct_ind,2),MP_Posi(slct_ind,1),100,PC_resp(i,slct_ind),...
    'filled');
hold off;
camroll(90); 
c_h = colorbar('Location','north','Position',[0.65,0.8,0.22,0.03]);
c_h.Label.String = 'Average PC Firing Rate (Hz)';
title(['Sine',tavgData{i,2}],'Position', [100, -200]);   
set(gca,'FontSize',16);

subplot('Position',[0.05 0.25 0.48 0.5]);
hold on;
for j = 1:remain_num
    spike_num = length(PC_spikes{i,remain_ind(reorder_ind(j))});
    plot(PC_spikes{i,remain_ind(reorder_ind(j))},j*ones(1,spike_num),'.k');
end
hold off;
xlabel('Time (secs)');
ylabel('PC#');

set(gca,'FontSize',16, 'YDir','reverse','XAxisLocation','top');

% set(fig_h,'PaperOrientation','landscape');
print(fig_h,sprintf('RasterPlots/Digit2PCResp_Sin%s',tavgData{i,2}),...
    '-djpeg');
% print(fig_h,sprintf('RasterPlots/Digit2PCResp_Sin%s',tavgData{i,2}),...
%     '-dpdf','-fillpage','-painters');
pause(0.01); close(fig_h);
end

%% -----------------Average Firing Rate Whole Hand Pattern ----------------
if 1 %---------------------------------------------------------------------

cRange = [min(PC_resp(:)), max(PC_resp(:))];

figure('Position',[60,60,1840,860],'Color','w');
colormap(jet(1000));
rowN = 3;
colN = 6;

for i = 1:sinNum
    r_i = ceil(i/6);
    c_i = mod(i-1,6);
    subplot('Position',[(c_i/colN),1-(r_i/rowN),1/colN,1/rowN]);

    imshow(imread([MapPath,'DigitIISinusoid_HandContour.jpg']));

    hold on;
%     scatter(MP_Posi(:,2),MP_Posi(:,1),20,tavgData{i,1},'filled');
    scatter(MP_Posi(:,2),MP_Posi(:,1),25,PC_resp(i,:),'filled');
    hold off;
    
    xlabel(tavgData{i,2});
%     c_h = colorbar;
    caxis(cRange);
    set(gca,'FontSize',12);
end

subplot('Position',[(5/colN),1-(3/rowN),1/colN,1/rowN]);
% imshow(imread([MapPath,'DigitIISinusoid_MP.jpg']));
axis off;
c_h = colorbar('Location','north');
c_h.Label.String = 'PC Firing Rate (Hz)';
caxis(cRange);

set(gca,'FontSize',16);
end %----------------------------------------------------------------------