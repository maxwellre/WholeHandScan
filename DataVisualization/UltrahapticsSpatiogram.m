%% Visualize PDV data 
% Created on 02/02/2019 based on 'UltrahapticsAllFreqBandVisualization.m'
% -------------------------------------------------------------------------
% clear all
% -------------------------------------------------------------------------

dataName = {'Greg_MovingSpot_1ms_Dir1','Greg_MovingSpot_1ms_Dir2',...
    'Greg_MovingSpot_2ms_Dir1','Greg_MovingSpot_2ms_Dir2',...
    'Greg_MovingSpot_4ms_Dir1','Greg_MovingSpot_4ms_Dir2',...
    'Greg_MovingSpot_7ms_Dir1','Greg_MovingSpot_7ms_Dir2',...
    'Greg_MovingSpot_11ms_Dir1','Greg_MovingSpot_11ms_Dir2',...
    'Greg_MovingSpot_15ms_Dir1','Greg_MovingSpot_15ms_Dir2'};
dataNum = length(dataName);
TrialNum = 2;

for d_i = 1:dataNum
    for t_i = 1:TrialNum
    
% -------------------------------------------------------------------------
MapPath = './UltrahapticsMeasurement/';
DataPath = sprintf('../Data_Ultrahaptics/%s_1.svd',dataName{d_i});
OutputPath = './UltrahapticsVideo/';

cleanDataPath = sprintf('../Data_Ultrahaptics/%s_%d.mat',dataName{d_i},t_i);
load(cleanDataPath);

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


%% Display selected measurement points
slct_ind = [297,295,290,284,277,269,259,249,237,227,215,201,185,170,154,...
    138,124,111,99,88,77,66,56,46,36,27,18];

% temp_fig = figure('Position',[60,60,1840,880],'Color','w');
% imshow(imread([MapPath,'Greg_MovingSpot_1ms_Dir1_1_MP.jpg']))
% hold on;
% for i = 1:locator_num
%     text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','c')
% end
% scatter(MP_Posi(slct_ind,2),MP_Posi(slct_ind,1),18,'r','filled');
% hold off;
% close(temp_fig);

MP_dist = 1000*((XYZ(slct_ind,1)-XYZ(slct_ind(1),1)).^2 +...
    (XYZ(slct_ind,2)-XYZ(slct_ind(1),2)).^2 +...
    (XYZ(slct_ind,3)-XYZ(slct_ind(1),3)).^2).^0.5; % m to mm

y_slct = y_vib_sync(slct_ind,:)';

[ spatiogram, freq ] = spectr( y_slct, Fs, [40 1000] );

%% Plot 
curr_fig = figure('Position',[60,360,1840,240],'Color','w');
colormap(jet(1000));
imagesc(20*log10(spatiogram'))
caxis([-160 -60]);
x_h = xticks();
xticklabels(round(freq(x_h))); xlabel('Frequency (Hz)');
y_h = yticks();
yticklabels(round(MP_dist(y_h))); ylabel('Distance (mm)');
title([dataName{d_i},sprintf(' Trial%d',t_i)],'Interpreter', 'none');
cb_h = colorbar;
cb_h.Label.String = 'DB';
drawnow;
    end
end