%% Visualize PDV data 2
% Created on 10/05/2018
% -------------------------------------------------------------------------
clear all
% -------------------------------------------------------------------------
MapPath = './TapDigitII_SBJ2/';

Alpha = 25.5; % (mm)
C = 0.4;
px2mm = 0.289;  % (mm/pixel)
interp_radius = 200;
maskThreshold = 200;
% -------------------------------------------------------------------------
%% Dorsal
DataPath = 'TapDigitII_Dorsal_SBJ2.svd';

if ~exist('data_info','var')
    [t,y,data_info] = GetPointData(DataPath, 'Time', 'Vib', 'Velocity',...
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

tavgData = rms(y,2)';

maskImg = imread([MapPath,'Dorsal_Mask.jpg']);
maskImg = rgb2gray(maskImg);
pointImg = imread([MapPath,'Dorsal_MP.jpg']);
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

figure('Position',[60,20,820,940],'Color','w')

subplot(2,2,1)
colormap(jet(1000));
imshow(imread([MapPath,'Dorsal_MP.jpg']))
% for i = 1:locator_num
%     text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','r')
% end
hold on
scatter(MP_Posi(:,2),MP_Posi(:,1),20,tavgData,'filled');
hold off
title('PDV Measurement (Dorsal)')

interpImg = interpMP(maskImg, MP_Posi, tavgData,...
    maskThreshold, interp_radius, Alpha, C, px2mm);

subplot(2,2,3)
surf(flipud(interpImg),'EdgeColor','none');
view(2)
axis equal; axis off;
title('Interpolated Measurement (Dorsal)')

%% Volar
clearvars -except MapPath Alpha C px2mm interp_radius maskThreshold

DataPath = 'TapDigitII_Volar_SBJ2.svd';

if ~exist('data_info','var')
    [t,y,data_info] = GetPointData(DataPath, 'Time', 'Vib', 'Velocity',...
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

tavgData = rms(y,2)';

maskImg = imread([MapPath,'Volar_Mask.jpg']);
maskImg = rgb2gray(maskImg);
pointImg = imread([MapPath,'Volar_MP.jpg']);
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

subplot(2,2,2)
colormap(jet(1000));
imshow(imread([MapPath,'Volar_MP.jpg']))
% for i = 1:locator_num
%     text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','r')
% end
hold on
scatter(MP_Posi(:,2),MP_Posi(:,1),20,tavgData,'filled');
hold off
title('PDV Measurement (Volar)')

interpImg = interpMP(maskImg, MP_Posi, tavgData,...
    maskThreshold, interp_radius, Alpha, C, px2mm);

subplot(2,2,4)
surf(flipud(interpImg),'EdgeColor','none');
view(2)
axis equal; axis off;
title('Interpolated Measurement (Volar)')