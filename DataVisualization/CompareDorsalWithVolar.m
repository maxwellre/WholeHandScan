%% Visualize PDV data
% Created on 09/29/2018
% -------------------------------------------------------------------------
MapPath = './VisualMap/';

Alpha = 30;
Range = 200;
interp_radius = 300;
maskThreshold = 200;

%% Dorsal
maskImg_Dorsal = imread([MapPath,'Dorsal_Mask.jpg']);
maskImg_Dorsal = rgb2gray(maskImg_Dorsal);
pointImg_Dorsal = imread([MapPath,'Dorsal_MP.jpg']);

% Compute the optical flow vector locations based on dot maps--------------
locator_num = 348; % Stream line locator number (74)
pointImg_Dorsal = pointImg_Dorsal(:,:,2); % Green points

% Low pass filter the chosen grid and computable area
lpF=ones(3)/9; % 3-by-3 mean filter
ptImg_PFilter=uint8(filter2(lpF,pointImg_Dorsal));
MP_Posi = findMP(ptImg_PFilter, locator_num);
MP_Posi = round(MP_Posi); % Round the position to get index
xDiff = diff(MP_Posi(:,2));
xDiff(abs(xDiff)<5) = 0;
MP_Score = cumsum([0; xDiff]).*10000 + MP_Posi(:,1);
[~,ind] = sort(MP_Score);
MP_Posi = MP_Posi(ind,:);

load('TapDigitII_Dorsal.mat');
tavgData = rms(PDV_data.Data,1);
MP_ID = PDV_data.Posi(1,:);
[~, ID_ind] = sort(MP_ID);
tavgData = tavgData(ID_ind);

pointImg_Dorsal = imread([MapPath,'Dorsal_MP.jpg']);
figure('Position',[60,80,1800,640])
subplot(2,2,1)
colormap(jet(1000));
imshow(pointImg_Dorsal)
% for i = 1:locator_num
%     text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','r')
% end
tavgData = tavgData./max(tavgData);
hold on
scatter(MP_Posi(:,2),MP_Posi(:,1),30,tavgData,'filled');
hold off
title('PDV Measurement')

interpImg = interpMP_mex(maskImg_Dorsal, MP_Posi, tavgData,...
    maskThreshold, interp_radius, Alpha, Range);

subplot(2,2,3)
surf(flipud(interpImg),'EdgeColor','none');
view(2)
axis equal; axis off;
title('Interpolated Measurement')

%% Volar 
maskImg_Volar = imread([MapPath,'Volar_Mask.jpg']);
maskImg_Volar = rgb2gray(maskImg_Volar);
pointImg_Volar = imread([MapPath,'Volar_MP.jpg']);

% Compute the optical flow vector locations based on dot maps--------------
locator_num = 348; % Stream line locator number (74)
pointImg_Volar = pointImg_Volar(:,:,2); % Green points

% Low pass filter the chosen grid and computable area
lpF=ones(3)/9; % 3-by-3 mean filter
ptImg_PFilter=uint8(filter2(lpF,pointImg_Volar));
MP_Posi = findMP(ptImg_PFilter, locator_num, 128);
MP_Posi = round(MP_Posi); % Round the position to get index
xDiff = diff(MP_Posi(:,2));
xDiff(abs(xDiff)<5) = 0;
MP_Score = cumsum([0; xDiff]).*10000 + MP_Posi(:,1);
[~,ind] = sort(MP_Score);
MP_Posi = MP_Posi(ind,:);

load('TapDigitII_Volar.mat');
tavgData = rms(PDV_data.Data,1);
MP_ID = PDV_data.Posi(1,:);
[~, ID_ind] = sort(MP_ID);
tavgData = tavgData(ID_ind);

pointImg_Volar = imread([MapPath,'Volar_MP.jpg']);
subplot(2,2,2)
colormap(jet(1000));
imshow(pointImg_Volar)
% for i = 1:locator_num
%     text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','r')
% end
tavgData = tavgData./max(tavgData);
hold on
scatter(MP_Posi(:,2),MP_Posi(:,1),30,tavgData,'filled');
hold off
title('PDV Measurement')

interpImg = interpMP_mex(maskImg_Volar, MP_Posi, tavgData,...
    maskThreshold, interp_radius, Alpha, Range);

% figure('Position',[20,50,1800,840])
subplot(2,2,4)
surf(flipud(interpImg),'EdgeColor','none');
view(2)
axis equal; axis off;
title('Interpolated Measurement')

