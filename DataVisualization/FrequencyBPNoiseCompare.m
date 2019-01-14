%% Visualize PDV data 
% Created on 01/11/2019
% -------------------------------------------------------------------------
clear all
% -------------------------------------------------------------------------
MapPath = './DigitIIBPNoise/';

Alpha = 25.5; % (mm)
C = 0.4;
px2mm = 0.289;  % (mm/pixel)
interp_radius = 200;
maskThreshold = 200;

Fs = 20000;
% -------------------------------------------------------------------------
%% Plot 1
DataPath = '../Data_DigitIIBPNoise/SubjectBha_Session02.svd';

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

maskImg = imread([MapPath,'DigitIIBPNoise_Mask.jpg']);
maskImg = rgb2gray(maskImg);
pointImg = imread([MapPath,'DigitIIBPNoise_MP.jpg']);
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

% -------------------------------------------------------------------------
tavgData = cell(3,2);
figure('Position',[60,60,1840,880],'Color','w')
slct_ind = 14001:20800;
[ seg_spect, f ] = spectr(ref(1,slct_ind), Fs, [0 800]);
plot(f,seg_spect); xlabel('Frequency (Hz)'); hold on;
tavgData{1,2} = sprintf('%.0fHz', f*(seg_spect./sum(seg_spect)));
tavgData{1,1} = rms(y(:,slct_ind),2)';

slct_ind = 29281:32340;
[ seg_spect, f ] = spectr(ref(1,slct_ind), Fs, [0 800]);
plot(f,seg_spect);
tavgData{2,2} = sprintf('%.0fHz', f*(seg_spect./sum(seg_spect)));
tavgData{2,1} = rms(y(:,slct_ind),2)';

slct_ind = 68001:71000;
[ seg_spect, f ] = spectr(ref(1,slct_ind), Fs, [0 800]);
plot(f,seg_spect);
tavgData{3,2} = sprintf('%.0fHz', f*(seg_spect./sum(seg_spect)));
tavgData{3,1} = rms(y(:,slct_ind),2)';

legend(tavgData{1,2},tavgData{2,2},tavgData{3,2});

% -------------------------------------------------------------------------
figure('Position',[60,60,1840,880],'Color','w')
for i = 1:3
    subplot(2,3,i)
    colormap(jet(1000));
    imshow(imread([MapPath,'DigitIIBPNoise_MP.jpg']))
    % for i = 1:locator_num
    %     text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','r')
    % end
    hold on
    scatter(MP_Posi(:,2),MP_Posi(:,1),20,tavgData{i,1},'filled');
    hold off
    title(['PDV Measurement: ',tavgData{i,2}])

    interpImg = interpMP(maskImg, MP_Posi, tavgData{i,1},...
        maskThreshold, interp_radius, Alpha, C, px2mm);

    subplot(2,3,3+i)
    surf(flipud(interpImg),'EdgeColor','none');
    view(2)
    axis equal; axis off;
    title(['Interpolated Measurement: ',tavgData{i,2}])
end
