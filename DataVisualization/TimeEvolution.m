%% Visualize PDV data 
% Created on 01/17/2019 based on 'FrequencySinCompare.m'
% -------------------------------------------------------------------------
clear all
% -------------------------------------------------------------------------
MapPath = './DigitIISinusoid/';

Alpha = 25.5; % (mm)
C = 0.4;
px2mm = 0.289;  % (mm/pixel)
interp_radius = 200;
maskThreshold = 200;

Fs = 20000;
% -------------------------------------------------------------------------
%% Plot 1
DataPath = '../Data_DigitIISinusoid/SubjectBha_Session04.svd';

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

% y_filt = movmedian(y',20)';
y_filt = y;

% -------------------------------------------------------------------------
% sinNum = 17;
% tavgData = cell(sinNum,2);
% 
% zero_ind = (abs(ref(1,:)) < 5e-4);
% onsite_ind = movmedian(double(zero_ind),19);
% 
% onsite_ind = [0,diff(onsite_ind)];
% onsite_ind = (onsite_ind == -1);

% figure('Position',[60,60,1840,880],'Color','w')
% subplot(2,1,1)
% plot(onsite_ind);
% hold on
% plot(ref(1,:)./max(ref(1,:)));
% plot(y_filt(end,:)./max(abs(y_filt(end,:))));
% 
% seg_ind = find(onsite_ind == 1);
% seg_ind = [seg_ind, (2*seg_ind(end)-seg_ind(end-1))];
% hold off
% 
% subplot(2,1,2)
% hold on
% for i = 1:sinNum
%     slct_ind = seg_ind(i):seg_ind(i+1);
%     [ seg_spect, f ] = spectr(ref(1,slct_ind), Fs, [0 800]);
% %     [ seg_spect, f ] = spectr(y(end,slct_ind), Fs, [0 800]);
%     plot(f,seg_spect); xlabel('Frequency (Hz)'); 
%     tavgData{i,2} = sprintf('%.0fHz', f*(seg_spect./sum(seg_spect)));
%     tavgData{i,1} = rms(y_filt(:,slct_ind),2)';
% end
% legend(tavgData(:,2));
% hold off

%% ------------------------------------------------------------------------
t = (0:length(y_filt)-1)/Fs;
t_display = [1.1477 1.16005];
ind_display = find(((t > t_display(1)) & (t < t_display(2)))==1);

% % ind_display = [4.1995e4, 4.2025e4];
% ind_display = [4.706e4, 4.7082e4];
t_display = ind_display/Fs;

frame_num = 7;

frame_ind = round(linspace(ind_display(1),ind_display(end),frame_num));

figure('Position',[60,260,1640,480],'Color','w')
plot(t,y_filt(end,:));
hold on;
for k = 1:frame_num
    plot([frame_ind(k),frame_ind(k)]/Fs,[0,0.02],'r');
end
hold off;
% xlim([t_display(1)-0.002 t_display(2)+0.002])

%%
y_display = y_filt(:,frame_ind)';
y_range = [min(y_display(:)),max(y_display(:))];

fig_h = figure('Position',[60,60,1840,780],'Color','w');
colormap(jet(1000));
for k = 1:frame_num
    interpImg = interpMP(maskImg, MP_Posi, y_display(k,:),...
        maskThreshold, interp_radius, Alpha, C, px2mm);
    subplot(2,4,k);
    surf(flipud(interpImg),'EdgeColor','none');
    caxis(y_range);
    view(2);
    axis equal; axis off;   
    drawnow;
end
subplot(2,4,frame_num+1)
caxis(y_range);
colorbar

%% print(gcf,'NewFigure','-dpdf','-bestfit','-r600','-painters')