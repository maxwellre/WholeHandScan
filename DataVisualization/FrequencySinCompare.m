%% Visualize PDV data 
% Created on 01/14/2019 based on 'FrequencyBPNoiseCompare.m'
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
sinNum = 17;
tavgData = cell(sinNum,2);

zero_ind = (abs(ref(1,:)) < 5e-4);
onsite_ind = movmedian(double(zero_ind),19);

onsite_ind = [0,diff(onsite_ind)];
onsite_ind = (onsite_ind == -1);

figure('Position',[60,60,1840,880],'Color','w')
subplot(2,1,1)
plot(onsite_ind);
hold on
plot(ref(1,:)./max(ref(1,:)));
plot(y_filt(end,:)./max(abs(y_filt(end,:))));

seg_ind = find(onsite_ind == 1);
seg_ind = [seg_ind, (2*seg_ind(end)-seg_ind(end-1))];
hold off

subplot(2,1,2)
hold on
for i = 1:sinNum
    slct_ind = seg_ind(i):seg_ind(i+1);
    [ seg_spect, f ] = spectr(ref(1,slct_ind), Fs, [0 800]);
%     [ seg_spect, f ] = spectr(y(end,slct_ind), Fs, [0 800]);
    plot(f,seg_spect); xlabel('Frequency (Hz)'); 
    tavgData{i,2} = sprintf('%.0fHz', f*(seg_spect./sum(seg_spect)));
    tavgData{i,1} = rms(y_filt(:,slct_ind),2)';
end
legend(tavgData(:,2));
hold off

%% ------------------------------------------------------------------------
fig_h = figure('Position',[60,60,1640,580],'Color','w');
colormap(jet(1000));
slct_sin = [1,2,3,5,9,13,17];
for i = 1:length(slct_sin)
    interpImg = interpMP(maskImg, MP_Posi, tavgData{i,1},...
        maskThreshold, interp_radius, Alpha, C, px2mm);

%     subplot(1,length(slct_sin),i)
    surf(flipud(interpImg),'EdgeColor','none');
    view(2)
    axis equal; axis off;
    print(fig_h,sprintf('./FreqCompareSine/FreqCompareSine%s',...
        tavgData{i,2}),'-dpng')
    pause(0.1);
end


% figure('Position',[60,60,1840,880],'Color','w')
% for i = 1:sinNum
%     subplot(5,7,i)
%     colormap(jet(1000));
%     imshow(imread([MapPath,'DigitIISinusoid_MP.jpg']))
% %     for i = 1:locator_num
% %         text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','r')
% %     end
%     hold on
%     scatter(MP_Posi(:,2),MP_Posi(:,1),5,tavgData{i,1},'filled');
%     hold off
%     if i == 1
%         ylabel('PDV Measurement')
%     end
%     title(tavgData{i,2})
% 
%     interpImg = interpMP(maskImg, MP_Posi, tavgData{i,1},...
%         maskThreshold, interp_radius, Alpha, C, px2mm);
% 
%     subplot(5,7,sinNum+i)
%     surf(flipud(interpImg),'EdgeColor','none');
%     view(2)
%     axis equal; axis off;
%     if i == 1
%         ylabel('Interpolated Measurement')
%     end
% end
% set(gcf,'PaperOrientation','landscape');
% print(fig_h,'FreqSinCompare','-dpdf',...
%     '-fillpage','-painters')