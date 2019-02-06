%% Estimate wave speed on volar hand
% Created on 02/05/2019 based on 'WaveSpeed.m'
% -------------------------------------------------------------------------
Data_Path = '../Data_DigitIIBPNoise/';
pathInfo = dir([Data_Path,'*Session02.svd']);
sbj_num = length(pathInfo);

% -------------------------------------------------------------------------
BPNoiseFreq = [20,40,80,120,160,200,240,280,320,360,400,440,480,520,560,...
    600,640]; 
BPNoiseNum = length(BPNoiseFreq); % Number of BP noise in each measurement

% Data preprocessing ------------------------------------------------------
if ~(exist('DataInfo','var') && exist('BPNoiseY','var'))
DataInfo = cell(sbj_num,3); % Label (ID) of subject, info and MP coordinates
BPNoiseY = cell(sbj_num,BPNoiseNum); % BP Noise segmented data
for sbj_i = 1:sbj_num
    a_name = pathInfo(sbj_i).name;
    DataInfo{sbj_i,1} = a_name(8:10);
    
    dataName = [Data_Path,a_name];
    
    [t,y,DataInfo{sbj_i,2}] = GetPointData(dataName, 'Time', 'Vib', 'Velocity',...
        'Samples', 0, 0);
    [t_ref,ref,~] = GetPointData(dataName, 'Time', 'Ref1', 'Voltage',...
        'Samples', 0, 0);
    
    dt = mean(diff(t));
    Fs = 1/dt;
    fprintf('Sampling Frequency = %.1f Hz\n',Fs);
    
    % Getting Measurement Point (MP) coordinates
    DataInfo{sbj_i,3} = GetXYZCoordinates(dataName, 0);
    
    locator_num = size(y,1);
    
%     y_filt = bandpass(y',[10 800],Fs)';
    y_filt = lowpass(y',8000,Fs)';

    % Segment the measurement using the reference signal
    zero_ind = (abs(ref(end,:)) < 9e-4);
    onsite_ind = movmedian(double(zero_ind),1200);
    onsite_ind = [0,diff(onsite_ind)];
    onsite_ind = (onsite_ind < -0.1);

    seg_ind = find(onsite_ind == 1);
    seg_ind = seg_ind(1:2:end);
    seg_ind = [seg_ind, (2*seg_ind(end)-seg_ind(end-1))];
    
    % Display data segmentation result
    figure('Position',[60,60,1840,880],'Color','w');
    plot(onsite_ind);
    title(sprintf('Subject: %s - Segment Number = %d',DataInfo{sbj_i,1},...
        length(seg_ind)-1));
    hold on
    plot(ref(end,:)./max(ref(end,:)));
    plot(y(end,:)./max(abs(y(end,:))));
    plot(y_filt(end,:)./max(abs(y(end,:))));
    plot(seg_ind,ones(size(seg_ind)),'*r')
    hold off
    drawnow;
    
    if (length(seg_ind)-1) == BPNoiseNum
        for BPN_i = 1:BPNoiseNum     
%             wb_h = waitbar(0, 'O', 'Name',sprintf('SBJ-%s %dHz',...
%             DataInfo{sbj_i,1}, BPNoiseFreq(BPN_i)));
%             for MP_i = 1:locator_num
%                 waitbar(MP_i/locator_num,wb_h,...
%                     sprintf('MP %d/%d',MP_i,locator_num));
%                 [r,~] = xcorr(ref(end,:),ref(MP_i,:));
%                 [~,max_i] = max(r);
%                 delay_i = max_i-length(ref(end,:));
%                 if (delay_i ~= 0)
%                     fprintf('Subject: %s - %dHz - MP%d Sync Error\n',...
%                         DataInfo{sbj_i,1},BPNoiseFreq(BPN_i),MP_i);
%                 end
%                 BPNoiseY{sbj_i,BPN_i}(MP_i,:) =...
%                     y_filt(MP_i,(seg_ind(BPN_i):seg_ind(BPN_i+1))+delay_i);
%             end
%             close(wb_h);
            BPNoiseY{sbj_i,BPN_i}=y_filt(:,seg_ind(BPN_i):seg_ind(BPN_i+1));
        end
    else
        error('Data segmentation failed!')
    end
end
close all;
end

%% Estimate wave speed
MP_radius = 0.02; % (Unit: m) estimate speed from nearest neighboor within 1 cm
Estim_radius = 0.03; % (Unit: m) area of estimation

estimate_num = cell(sbj_num,BPNoiseNum);
wavespeed = cell(sbj_num,BPNoiseNum);
medianSpeed = cell(sbj_num,BPNoiseNum);
avgSpeed = zeros(sbj_num,BPNoiseNum);
for BPN_i = 1:BPNoiseNum
    for sbj_i = 1:sbj_num
        y = BPNoiseY{sbj_i,BPN_i};
        XYZ = DataInfo{sbj_i,3};
        MP_num = size(y,1);
        
        estimate_num{sbj_i,BPN_i} = zeros(MP_num,1);
        wavespeed{sbj_i,BPN_i} = NaN(MP_num,100);
        
        wb_h = waitbar(0, 'O', 'Name',sprintf('SBJ-%s %dHz',...
            DataInfo{sbj_i,1}, BPNoiseFreq(BPN_i)));
        
        tip_ind = find((XYZ(:,1) < (XYZ(end,1)+Estim_radius)) &...
                (XYZ(:,1) > (XYZ(end,1)-Estim_radius)) &...
                (XYZ(:,2) < (XYZ(end,2)+Estim_radius)) &...
                (XYZ(:,2) > (XYZ(end,2)-Estim_radius)) &...
                (XYZ(:,3) < (XYZ(end,3)+Estim_radius)) &...
                (XYZ(:,3) > (XYZ(end,3)-Estim_radius)))';           
            
        for i = tip_ind
            waitbar(i/MP_num,wb_h,sprintf('MP %d/%d',i,MP_num));
            ind = find((XYZ(:,1) < (XYZ(i,1)+MP_radius)) &...
                (XYZ(:,1) > (XYZ(i,1)-MP_radius)) &...
                (XYZ(:,2) < (XYZ(i,2)+MP_radius)) &...
                (XYZ(:,2) > (XYZ(i,2)-MP_radius)) &...
                (XYZ(:,3) < (XYZ(i,3)+MP_radius)) &...
                (XYZ(:,3) > (XYZ(i,3)-MP_radius)));
            ind = ind(ind~=i); % Indcies of neighboor MP
            
            if ~isempty(ind)
                for j = 1:length(ind)
                    [r,lags] = xcorr(y(i,:),y(ind(j),:));
                    [~,max_i] = max(r);
                    delay_t = lags(max_i)*dt; % (sec)
                    dist_3D = sum((XYZ(i,:) - XYZ(ind(j),:)).^2).^0.5; % (m)
                    
                    if delay_t ~= 0
                        estimate_num{sbj_i,BPN_i}(i) =...
                            estimate_num{sbj_i,BPN_i}(i) +1;
                        wavespeed{sbj_i,BPN_i}(i,...
                            estimate_num{sbj_i,BPN_i}(i))=...
                            abs(dist_3D/delay_t); % (m/sec)
                    end
                end
            end
        end
        close(wb_h);
        
        medianSpeed{sbj_i,BPN_i} = median(wavespeed{sbj_i,BPN_i},2,'omitnan');
        
        avgSpeed(sbj_i,BPN_i) = mean(medianSpeed{sbj_i,BPN_i},'omitnan');       
        
        fprintf('SBJ-%s %dHz Avg. Speed = %.0f m/s \n',DataInfo{sbj_i,1},...
            BPNoiseFreq(BPN_i), avgSpeed(sbj_i,BPN_i));
    end
end

%% Plot the Average Speed
if 1 %---------------------------------------------------------------switch
    
slct_sbj = [1,2,3,4,5];
slct_ind = 2:BPNoiseNum;
figure('Position',[160,80,1620,600],'Color','w')
plot(BPNoiseFreq(slct_ind),avgSpeed(slct_sbj,slct_ind)','.-','MarkerSize',20); 
% boxplot(avgSpeed(:,slct_ind))
box off;
xlim([0 690])
xlabel('Frequency (Hz)')
ylabel('Avg. Speed (m/s)')
yRange = ylim();

% Curve fitting
fInd = repmat(BPNoiseFreq(slct_ind),[1,length(slct_sbj)])';
slctSpeed = avgSpeed(slct_sbj,slct_ind)'; slctSpeed = slctSpeed(:);

fitOrderNum = 50;
fit_h = cell(fitOrderNum,1);
fit_err = NaN(fitOrderNum,1);
for i = 1:fitOrderNum
[fit_h{i}, errInfo] = polyfit(fInd, slctSpeed, i);
fit_err(i) = errInfo.normr;
end
[lowestErr,lowErrInd] = min(fit_err);
fprintf('Order of fitting with lowest error (%.1f) = %d\n',lowestErr,...
    lowErrInd);

slct_fit = 6;
f_intp = BPNoiseFreq(slct_ind(1)):BPNoiseFreq(slct_ind(end));
fit_speed = polyval(fit_h{slct_fit},f_intp);

hold on
% plot(fInd,slctSpeed,'ok')
plot(f_intp,fit_speed,'r','LineWidth',4)
hold off
ylim(yRange);

legend([DataInfo(slct_sbj,1);'Fitting'],'box','off');
set(gca,'FontSize',16);

end %------------------------------------------------------------switch end

%% Plot the speed distribution on the whole hand
Alpha = 25.5; % (mm)
C = 0.4;
px2mm = 0.289;  % (mm/pixel)
interp_radius = 100;
maskThreshold = 100;

sbj_path = 'DigitIIBPNoise_Bha/';
sbj_i = 4; % Subject: Bha

maskImg = imread([sbj_path,'DigitIIBPNoise_Mask.jpg']);
maskImg = rgb2gray(maskImg);
pointImg = imread([sbj_path,'DigitIIBPNoise_MP.jpg']);
pointImgGray = rgb2gray(pointImg); 

locator_num = length(medianSpeed{sbj_i,1});

% Low pass filter the chosen grid and computable area
lpF = ones(3)/9; % 3-by-3 mean filter
ptImg_PFilter = uint8(filter2(lpF,pointImgGray));
MP_Posi = findMP(ptImg_PFilter, locator_num, 128);
MP_Posi = round(MP_Posi); % Round the position to get index
xDiff = diff(MP_Posi(:,2));
xDiff(abs(xDiff)<5) = 0;
MP_Score = cumsum([0; xDiff]).*10000 + MP_Posi(:,1);
[~,ind] = sort(MP_Score);
MP_Posi = MP_Posi(ind,:);

handRangeX = [min(MP_Posi(:,2))-20,max(MP_Posi(:,2))+20];
handRangeY = [min(MP_Posi(:,1))-20,max(MP_Posi(:,1))+20];

% % % Validate MP position estimation
% imshow(pointImg);
% for i = 1:locator_num
%     text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','r')
% end

sbp_w = 0.24; % Subplot width (0 - 1)
sbp_h = 0.24; % Subplot height (0 - 1)

slct_ind = 2:BPNoiseNum;
figure('Position',[160,80,1620,860],'Color','w');
colormap(jet(1000))
for i = 1:length(slct_ind)
    r_n = floor((i-1)/4)+1;
    c_n = mod(i-1,4);
    subplot('Position',[0.005+c_n*sbp_w,0.98-r_n*sbp_h,sbp_w,sbp_h-0.01]);
    sbj_avgSpeed = medianSpeed{sbj_i,slct_ind(i)};
    
    discard_ind = (sbj_avgSpeed > 80);
    sbj_avgSpeed(discard_ind) = NaN;
    
    imshow([sbj_path,'DigitIIBPNoise_Contour.jpg']);
    xlim(handRangeX); ylim(handRangeY)
    hold on
    scatter(MP_Posi(:,2),MP_Posi(:,1),30,sbj_avgSpeed,'filled');
    hold off
%     caxis(yRange)
%     caxis([0 1])
    c_h = colorbar('Ticks',[]);
    c_h.Label.String = sprintf('%.0f - %.0f (m/s)',caxis);

%     
%     interpImg = interpMP(maskImg, MP_Posi,sbj_avgSpeed',maskThreshold,...
%         interp_radius, Alpha, C, px2mm);      
%     sc_h = imagesc(interpImg); 
%     set(sc_h,'AlphaData',~isnan(interpImg));
%     view(2);
%     xlim(handRangeX); ylim(handRangeY)
%     axis equal; axis off;
%     title(sprintf('%d Hz',BPNoiseFreq(slct_ind(i))));
%     c_h = colorbar('Ticks',[]);
%     c_h.Label.String = sprintf('%.0f - %.0f (m/s)',caxis);

    title(sprintf('%d Hz %.0f - %.0f (m/s)',BPNoiseFreq(slct_ind(i)),caxis));
    drawnow;
end
