%% Visualize PDV data 
% Created on 02/05/2019 based on 'FrequencySinCompare.m'
% -------------------------------------------------------------------------
clear all
% -------------------------------------------------------------------------
MapPath = './DigitIIChirp/';

Alpha = 25.5; % (mm)
C = 0.4;
px2mm = 0.289;  % (mm/pixel)
interp_radius = 200;
maskThreshold = 200;

Fs = 20000;
% -------------------------------------------------------------------------
%% Plot 1
DataPath = '../Data_DigitIIChirp/SubjectBha_Session01.svd';

if ~exist('data_info','var')
    [f,y,data_info] = GetPointData(DataPath, 'FFT', 'Vib', 'Velocity',...
        'Magnitude', 0, 0);
    
    [f_ref,ref,~] = GetPointData(DataPath, 'FFT', 'Ref1', 'Voltage',...
        'Magnitude', 0, 0);
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

maskImg = imread([MapPath,'DigitIIChirp_Mask.jpg']);
maskImg = rgb2gray(maskImg);
pointImg = imread([MapPath,'DigitIIChirp_MP.jpg']);
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

handRangeX = [min(MP_Posi(:,2))-24,max(MP_Posi(:,2))+24];
handRangeY = [min(MP_Posi(:,1)),max(MP_Posi(:,1))];

%% Frequency band analysis
% freqIntrst=logspace(log10(40),log10(640),13);
freqIntrst = [40, 80, 160, 320, 480, 640];
numFreqs=length(freqIntrst);
freqVector=[];
for i=1:numFreqs
    freqVector=[freqVector; [freqIntrst(i)*.9,freqIntrst(i)*1.1]];
end

%% ------------------------------------------------------------------------
% Validate the data input
temp_fig = figure('Position',[60,60,1840,880],'Color','w');
subplot('Position',[0.535 0.005 0.46 0.9]);
imshow(imread([MapPath,'DigitIIChirp_MP.jpg']));
for i = 1:locator_num
    text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','r');
end
drawnow;
subplot('Position',[0.025 0.2 0.5 0.6]);
yyaxis right;
plot(f_ref,ref(end,:),'b');
ylabel('Reference Magnitude (m/s)');
xlabel('Frequency (Hz)')
yyaxis left;
plot(f,y(end,:),'g');
ylabel('Measurement Magnitude (m/s)');
hold on;
seg_En = NaN(numFreqs,locator_num);
for i = 1:numFreqs
    seg_ind = (f > freqVector(i,1)) & (f < freqVector(i,2));
    seg_spectr = y(:,seg_ind);
    seg_En(i,:) = rms(seg_spectr,2)';
    plot(f(seg_ind),seg_spectr(end,:)','-r')
    drawnow;
end
hold off;

close(temp_fig);

%% ------------------------------------------------------------------------
fig_h = figure('Position',[60,60,1840,880],'Color','w');
colormap(jet(1000));
subplot_width = 0.16;
for i = 1:numFreqs
    subplot('Position',...
        [0.005+(i-1)*subplot_width,0.2,subplot_width-0.02,0.6]); 
    
    interpImg = interpMP(maskImg, MP_Posi, seg_En(i,:),...
        maskThreshold, interp_radius, Alpha, C, px2mm);
    
    sc_h = imagesc(interpImg); 
    set(sc_h,'AlphaData',~isnan(interpImg));
    view(2); xlim([handRangeX(1)+500, handRangeX(2)]); ylim(handRangeY);
    axis equal; axis off;
    
    fprintf('%d Hz [%.1f - %.1f mm/s]\n',freqIntrst(i),1000*caxis);
    
    set(gca, 'units', 'normalized'); 
    title(sprintf('%d Hz',freqIntrst(i)))
    
    if i == round(0.5*numFreqs)
        c_h = colorbar('Location','south','Ticks',[]);
        c_h.Label.String = sprintf('Velocity %.2f - %.2f (mm/s)\n',...
            1000*caxis);
        c_h.Label.Position = [0 -2 0];
        disp(c_h.Label.String{1});
    end
    
    drawnow;
end
print(fig_h,'ChirpSpectrEn','-dpdf','-bestfit','-r600','-painters')