dataName = {'Greg_MovingSpot_1ms_Dir1','Greg_MovingSpot_1ms_Dir2',...
    'Greg_MovingSpot_11ms_Dir1','Greg_MovingSpot_11ms_Dir2'};
dataNum = length(dataName);
TrialNum = 2;

slct_MP_ind = [249, 111];
MP_label = {'Fingertip','Wrist'};

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

curr_fig = figure('Position',[60,160,1840,740],'Color','w');
colormap(jet(1000));
for MP_i = 1:length(slct_MP_ind)
    subplot(2,1,MP_i)
    y_slct = y_vib_sync(slct_MP_ind(MP_i),:)';
    spectrogram(y_slct,2048,2000,2048,Fs,'yaxis');
    caxis([-150 -90])
    ylim([0 1])
    
    xlabel({'Time (ms)';...
        [dataName{d_i},sprintf(' Trial%d -- ',t_i), MP_label{MP_i}]},...
        'Interpreter', 'none');
end
    drawnow;
    end
end