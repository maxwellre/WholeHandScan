%% Convert PDV data from ASCII to .mat file
% Created on 09/28/2018
% -------------------------------------------------------------------------
% Path_Source = 'TapDigitII_Volar_SBJ1/';
Path_Source = 'TapDigitII_Dorsal_SBJ1/';

Fs = 20000; % Expected Sampling Frequency (Hz)
dt = 1/Fs; % Expected Sampling duration (secs)

strInd = regexp(Path_Source,'_+');
condiLabel = Path_Source(1:strInd(1)-1);
handDirection = Path_Source(strInd(1)+1:strInd(2)-1);
disp([condiLabel,' - ',handDirection])

path_info = dir([Path_Source,'*.txt']);
path_len = length(path_info);

PDV_data.Data = [];
PDV_data.Posi = zeros(4,path_len);

wb_h = waitbar(0, 'O', 'Name','Converting PDV data from ASCII to Matrix');

for i = 1:path_len
    a_name = path_info(i).name;
    
    numInd = regexp(a_name,'\d');
    
    p_i = str2double(a_name(numInd));
    
    fileID = fopen([Path_Source,a_name],'r');
    dataInfo.Source = fgetl(fileID);
    dataInfo.Position = fgetl(fileID);
    dataInfo.Component = fgetl(fileID);
    dataInfo.SignalType = fgetl(fileID);
    fgetl(fileID);
    dataInfo.Column = fgetl(fileID);
    dataInfo.Unit = fgetl(fileID);
    a_data = textscan(fileID,'%s');

    posiStr = regexp(dataInfo.Position,'-?+\d+\.?\d*','match');
    for j = 1:4
        PDV_data.Posi(j,i) = str2double(posiStr{j});
    end
    
    waitbar(i/path_len,wb_h,sprintf('Sample %d [%d match %d]',...
        i,p_i,PDV_data.Posi(1,i)));
    
    PDV_data.Data = [PDV_data.Data,cellfun(@str2double,a_data{1}(2:2:end))];
    
    fclose(fileID);
end
close(wb_h);

PDV_data.t = cellfun(@str2double,a_data{1}(1:2:end));

actual_dt = mean(diff(PDV_data.t));
if (actual_dt - dt) > 1e-10
    error('Sampling Frequency Error > 1e-10!')
end

save(sprintf('%s_%s.mat',condiLabel,handDirection),'PDV_data','dataInfo');

