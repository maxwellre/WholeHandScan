%% Display signals used for pilot study
% Created on 01/17/2019
% -------------------------------------------------------------------------
if ~exist('outSig','var')
    load('WHC2019FreqExpSignals.mat');
end


%% Spectrogram
discardLen = round(pauseTimeInSec*Fs)-2000;

fig_h = figure('Position',[60 120 800 860],'Color','w');
colormap(flipud(hot(1000)));
for i = 1:4
    subplot(4,1,i)
    switch i
        case 1
            temp = outSig.sigA(discardLen:(end-discardLen));
        case 2
            temp = outSig.sigB(discardLen:(end-discardLen));
        case 3
            temp = outSig.rsA(discardLen:(end-discardLen));
        case 4
            temp = outSig.rsB(discardLen:(end-discardLen));
    end
    spectrogram(temp,1024,1000,1024,Fs,'yaxis'); ylim([0 0.8]);
    box off;
end
fprintf('Colorbar range: %.1f - %.1f\n',caxis);

%% Compact plot of signals
% zeropadLen = 100; % Padding 100 zeros on each side of the signal
% t_break = 0.01; % Time break between signals
% 
% fig_h = figure('Position',[60 120 1800 860],'Color','w');
% hold on
% t0 = 0;
% for i = 1:numSigs     
%     temp = [zeros(1,zeropadLen),sigSeg{i,1},zeros(1,zeropadLen)];
%     sigLen = length(temp);
%     t = t0 +(1:sigLen)/Fs;
%     plot(t,temp);
%     
%     t0 = t(end) +t_break;
% end
% xlabel('Time (s)'); ylabel('Amplitude (V)');
% set(gca,'FontSize',16,'LineWidth',0.75);

%%
% % % Save your figure as .pdf file in order to import into illustrator
% print(fig_h,'NewFigure','-dpdf','-bestfit','-r600','-painters')