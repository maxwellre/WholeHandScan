%% Prepare sawtooth and impulse train for the actuator array
%--------------------------------------------------------------------------
function [ST_sig, Imp_sig] = sawtoothImpulse(fn,Fs,P2PAmp,trial_T,win_len,...
    is_disp)
%--------------------------------------------------------------------------
% Author: Yitian Shao (syt@ece.ucsb.edu)
% Created on 06/22/2018 
%--------------------------------------------------------------------------
if nargin < 6
    is_disp = 0;
end
%--------------------------------------------------------------------------
Ts = 1/Fs; % Sampling duration
f_num = length(fn);
t = 0:Ts:trial_T;

trial_len = length(t);

ST_sig = []; % Sawtooth
for i = 1:f_num
    x = P2PAmp * (sawtooth(2*pi*fn(i)*t))';
    x((win_len+1):end) = 0;
    ST_sig = [ST_sig;x];
end

pulseRatio = 0.95*2-1; % Pulse width: Sawtooth >95% 
Imp_sig = zeros(size(ST_sig));
Imp_sig(ST_sig > (pulseRatio*P2PAmp)) = P2PAmp;
Imp_sig(ST_sig <= (pulseRatio*P2PAmp)) = -P2PAmp;

if is_disp
    figure('Position',[80,80,1600,800]);
    plot(Ts:Ts:length(ST_sig)*Ts,ST_sig);
    hold on
    plot(Ts:Ts:length(Imp_sig)*Ts,Imp_sig);
    hold off
    xlabel('Time (secs)')
    ylabel('Amplitude (Volt)')
    xlim([0 length(ST_sig)*Ts])
end