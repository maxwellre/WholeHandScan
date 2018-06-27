%% Prepare discrete sine sweep for the actuator array
%--------------------------------------------------------------------------
function DSW_sig = discreteSineSweep(fn,Fs,sineP2PAmp,trial_T,win_len,...
    is_disp)
%--------------------------------------------------------------------------
% Author: Yitian Shao (syt@ece.ucsb.edu)
% Created on 06/19/2018 
%--------------------------------------------------------------------------
if nargin < 6
    is_disp = 0;
end
%--------------------------------------------------------------------------
Ts = 1/Fs; % Sampling duration
f_num = length(fn);
t = 0:Ts:trial_T;

bin_num = win_len/2+1;
trial_len = length(t);

DSW_sig = []; % Discrete sine sweep signal
for i = 1:f_num
    x = (sineP2PAmp*sin(2*pi*fn(i)*t))';
    cross_ind = find(abs(x(2:win_len)) < 0.01);
    x((cross_ind(end)+1):end) = 0;
    DSW_sig = [DSW_sig;x];
end

%--------------------------------------------------------------------------
% Validate the discrete sine sweep stimuli
if is_disp
figure('Position',[80,80,1600,800]);
subplot(2,1,1)
plot(0:Ts:Ts*(length(DSW_sig)-1),DSW_sig)
xlim([0 Ts*length(DSW_sig)])
xlabel('Time (secs)')
ylabel('Amplitude (Volt)')
subplot(2,1,2)
hold on
for i = 1:f_num
    ind0 = (i-1)*trial_len;
    
    Y_fft = fft(DSW_sig(ind0+(1:win_len)),win_len)/win_len; % FFT 
    f = Fs/2*linspace(0,1,bin_num); % Frequency scale
    abs_Y_fft = 2*abs(Y_fft(1:bin_num)); % Spectrum 
    
    plot(f,abs_Y_fft)
end
xlim([fn(1) fn(end)+1])
xlabel('Frequency (Hz)')
ylabel('Amplitude (Volt)')
end
end