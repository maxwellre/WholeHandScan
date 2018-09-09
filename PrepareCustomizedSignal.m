
Fs = 5000;
dt = 1/Fs;

t = dt:dt:10;

sine_freq = 100;

sinewave = sin(2*pi*sine_freq.*t);

sinewave = (sinewave - min(sinewave));

[ spect, freq ] = spectr(sinewave, Fs);

subplot(2,1,1);
plot(t,sinewave);
subplot(2,1,2);
plot(freq,spect);

% Output signal as .csv file
file_id = fopen(sprintf('sine%dhz_Uni.txt',sine_freq), 'w');
fprintf(file_id, '%f\n', sinewave);
fclose(file_id);