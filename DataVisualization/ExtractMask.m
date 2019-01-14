img = imread('DigitIISinusoid/DigitIISinusoid_Mask.jpg');

lab_he = rgb2lab(img);

ab = lab_he(:,:,2:3);
ab = im2single(ab);
nColors = 2;
% repeat the clustering 10 times to avoid local minima
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',10,'Threshold',1e-10);


J = medfilt2(pixel_labels,[30,30]);

figure('Position',[60,160,1800,600]);
subplot(1,2,1)
imshow(pixel_labels,[])
title('Image Labeled by Cluster Index');
subplot(1,2,2)
imshow(J,[])