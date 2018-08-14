%% Get location of the laser point from a camera shot
% Note: the unit of estimated location of the laser point is determined by
% the perspective mapping matrix.
%--------------------------------------------------------------------------
% laserLoc = getLaserPosi(img, laserPointSize, mapping, dispImg)
%
% Inputs:
% 1. 'img' - a string contain the path of the camera shot image, or a
%            tensor containing the image.
%                   
% 2. 'laserPointSize' - size of the laser point in the image (pixel^2).
% 
% 3. 'mapping' - a 3-by-3 matrix defining perspective mapping of the
%                camera.                
%
% Optional inputs:
% 4. 'dispImg' - display the location estimation result
%
% Outputs:
% 1. 'laserLoc' - a 2-by-1 vector containing the location of the laser 
%                 point.
%--------------------------------------------------------------------------
function laserLoc = getLaserPosi(img, laserPointSize, mapping, dispImg)
% Author: Yitian Shao
% Created on 05/01/2018
%==========================================================================
if nargin < 4
    dispImg = 0;
end
if ischar(img)    
    J1 = imread(img);   
else
    J1 = img;
end

J2 = rgb2gray(J1);   
[~,ind] = sort(J2(:),'descend');
[y_i,x_i] =ind2sub(size(J2),ind);

pixelNum = length(x_i);
k = 0;
nearDist = 0;
while (k < pixelNum) && (nearDist < laserPointSize)
    k = k +1;
    nearDist = (y_i(k)-y_i(k+1))^2+(x_i(k)-x_i(k+1))^2;   
end

Laser_pix_loc = [mean(x_i(1:k));mean(y_i(1:k));1];
laserLoc = mapping*Laser_pix_loc;
laserLoc = laserLoc(1:2);

if dispImg
    imshow(J2);
    hold on;
    plot((x_i(1:k)),(y_i(1:k)),'.g');
    plot(Laser_pix_loc(1),Laser_pix_loc(2),'or');
    hold off
end