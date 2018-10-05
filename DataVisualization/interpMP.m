function interpImg = interpMP(maskImg, MP_Posi, tavgData,...
    maskThreshold, interp_radius, Alpha, C, px2mm)
% Created on 10/01/2018
% -------------------------------------------------------------------------
if nargin < 4
    maskThreshold = 255;
end

if nargin < 5
    interp_radius = 30; % (pixel)
end

if nargin < 6
    Alpha = 25.5; % (mm)
end

if nargin < 7
    C = 0.087;
end

if nargin < 8
    px2mm = 0.289;  % (mm/pixel)
end

[r,c] = find(maskImg < maskThreshold);
interpImg = NaN(size(maskImg));

interpNum = length(r);
for i = 1:interpNum
    ind = find((MP_Posi(:,1) < (r(i)+interp_radius)) &...
        (MP_Posi(:,1) > (r(i)-interp_radius)) &...
        (MP_Posi(:,2) < (c(i)+interp_radius)) &...
        (MP_Posi(:,2) > (c(i)-interp_radius)));
    if ~isempty(ind)
        d = px2mm*(((MP_Posi(ind,1) - r(i)).^2 +...
            (MP_Posi(ind,2) - c(i)).^2).^0.5);
           
        Phi=17./(d+Alpha) - C;  
        validInd = find(Phi>0);
        Phi = Phi(validInd);
        
        interpImg(r(i),c(i))=(tavgData(ind(validInd))*Phi)./sum(Phi);        
    end
end
