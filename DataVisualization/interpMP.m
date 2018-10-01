function interpImg = interpMP(maskImg, MP_Posi, tavgData,...
    maskThreshold, interp_radius, Alpha, Range)
% Created on 10/01/2018
% -------------------------------------------------------------------------
if nargin < 4
    maskThreshold = 255;
end

if nargin < 5
    interp_radius = 30;
end

if nargin < 6
    Alpha = 50;
end

if nargin < 7
    Range = 100;
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
        d = ((MP_Posi(ind,1) - r(i)).^2 + (MP_Posi(ind,2) - c(i)).^2).^0.5;
        
%         Phi=1./(d+Alpha);    
        Phi=1./(d+Alpha) - 1./(Range+Alpha);  
        validInd = find(Phi>0);
        Phi = Phi(validInd);
        
        interpImg(r(i),c(i))=(tavgData(ind(validInd))*Phi)./sum(Phi);        
    end
end
