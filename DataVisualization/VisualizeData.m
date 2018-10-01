%% Visualize PDV data
% Created on 09/29/2018
% -------------------------------------------------------------------------
MapPath = './VisualMap/';

maskImg = imread([MapPath,'Dorsal_Mask.jpg']);
maskImg = rgb2gray(maskImg);
pointImg = imread([MapPath,'Dorsal_MP.jpg']);

% Compute the optical flow vector locations based on dot maps--------------
locator_num = 348; % Stream line locator number (74)
pointImg = pointImg(:,:,2); % Green points

% Low pass filter the chosen grid and computable area
lpF=ones(3)/9; % 3-by-3 mean filter
ptImg_PFilter=uint8(filter2(lpF,pointImg));
% OF_image(2:end-1,2:end-1)=OF_PFilter(2:end-1,2:end-1);
MP_Posi = findMP(ptImg_PFilter, locator_num);

MP_Posi = round(MP_Posi); % Round the position to get index

xDiff = diff(MP_Posi(:,2));
xDiff(abs(xDiff)<5) = 0;
MP_Score = cumsum([0; xDiff]).*10000 + MP_Posi(:,1);

[~,ind] = sort(MP_Score);
MP_Posi = MP_Posi(ind,:);

% figure
% imshow(pointImg)
% for i = 1:locator_num
%     text(MP_Posi(i,2),MP_Posi(i,1),num2str(i),'Color','r')
% end

tavgData = rms(PDV_data.Data,1);

Alpha = 50;

% figure
interp_radius = 100;
[r,c] = find(maskImg < 255);
interpImg = zeros(size(maskImg));

tic
wb_h = waitbar(0, 'O', 'Name','Interpolating Measurement');
interpNum = length(r);
for i = 1:interpNum
    ind = find((MP_Posi(:,1) < (r(i)+interp_radius)) &...
        (MP_Posi(:,1) > (r(i)-interp_radius)) &...
        (MP_Posi(:,2) < (c(i)+interp_radius)) &...
        (MP_Posi(:,2) > (c(i)-interp_radius)));
    if ~isempty(ind)
         waitbar(i/interpNum,wb_h,num2str(i));
%         imshow(pointImg)
%         hold on
%         plot(c(i),r(i),'*r')
%         plot(MP_Posi(ind,2),MP_Posi(ind,1),'og')
%         hold off
%         drawnow 

        d = ((MP_Posi(ind,1) - r(i)).^2 + (MP_Posi(ind,2) - c(i)).^2).^0.5;
        
        Phi=1./(d+Alpha);      
        interpImg(r(i),c(i))=(tavgData(ind)*Phi)./sum(Phi);        
    end
end
close(wb_h);

imagesc(interpImg)

toc



% pX = PDV_data.Posi(2,:);
% pY = PDV_data.Posi(3,:);
% pZ = PDV_data.Posi(4,:);
% 
% validInd = (abs(pX)<1) & (abs(pY)<1) & (abs(pZ)<1);
% pX = pX(validInd);
% pY = pY(validInd);
% pZ = pZ(validInd);
% 
% figure('Position',[20,50,1600,840])
% colormap(jet(1000))
% 
% tavgData = rms(PDV_data.Data(:,validInd),1);
% 
% scatter3(pX,pY,pZ,100,tavgData,'filled');
% view([80 70])




% % % % for i = 1:length(PDV_data.t)
% % % % currFrame = PDV_data.Data(1,validInd);
% % % % scatter3(pX,pY,pZ,100,currFrame,'filled');
% % % % axis equal
% % % % view([96 61])
% % % % 
% % % % text(0.12,0.02,0.02,sprintf('t = %.1f ms',PDV_data.t(i)*1000),...
% % % %     'FontSize',20,'Color','k')
% % % %     
% % % % drawnow
% % % % end

