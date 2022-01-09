function [mse,rmse,nmse,psnr,stdMSE,stdRMSE, stdNMSE, stdPSNR] = MSE_RMSE_NMSE_PSNR(ref, X)
% Function MSE_RMSE_NMSE_PSNR is to calculate the cofficients: MSE(Mean Squared Error), RMSE (Root
% Mean Squared Error), NMSE (Normalized Mean Squared Error)and PSNR (Peak Sginal-to-Noise Ratio), for
% qualitatively evaluating the image quality.
%
% Input Variables:
% ref - reference image , or called the ground truth 
% X - image to be evaluated
%
% Output Variables:
% mse - Mean Square Error
% rmse - Root Mean Square Error
% nmse - Normalized Mean Square Error
% psnr - Peak Sginal-to-Noise Ratio
% stdMSE - standard deviation of MSE (only for 3D case)
% stdRMSE - standard deviation of RMSE (only for 3D case)
% stdNMSE - standard deviation of NMSE (only for 3D case)
% stdPSNR - standard deviation of PSNR (only for 3D case)

%
% Record of Revisions:
% Jul-31-2020===Zhao He===original code
% Aug-04-2020===Zhao He===add NMSE coffficient


n = ndims(ref); 
ref = normabs(ref)*255;
X = normabs(X)*255;

if n==2 % if 2D image
     mse = sum(sum((ref-X).^2))/size(ref,1)/size(ref,2); 
     rmse = sqrt(mse); 
     nmse = sum(sum((ref-X).^2))/sum(sum(ref.^2));
     psnr = 10*log10(255^2/mse);
elseif n==3 % if 3D image 
     for ii = 1:size(ref,3)
         E(ii)= sum(sum((ref(:,:,ii)-X(:,:,ii)).^2))/size(ref,1)/size(ref,2); 
         NE(ii) = sum(sum((ref(:,:,ii)-X(:,:,ii)).^2))/sum(sum(ref(:,:,ii).^2));
     end
     Ep = 10*log10(255^2 ./E);  
     mse = mean(E); stdMSE = std(E);
     rmse = mean(sqrt(E)); stdRMSE = std(sqrt(E)); 
     nmse = mean(NE); stdNMSE = std(NE);
     psnr = mean(Ep); stdPSNR = std(Ep);
else
     disp('Input must be 2D or 3D data!');return;
end
    
end


function imgNorm = normabs (imgIn)
imgIn = abs(imgIn);
imgNorm = (imgIn - min(imgIn(:))) / (max(imgIn(:))-min(imgIn(:)));
end
