function savegif(savename, IM, fps)
% Function savegif is to save 3D imgags as gif animation
% 
% Input Variables:
% savename - gif name
% IM - input image frames (3D images: [row,cloumn,frameN]) 
% fps - frames per second
%
% Example:
% savename('cardiac_DMRI.gif',IM,10)     
% 
% Record of Revision
% Jun-05-2020===Zhao He===Original Code
% Jun-26-2020===Zhao He===normalized frame before imshow

% get 3D images dimension
[row,clom,frameN] = size(IM);

% show and save 3D images as gif animation
f1 = figure;
for i = 1:frameN
    
    % show per frame
    frame = normabs(squeeze(IM(:,:,i)));
    imagesc(frame);  axis image; axis off; colormap gray;   
    
    % show the number of current frame
    px = clom - 30; py = row - 10;
    text(px,py,[num2str(i)],'FontSize',12,'color','w'); 
    
    % get current frame     
    frame = getframe();
    im = frame2im(frame); 
    [I,map] = rgb2ind(im,256);  
    
    % save gif
    dt = 1/fps; % seconds per frame
	if i == 1
		imwrite(I,map,savename,'gif','WriteMode','overwrite','DelayTime',...
                 dt,'LoopCount',inf); % save first frame
	else
		imwrite(I,map,savename,'WriteMode','append','DelayTime',dt); % save other frames
    end
    
end

% close figure
close(f1);

end