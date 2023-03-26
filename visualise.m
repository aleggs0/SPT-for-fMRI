%generates video of a slice having preprocessed
figure; set(gcf,'position',[20,60,900,600]);
axis equal; xlabel('y'); ylabel('x'); zlabel('z');
x = 9:56; y = 5:56; z = 1:33;
xslice = [32];
yslice = [30];
zslice = [15];
for t=1:225
    slice(y,x,z,XXorig(x,y,z,t),yslice,xslice,zslice); clim([0 1500]); drawnow
    F(t) = getframe(gcf);
end

% create the video writer with 1 fps
writerObj = VideoWriter('visualised.avi');
writerObj.FrameRate = 10;

% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

% layer=permute(XXorig(:,:,15,:),[1 2 4 3]);
% figure; %%%
% for i = 1:size(layer,3)
%     imagesc([layer(:,:,i)],[0,1345]); axis equal; colorbar;
%     F(i) = getframe(gcf);
%     drawnow
% end