function [XX, XXraw,XXuncropped] = preprocess(filename)
%PREPROCESS given a motion-corrected nii file, perform the remainder of the preprocessing
% (cubic trend removal, cropping, standardisation)

XXraw = double(niftiread(niftiinfo(filename)));
XXuncropped=XXraw;
% filename='C:\Users\alexy\OneDrive\Summer Project\Code\data\othersources\ExBox11\fmri.nii';
% filename='C:\Users\alexy\OneDrive\Summer Project\Code\data\sub01018 (change)\restmcf.nii';
% files = gunzip('*.gz')

P=zeros(size(XXuncropped,1),size(XXuncropped,2),size(XXuncropped,3),4);
for i=1:size(XXuncropped,1)
    for j=1:size(XXuncropped,2)
        for k=1:size(XXuncropped,3)
            line=permute(XXuncropped(i,j,k,:),[4 1 2 3]);
            nums=(1:length(line))';
            p = polyfit(nums,double(line),3);
            P(i,j,k,:)=p;
            XXuncropped(i,j,k,:)=XXuncropped(i,j,k,:)-permute(p(end)+p(end-1)*nums+p(end-2)*nums.^2+p(end-3)*nums.^3, [4 3 2 1]);
        end
    end
end
%visualise coefficients in cubic trend removal
% figure; tiledlayout(2,2);
% nexttile; imagesc(P(:,:,15,4)); axis square; colorbar;
% nexttile; imagesc(P(:,:,15,3)); axis square; colorbar;
% nexttile; imagesc(P(:,:,15,2)); axis square; colorbar;
% nexttile; imagesc(P(:,:,15,1)); axis square; colorbar;

%%%standardisation step
% layer=permute(Y(:,:,ceil(size(Y,3)/2),:), [1 2 4 3]);
% figure; imagesc(mean(layer,3),[0,1345]); axis square; colorbar;
% figure; imagesc(std(layer,0,3), [0 100]); axis square; colorbar;
% figure; imagesc(std(Y(:,:,15,:),0,4), [0 100]); axis square; colorbar;
XXuncropped=XXuncropped./std(XXuncropped,1,4);
XXuncropped=XXuncropped-mean(XXuncropped,4);
XXuncropped(isnan(XXuncropped))=0; %NaN arises from dividing by zero standard deviation

%%%Cropping
XX=XXuncropped(5:60,5:60,3:31,:);
end