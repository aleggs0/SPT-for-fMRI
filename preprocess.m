function [XX, XXorig] = preprocess(filename)
%PREPROCESS given a motion-corrected nii file, perform the remainder of the preprocessing
% (cubic trend removal, cropping, standardisation)

XXorig = double(niftiread(niftiinfo(filename)));
% filename='C:\Users\alexy\OneDrive\Summer Project\Code\data\othersources\ExBox11\fmri.nii';
% filename='C:\Users\alexy\OneDrive\Summer Project\Code\data\sub01018 (change)\restmcf.nii';
% files = gunzip('*.gz')

P=zeros(size(XXorig,1),size(XXorig,2),size(XXorig,3),4);
for i=1:size(XXorig,1)
    for j=1:size(XXorig,2)
        for k=1:size(XXorig,3)
            line=permute(XXorig(i,j,k,:),[4 1 2 3]);
            nums=(1:length(line))';
            p = polyfit(nums,double(line),3);
            P(i,j,k,:)=p;
            XXorig(i,j,k,:)=XXorig(i,j,k,:)-permute(p(end)+p(end-1)*nums+p(end-2)*nums.^2+p(end-3)*nums.^3, [4 3 2 1]);
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
XXorig=XXorig./std(XXorig,1,4);
XXorig=XXorig-mean(XXorig,4);
XXorig(isnan(XXorig))=0; %NaN arises from dividing by zero standard deviation

%%%Cropping
XX=XXorig(5:60,5:60,3:31,:);
end