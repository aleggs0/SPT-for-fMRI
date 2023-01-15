function [eta_hat,sortedevals] = pcaprojections3d(XX,d)
%PCAPROJECTIONS3D Project data onto the first d principal components
if nargin<2
    d=15;
end
XX=XX-mean(XX,4);

Z=reshape(XX,[size(XX,1)*size(XX,2)*size(XX,3),size(XX,4)]);
[Vtime,~] = eig(Z'*Z);
    %outer product in time dimension, summing over both spatial dimensions
efns=Z*Vtime./vecnorm(Z*Vtime);
evals=vecnorm(Z*(Z'*efns))./size(Z,2);
[sortedevals,idx]=sort(evals,'descend');
sortedefns=zeros(size(efns));
for k=1:length(evals)
    sortedefns(:,k)=efns(:,idx(k));
end
firstefns=zeros(size(XX,1),size(XX,2),size(XX,3),d);
firstevals=zeros(d,1);
for k=1:d
    firstevals(k)=sortedevals(k);
    firstefns(:,:,:,k)=reshape(sortedefns(:,k),[size(XX,1) size(XX,2) size(XX,3)]);
end
eta_hat=einsum(XX,firstefns,[1 2 3],[1 2 3]);