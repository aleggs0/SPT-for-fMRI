function [XX,c1,c2,c3] = simulate(sizeXX,p,tau,first_change,last_change,change_size,change_direction)
%SIMULATE data XX of size sizeXX
%   XX is the sum of a separable part with covariance c1,c2,c3, a banded
%   part with bandwidth d=2p+1, and an epidemic change
if nargin<1
    length_x=56; length_y=56; length_z=29; n=225;
else
    length_x=sizeXX(1);
    length_y=sizeXX(2);
    length_z=sizeXX(3);
    n=sizeXX(4);
end
if nargin<3
    p=0; tau=0;
end
if nargin<6
    first_change=1; last_change=1; change_size=0;
end
if nargin<7
    changex=normrnd(0,1,length_x+4,1);
    changey=normrnd(0,1,1,length_y+4);
    changez=normrnd(0,1,1,1,length_z+4);
    changex=changex(1:end-4)+changex(2:end-3)+changex(3:end-2)+changex(4:end-1)+changex(5:end);
    changey=changey(1:end-4)+changey(2:end-3)+changey(3:end-2)+changey(4:end-1)+changey(5:end);
    changez=changez(1:end-4)+changez(2:end-3)+changez(3:end-2)+changez(4:end-1)+changez(5:end);
    change_direction=changex(1:length_x).*changey(1:length_y).*changez(1:length_z);
%     change_direction=normrnd(0,1,length_x,length_y,length_z);
end
d1=13;
d2=13;
d3=13;

degs1=0:d1-1;
degs1 = degs1(ones(1,length_x),:);
points1=(-1:2/(length_x-1):1)';
points1 = points1(:,ones(1,d1));
trueV1=legendreP(degs1,points1);
degs2=0:d2-1;
degs2 = degs2(ones(1,length_y),:);
points2=(-1:2/(length_y-1):1)';
points2 = points2(:,ones(1,d2));
trueV2=legendreP(degs2,points2);
degs3=0:d3-1;
degs3 = degs3(ones(1,length_z),:);
points3=(-1:2/(length_z-1):1)';
points3 = points3(:,ones(1,d3));
trueV3=legendreP(degs3,points3);
trueevals1=(1:d1).^(-2);%
trueevals2=(1:d2).^(-1.8);%
trueevals3=(1:d3).^(-1.6);%
trueV1=trueV1./vecnorm(trueV1);
trueV2=trueV2./vecnorm(trueV2);
trueV3=trueV3./vecnorm(trueV3);
is=(-p:p)'; js=-p:p; ks=permute(-p:p,[1 3 2]);
% Q=(-1).^abs(is+js+ks);
Q=9/16*(1-abs(is)./(p+1)).*(1-abs(js)./(p+1)).*(1-abs(ks)./(p+1));

epsilon=normrnd(0,1,d1,d2,d3,n);
XX=einsum( einsum( einsum((trueevals1.^0.5).*trueV1,epsilon,'ip,pqrt->iqrt')...
        ,trueevals2.^0.5.*trueV2,'iqrt,jq->ijrt'),...
        trueevals3.^0.5.*trueV3,'ijrt,kr->ijkt');
epsilon=normrnd(0,1,length_x+2*p,length_y+2*p,length_z+2*p,n);
W=zeros(length_x,length_y,length_z,n);
for i=1:length_x
    for j=1:length_y
        for k=1:length_z
            W(i,j,k,:)=sum(Q.*epsilon(i:i+2*p,j:j+2*p,k:k+2*p,:),[1 2 3]);
        end
    end
end
XX=XX+tau*W;
change=change_size*change_direction./norm(change_direction(:));
XX(:,:,:,first_change:last_change)=XX(:,:,:,first_change:last_change)+change;

c1=trueV1*diag(trueevals1)*trueV1'; 
c2=trueV2*diag(trueevals2)*trueV2';
c3=trueV3*diag(trueevals3)*trueV3';
end