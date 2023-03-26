function [eta_hat123,firstevals123,c1_est,c2_est,c3_est] = sptprojections3d(XX,D,d1)
%SPTPROJECTIONS3D Project data XX onto the first d1*d1*d1 eigenfunctions found
%by SPT with bandwidth D
if nargin<2
    D=0;
end
if nargin<3
    d1=5;
end

XX=XX-mean(XX,4);
n=size(XX,4);
st = einsum(XX(1:end-D,1:end-D,1:end,:),XX(1+D:end,1+D:end,1:end,:),[1 2 3 4],[1 2 3 4])/n;

c1_est = einsum(XX(:,1:end-D,1:end,:),XX(:,1+D:end,1:end,:),[2 3 4],[2 3 4])/n;
c2_est = einsum(XX(1:end-D,:,1:end,:),XX(1+D:end,:,1:end,:),[1 3 4],[1 3 4])/n/st;
c3_est = einsum(XX(1:end-D,1:end-D,:,:),XX(1+D:end,1+D:end,:,:),[1 2 4],[1 2 4])/n/st;

sym1size=norm(c1_est+c1_est','fro');antisym1size=norm(c1_est-c1_est','fro'); c1symscore=(sym1size-antisym1size)/(sym1size+antisym1size);

c1_est = (c1_est+c1_est')/2;
[V1,D1] = eig(c1_est);
D1pos=max(D1,0);
evals1=diag(D1);
c2_est = (c2_est+c2_est')/2;
[V2,D2] = eig(c2_est);
D2pos=max(D2,0);
evals2=diag(D2);
c3_est = (c3_est+c3_est')/2;
[V3,D3] = eig(c3_est);
D3pos=max(D3,0);
evals3=diag(D3);

c1_est=V1*D1pos*V1';c2_est=V2*D2pos*V2';c3_est=V3*D3pos*V3';

[sortedevals1,idx1] = sort(evals1(:),'descend'); firstevals1=sortedevals1(1:d1);
firstefns1=V1(:,idx1(1:d1));
[sortedevals2,idx2] = sort(evals2(:),'descend'); firstevals2=sortedevals2(1:d1);
firstefns2=V2(:,idx2(1:d1));
[sortedevals3,idx3] = sort(evals3(:),'descend'); firstevals3=sortedevals3(1:d1);
firstefns3=V3(:,idx3(1:d1));
firstevals123=firstevals1.*firstevals2'.*permute(firstevals3,[3 2 1]);
firstefns123=permute(firstefns1,[1 3 4 2]).*permute(firstefns2,[3 1 4 5 2]).*permute(firstefns3,[3 4 1 5 6 2]);
eta_hat123=einsum(XX,firstefns123,[1 2 3],[1 2 3]);
end