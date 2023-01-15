function [start_change,end_change,TA,TB] = changepoints(eta_hat,Sigma_hat_inv)
%CHANGEPOINTS Given projected data, perform change-point detection
if nargin==1
    Sigma_hat_inv=diag(1./var(eta_hat,[],1));
end
[n,d]=size(eta_hat);
eta_hat_mean=1/n*sum(eta_hat,1);
S_summand=eta_hat-eta_hat_mean; %the j row is the jth summand

S_mat=zeros(n,n,d); %S((k1-1)/n,k2/n)=S_mat(k1,k2,:) permuted,
    %so that the diagonal entries of S_mat contain one term
    %though the diagonal of S is zero
for k1=1:n
    S_mat(k1,k1:end,:)=cumsum(S_summand(k1:end,:),1);
end

if isdiag(Sigma_hat_inv)
    %%%more efficient if Sigma_hat is diagonal
    w=diag(Sigma_hat_inv);
    objmat=tensorprod(S_mat.^2,w,3,1);
    objective_sum=sum(objmat,'all');
    [objective,idx]=max(objmat,[],'all');
    [start_change,end_change]=ind2sub([n,n],idx);
else
    objective=0;
    objective_sum=0;
    for k1=1:n-1
        for k2=k1+1:n
            summand=tensorprod(tensorprod(Sigma_hat_inv,S_mat(k1,k2,:),2,3),S_mat(k1,k2,:),1,3);
            objective_sum=objective_sum+summand;
            if summand>objective
                start_change=k1; end_change=k2;
                objective=summand;
            end
        end
    end
end

TA=objective_sum/n^3;
TB=objective/n;
end