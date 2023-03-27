%process one subject, assumed already motion corrected and unzipped
showfigs=true;
filename = 'C:\Users\alexy\OneDrive\Summer Project\Code\data\sub69518 (no change)\rest_mcf.nii';
figurename = "sub69518";
disp(figurename)
num_compons1=5;
num_compons_pca=15;
Ds=[0 1 2 3 4 5 6 7 8];
[XX,XXorig]=preprocess(filename);
[length_x, length_y, length_z, n]=size(XX);
TAs=zeros(1,length(Ds)+1);
start_changes=zeros(1,length(Ds)+1);
end_changes=zeros(1,length(Ds)+1);
c1=zeros(size(XX,1));c2=zeros(size(XX,2));c3=zeros(size(XX,3));
%%%

[eta_hat,sortedevals] = pcaprojections3d(XX,num_compons_pca);
[start_change,end_change,TA,TB] = changepoints(eta_hat);
%changedirection=mean(XX(:,:,:,start_change:end_change),4)-mean(XX,4);
TAs(1)=TA;
start_changes(1)=start_change;
end_changes(1)=end_change; 
for Dno=1:length(Ds)
    D=Ds(Dno);
    [eta_hat123,sortedevals]=sptprojections3d(XX,D,num_compons1,showfigs&&(D==0||D==4));
    [start_change,end_change,TA,TB] = changepoints(reshape(eta_hat123,n,num_compons1^3));
    TAs(Dno+1)=TA;
    start_changes(Dno+1)=start_change;
    end_changes(Dno+1)=end_change;
end
disp(TAs)

%estimating delta takes time and is optional (can just choose D=2)
[objectives] = estdelta3d_mod(XX,Ds,num_compons1);
disp(objectives)