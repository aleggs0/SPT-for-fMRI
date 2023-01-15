%process one subject, assumed already motion corrected and unzipped
filename = 'C:\Users\alexy\OneDrive\Summer Project\Code\data\sub00440\rest_mcf.nii';
figurename = "sub00440";
disp(figurename)
num_compons1=5;
Ds=[0 1 2 3 4];
XX=preprocess(filename);
[length_x, length_y, length_z, n]=size(XX);
TAs=zeros(1,length(Ds)+1);
c1=zeros(size(XX,1));c2=zeros(size(XX,2));c3=zeros(size(XX,3));
%%%
figure('Name',figurename);
tiledlayout(2,3); set(gcf,'position',[20,60,900,450]);

eta_hat = pcaprojections3d(XX,num_compons1^3);
[start_change,end_change,TA,TB] = changepoints(eta_hat);
TAs(1)=TA;
for Dno=1:length(Ds)
    D=Ds(Dno);
    eta_hat123=sptprojections3d(XX,D,num_compons1);
    [start_change,end_change,TA,TB] = changepoints(reshape(eta_hat123,n,num_compons1^3));
    TAs(Dno+1)=TA;
end
disp(TAs)

%estimating delta takes time and is optional (can just choose D=2)
[objectives] = estdelta3d(XX,Ds,num_compons1);
disp(objectives)