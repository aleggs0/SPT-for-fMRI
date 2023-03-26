text = fileread('/home/jada2/Data_Store/Warwick_Data/aston1/Brain_Imaging/Beijing_Zang/Beijing_Zang_subjects.txt');
subjectnames = split(text); subjectnames=subjectnames(1:end-1);
TA_pca=zeros(length(subjectnames),1);
TA_pt=zeros(length(subjectnames),1);
TA_spt=zeros(length(subjectnames),1);
D=4;
destfolder="/mhome/damtp/r/ay343/code/";
destname=destfolder+"rest_mcf.nii";
delete(destname)
num_compons_pca=12;
num_compons1=5;
n=225;
% for subjectno=[25 30 32 57 107 119]
% 128 missing
% 25 skull noise
% 30,32,57 corrupted
% 31 somewhat corrupted
% 119 weird shake
for subjectno=1:length(subjectnames)
    subjectname=string(subjectnames(subjectno));
    disp(subjectname)
    pause(1) %prevent matlab from crashing
    sourcename="/home/jada2/Data_Store/Warwick_Data/aston1/Brain_Imaging/Beijing_Zang/"+subjectname+"/func/rest_mcf.nii";
    if isfile(sourcename)
        copyfile(sourcename,destname)
    elseif subjectname=="sub62966"
        continue
    else
        gunzip(sourcename+".gz",destfolder)
    end
    [XX, ~] = preprocess(destname);
    %%%
    [eta_hat,sortedevals] = pcaprojections3d(XX,num_compons_pca);
    [~,~,TA,~] = changepoints(eta_hat);
    TA_pca(subjectno)=TA;
    eta_hat123=sptprojections3d(XX,0,num_compons1);
    [~,~,TA,~] = changepoints(reshape(eta_hat123,n,num_compons1^3));
    TA_pt(subjectno)=TA;
    eta_hat123=sptprojections3d(XX,D,num_compons1);
    [~,~,TA,~] = changepoints(reshape(eta_hat123,n,num_compons1^3));
    TA_spt(subjectno)=TA;
    %%%
    delete(destname)
end
TAall=[TA_pca, TA_pt, TA_spt];
save("data_"+string(today("datetime")),"TAall");
%%%%TAall=TAall(setdiff(1:198,[128, 30,32,57, 31]),:)