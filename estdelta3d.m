function [objectives] = estdelta3d(XX,Ds,num_compons1)
%ESTDELTA3D Given data, finds cost for various bandwidths
%Set XX=data matrix; if data is fMRI, first apply load_preprocess to the
%motion-corrected nii file

if nargin<2
    Ds=0:6;
elseif length(Ds)<7
    Ds=0:6;
end
if nargin<3 
    num_compons1=5;
end
num_compons=num_compons1.^3;
[length_x, length_y, length_z, n]=size(XX);
if n~=225
    disp("Warning: this code only works for n=225")
end
objectives=zeros(size(Ds));
%diag_objectives=zeros(size(Ds));
sizes=zeros(size(Ds));
Asquarederrors=zeros(size(Ds));
%Bsquarederrors=zeros(size(Ds));
%Csquarederrors=zeros(size(Ds));
Achanges=zeros(size(Ds));
Bchanges=zeros(size(Ds));
Cchanges=zeros(size(Ds));
sepvar=zeros(size(Ds));

explainedvars=zeros(size(Ds));
explainedvar123s=zeros(size(Ds));
totalevals=zeros(size(Ds));
totaleval123s=zeros(size(Ds));
cumexplainedvars=zeros(num_compons,length(Ds));
allsortedevals=zeros(length_x*length_y*length_z,length(Ds));

spt1sum=zeros(length_x,length_x,length(Ds)); %these need to be divided by n
spt2sum=zeros(length_y,length_y,length(Ds));
spt3sum=zeros(length_z,length_z,length(Ds));
stsum=zeros(length(Ds),1);

previous_c1=zeros(length_x);previous_c2=zeros(length_y);previous_c3=zeros(length_z);
previous_symbol=zeros(length_x,length_y,length_z);

XX=XX-mean(XX,4);

disp('calculating symbols...')
symbol_times=zeros(length_x,length_y,length_z,n);
for t=1:n
    symbol_times(:,:,:,t)=fasttoeplitzsymbol(XX(:,:,:,t));
end
symbol_sum=sum(symbol_times,4);


disp('initial estimates...')
for Dno=1:length(Ds)
    D=Ds(Dno);
    stsum(Dno)=tensorprod(XX(1:end-D,1:end-D,1:end-D,:),XX(1+D:end,1+D:end,1+D:end,:),[1 2 3 4],[1 2 3 4]);
    spt1sum(:,:,Dno)=einsum(XX(:,1:end-D,1:end-D,:),XX(:,1+D:end,1+D:end,:),[2 3 4],[2 3 4]);
    spt2sum(:,:,Dno)=einsum(XX(1:end-D,:,1:end-D,:),XX(1+D:end,:,1+D:end,:),[1 3 4],[1 3 4]);
    spt3sum(:,:,Dno)=einsum(XX(1:end-D,1:end-D,:,:),XX(1+D:end,1+D:end,:,:),[1 2 4],[1 2 4]);
    
    st=stsum(Dno)/n;
    spt1=spt1sum(:,:,Dno)/n;
    spt2=spt2sum(:,:,Dno)/n;
    spt3=spt3sum(:,:,Dno)/n;
    c1_est=spt1;
    c1_est = (c1_est+c1_est')/2;
    [V1,D1] = eig(c1_est); evals1=diag(D1); D1pos=max(D1,0); c1_est=V1*D1pos*V1';
    c2_est=spt2/st;
    c2_est = (c2_est+c2_est')/2;
    [V2,D2] = eig(c2_est); evals2=diag(D2); D2pos=max(D2,0); c2_est=V2*D2pos*V2';
    c3_est=spt3/st;
    c3_est = (c3_est+c3_est')/2;
    [V3,D3] = eig(c3_est); evals3=diag(D3); D3pos=max(D3,0); c3_est=V3*D3pos*V3';
    

    %%%Now must use Toeplitz averaging
    Bsymbol=symbol_sum/n-toeplitzouterprod(c1_est,c2_est,c3_est);
    croppedBsymbol=Bsymbol(1:D,1:D,1:D); Bsymbol=zeros(length_x,length_y,length_z); Bsymbol(1:D,1:D,1:D)=croppedBsymbol;
    objectives(Dno)=objectives(Dno)+fastsquarednorm(c1_est,c2_est,c3_est,Bsymbol,...
        zeros(size(c1_est)),zeros(size(c2_est)),zeros(size(c3_est)),zeros(size(Bsymbol)));
    %diag_objectives(Dno)=norm(diag(c1_est).*diag(c2_est)'.*permute(diag(c3_est),[3 2 1])+Bsymbol(1,1,1),'fro');


    sizes(Dno)=objectives(Dno);
    Asquarederrors(Dno)=...
        fastsquarednorm(c1_est,c2_est,c3_est,zeros(length_x,length_y,length_z),zeros(length_x),zeros(length_y),zeros(length_z),zeros(length_x,length_y,length_z));
%     Bsquarederrors(Dno)=...
%         fastsquarednorm(zeros(length_x),zeros(length_y),zeros(length_z),Bsymbol,zeros(length_x),zeros(length_y),zeros(length_z),truesymbol*tau^2);
%     Csquarederrors(Dno)=...
%         fastsquarednorm(c1_est,c2_est,c3_est,Bsymbol,c1,c2,c3,truesymbol*tau^2);
    Achanges(Dno)=fastsquarednorm(c1_est,c2_est,c3_est,zeros(length_x,length_y,length_z),previous_c1,previous_c2,previous_c3,zeros(length_x,length_y,length_z));
    Bchanges(Dno)=fastsquarednorm(zeros(length_x),zeros(length_y),zeros(length_z),Bsymbol,previous_c1,previous_c2,previous_c3,previous_symbol);
    Cchanges(Dno)=fastsquarednorm(c1_est,c2_est,c3_est,Bsymbol,zeros(length_x),zeros(length_y),zeros(length_z),previous_symbol);
    sepvar(Dno)=trace(c1_est)*trace(c2_est)*trace(c3_est); %+length_x*length_y*length_z*Bsymbol(1,1,1);
    

    evals=evals1.*evals2'.*permute(evals3,[3 2 1]);
    [sortedevals,idx] = sort(evals(:),'descend');
    firstefns=zeros(size(XX,1),size(XX,2),size(XX,3),num_compons);
    firstevals=sortedevals(1:num_compons);
    for k=1:num_compons
        [i1,i2,i3] = ind2sub(size(XX(:,:,:,1)),idx(k));
        firstefns(:,:,:,k) = V1(:,i1).*V2(:,i2)'.*permute(V3(:,i3),[3 2 1]);
    end
    eta_hat=einsum(XX,firstefns,[1 2 3],[1 2 3]);
    [sortedevals1,idx1] = sort(evals1(:),'descend'); firstevals1=sortedevals1(1:num_compons1);
    firstefns1=V1(:,idx1(1:num_compons1));
    [sortedevals2,idx2] = sort(evals2(:),'descend'); firstevals2=sortedevals2(1:num_compons1);
    firstefns2=V2(:,idx2(1:num_compons1));
    [sortedevals3,idx3] = sort(evals3(:),'descend'); firstevals3=sortedevals3(1:num_compons1);
    firstefns3=V3(:,idx3(1:num_compons1));
    firstevals123=firstevals1.*firstevals2'.*permute(firstevals3,[3 2 1]);
    firstefns123=permute(firstefns1,[1 3 4 2]).*permute(firstefns2,[3 1 4 5 2]).*permute(firstefns3,[3 4 1 5 6 2]);
    eta_hat123=einsum(XX,firstefns123,[1 2 3],[1 2 3]);
    cumexplainedvars(:,Dno)=cumsum(var(eta_hat,[],1),2);

    explainedvars(Dno)=sum(var(eta_hat,[],1),2);
    explainedvar123s(Dno)=sum(var(eta_hat123,[],1),[2 3 4]);
    totalevals(Dno)=sum(firstevals(1:125));
    totaleval123s(Dno)=sum(firstevals123(1:5,1:5,1:5),'all');

    allsortedevals(:,Dno)=sortedevals;
    previous_c1=c1_est;previous_c2=c2_est;previous_c3=c3_est;previous_symbol=Bsymbol;
end
%%
disp("Performing 15-fold cross-validation:")
for m=1:15
    disp(m)
    start_leaveout=max(1,15*(m-1)-4); end_leaveout=min(15*m+5,225);
    leaveouts=start_leaveout:end_leaveout;
    symbolsumpartial=symbol_sum-sum(symbol_times(:,:,:,leaveouts),4);
    for Dno=1:length(Ds)
        D=Ds(Dno);
        st=(stsum(Dno)-einsum(XX(1:end-D,1:end-D,1:end-D,leaveouts),XX(1+D:end,1+D:end,1+D:end,leaveouts),...
            [1 2 3 4],[1 2 3 4]))/(n-length(leaveouts));
        spt1=(spt1sum(:,:,Dno)-einsum(XX(:,1:end-D,1:end-D,leaveouts),XX(:,1+D:end,1+D:end,leaveouts),...
            [2 3 4],[2 3 4]))/(n-length(leaveouts));
        spt2=(spt2sum(:,:,Dno)-einsum(XX(1:end-D,:,1:end-D,leaveouts),XX(1+D:end,:,1+D:end,leaveouts),...
            [1 3 4],[1 3 4]))/(n-length(leaveouts));
        spt3=(spt3sum(:,:,Dno)-einsum(XX(1:end-D,1:end-D,:,leaveouts),XX(1+D:end,1+D:end,:,leaveouts),...
            [1 2 4],[1 2 4]))/(n-length(leaveouts));
        c1_est=spt1; c1_est = (c1_est+c1_est')/2;
        c2_est=spt2; c2_est = (c2_est+c2_est')/2;
        c3_est=spt3/st/st; c3_est = (c3_est+c3_est')/2;
        [V1,D1] = eig(c1_est); D1pos=max(D1,0); c1_est=V1*D1pos*V1';
        [V2,D2] = eig(c2_est); D2pos=max(D2,0); c2_est=V2*D2pos*V2';
        [V3,D3] = eig(c3_est); D3pos=max(D3,0); c3_est=V3*D3pos*V3';
        %%%Now must use Toeplitz averaging
        Bsymbol=symbolsumpartial/(n-length(leaveouts))-toeplitzouterprod(c1_est,c2_est,c3_est);
        croppedBsymbol=Bsymbol(1:D,1:D,1:D); Bsymbol=zeros(length_x,length_y,length_z); Bsymbol(1:D,1:D,1:D)=croppedBsymbol;
        Qsymbol=zeros(length_x*2-1,length_y*2-1,length_z*2-1);
        Qsymbol(:,:,1:length_z)=[Bsymbol,flip(Bsymbol(:,2:end,:),2);...
            flip(Bsymbol(2:end,:,:),1), flip(flip(Bsymbol(2:end,2:end,:),1),2)];
        Qsymbol(:,:,length_z+1:end)=[flip(Bsymbol(:,:,2:end),3),flip(flip(Bsymbol(:,2:end,2:end),2),3);...
            flip(flip(Bsymbol(2:end,:,2:end),1),3), flip(flip(flip(Bsymbol(2:end,2:end,2:end,:),1),2),3)];
        for tt=15*(m-1)+1:15*m
            paddedXtt=zeros(2*length_x-1,2*length_y-1,2*length_z-1);
            paddedXtt(1:length_x,1:length_y,1:length_z)=XX(:,:,:,tt);
            Q_matX=ifftn(fftn(Qsymbol).*fftn(paddedXtt));
            B_matX=Q_matX(1:length_x,1:length_y,1:length_z);
            A_matX=tensorprod(tensorprod(tensorprod(XX(:,:,:,tt),c1_est,1,1),c2_est,1,1),c3_est,1,1);
            XCX=tensorprod(XX(:,:,:,tt),A_matX+B_matX,[1 2 3]);
            objectives(Dno)=objectives(Dno)-2*XCX/n;
            %diag_objectives(Dno)=diag_objectives(Dno)-2*tensorprod(XX(:,:,:,tt).^2,(diag(c1_est).*diag(c2_est)'.*permute(diag(c3_est),[3 2 1])+Bsymbol(1,1,1)),[1 2 3])/210; %n
        end
    end
end

%%
figure('Name',"Accuracy, true D="+"unknown" +", l_x="+length_x+", l_y="+length_y+", l_z="+length_z);
tiledlayout(2,2); set(gcf,'position',[20,60,1200,560]);

nexttile; plot(Ds,objectives); title('objectives'); ylim([min(objectives),2*median(objectives)-min(objectives)])
nexttile; plot(Ds,sizes); title('sizes'); ylim([min(sizes),2*median(sizes)-min(sizes)])
% nexttile; plot(Ds,sepvar); title('Variance from A, tr(A)'); ylim([min(sepvar),2*median(sepvar)-min(sepvar)])
% nexttile; plot(Ds(2:end),Achanges(2:end)); title('Change in A'); ylim([0,2*median(Achanges)])
% nexttile; plot(Ds(2:end),Bchanges(2:end)); title('Change in B'); ylim([0,2*median(Bchanges)])
% nexttile; plot(Ds(2:end),Cchanges(2:end)); title('Change in C'); ylim([0,2*median(Cchanges)])
nexttile; plot(Ds,explainedvars,Ds,explainedvar123s); title('Explained vars'); legend('first 125','5 in each dir')
nexttile; plot(Ds,totalevals,Ds,totaleval123s); title('Total evals'); legend('first 125','5 in each dir')
% nexttile; plot(Ds,Asquarederrors); title("A squared errors (true size "+...
%    norm(c1,'fro')*norm(c2,'fro')*norm(c3,'fro')+")"); ylim([min(Asquarederrors),2*median(Asquarederrors)-min(Asquarederrors)])
% nexttile; plot(Ds,Bsquarederrors); title("B squared errors (true size "+...
%     fastsquarednorm(zeros(length_x),zeros(length_y),zeros(length_z),tau^2*truesymbol,zeros(length_x),zeros(length_y),zeros(length_z),zeros(size(truesymbol)))+")");
% nexttile; plot(Ds,Csquarederrors); title("C squared errors (true size "+...
%     fastsquarednorm(c1,c2,c3,tau^2*truesymbol,zeros(length_x),zeros(length_y),zeros(length_z),zeros(size(truesymbol)))+")");
figure; plot(Ds,explainedvar123s); xlabel("D"); grid on;
figure; plot(allsortedevals(:,1:7));grid on; title('sorted evals');legend("D="+Ds(1:7)); ylim([0 20]); xlim([0 1000])
figure; plot(Ds,objectives); ylim([min(objectives),2*median(objectives)-min(objectives)]); xlabel("D");ylabel("objective values"); grid on
end