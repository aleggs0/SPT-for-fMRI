%inputs: XX, etas,efns,D,evals
figure('Name',"D="+D+""); tiledlayout(2,3); set(gcf,'position',[20,60,1200,560]);
k3=1;
for k2=1:2
    for k1=1:3
        nexttile;
        V=firstefns123(:,:,:,k1,k2,k3);

        s = isosurface(V,-0.006);
        p = patch(s);
        isonormals(V,p)
        view(3);
        set(p,'FaceColor',[0.5 1 0.5]); 
        set(p,'EdgeColor','none');
        camlight;

        s = isosurface(V,0.006);
        p = patch(s);
        isonormals(V,p)
        view(3);
        set(p,'FaceColor',[1 0.5 0.5]); 
        set(p,'EdgeColor','none');
        camlight;
        
        title("Efn ("+k1+","+k2+","+k3+") : \lambda="+sprintf('%#.3g',firstevals123(k1,k2,k3)));
        xlim([1 56]); ylim([1 56]); zlim([1 29]);
    end
end

[start_change,end_change,TA,TB] = changepoints(reshape(eta_hat123,n,d1^3));

figure('Name',"D="+D+" scores"); tiledlayout(3,5); set(gcf,'position',[20,60,1200,560]);
tiledlayout(2,3); set(gcf,'position',[20,60,1200,560]);
if max(abs(eta_hat123),[],'all')<30
    lim=20; lims=[-lim lim];
elseif max(abs(eta_hat123),[],'all')<120
    lim=80; lims=[-lim lim];
else
    lim=[];
end
k3=1;
for k2=1:2
    for k1=1:3
        nexttile; plot(eta_hat123(:,k1,k2,k3)); hold on;
        title("Efn ("+k1+","+k2+","+k3+") : \lambda="+sprintf('%#.3g',firstevals123(k1,k2,k3)));
        if isempty(lim)
            lims=ylim;
        else
            ylim([-lim lim]) %comment this if unprocessed data
        end
        if TA>13
            plot([start_change-0.5 start_change-0.5],[-lim lim],[end_change+0.5 end_change+0.5],[-lim lim])
        end
        xlim([0 225])
    end
end

figure('Name',"D="+D+" covariance matrices"); tiledlayout(1,3); set(gcf,'position',[20,60,1200,560]);
nexttile, imagesc(c1_est), axis square, title("C_1 (covariance in x-direction)"), maxc=max(c1_est,[],'all'); set(gca,'CLim',[-0.2*maxc maxc]);
nexttile, imagesc(c2_est), axis square, title("C_2 (covariance in y-direction)"), maxc=max(c2_est,[],'all'); set(gca,'CLim',[-0.2*maxc maxc]);
nexttile, imagesc(c3_est), axis square, title("C_3 (covariance in z-direction)"), maxc=max(c3_est,[],'all'); set(gca,'CLim',[-0.2*maxc maxc]);