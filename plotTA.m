TAorig=TAD0(keptindices); TAalt1=TAD2(keptindices); TAalt2=TAD4(keptindices);
[TAsorted,idx1]=sort(TAorig);
figure; plot(TAsorted); hold on; scatter(1:length(keptindices),TAalt1(idx1),'x'); scatter(1:length(keptindices),TAalt2(idx1),'x');
legend("D=0","D=2","D=4",'Location','southeast'); title("Test statistic values"); grid on
xlabel('subjects ordered by TA from D=0'), ylabel('TA');