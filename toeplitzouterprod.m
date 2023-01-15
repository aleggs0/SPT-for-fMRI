function symbol = toeplitzouterprod(C1,C2,C3)
%TOEPLITZOUTERPROD Toeplitz symbol for the outer product of C1,C2,C3
%   Works for 3D. For 2D set A3=1.
c1_averages=zeros(length(C1),1);
for h1=0:length(C1)-1
    c1_averages(h1+1)=mean(diag(C1,h1));
end
c2_averages=zeros(1,length(C2));
for h2=0:length(C2)-1
    c2_averages(h2+1)=mean(diag(C2,h2));
end
c3_averages=zeros(1,1,length(C3));
for h3=0:length(C3)-1
    c3_averages(h3+1)=mean(diag(C3,h3));
end
symbol=c1_averages.*c2_averages.*c3_averages;
end