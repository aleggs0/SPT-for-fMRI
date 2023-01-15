function squarednorm = fastsquarednorm(A1,A2,A3,Bsymbol,C1,C2,C3,Dsymbol)
%FASTNORM frob norm of A1*A2*A3+sym2B(Bsymbol)-C1*C2*C3-sym2B(Dsymbol)
%   the user should supply defaults if not used:
%   A3=C3=1 if 2D
%   Bsymbol=Dsymbol=zeros(length_x,length_y,length_z) if no banded part
%   C1=zeros(size(A1)), C2=zeros(size(A1)), C3=zeros(size(A1)) if no truth
[length_x, length_y, length_z]=size(Bsymbol);

Rs1=zeros(3,3,length_x);
for i=1:length_x
    diagB1i=ones(length_x+1-i,1);
    diagA1i=diag(A1,i-1);
    diagC1i=diag(C1,i-1);
    M=[diagB1i, diagA1i, diagC1i];
    Recon=qr(M,"econ");
    Rs1(1:size(Recon,1),:,i)=Recon;
end
Rs2=zeros(3,3,length_y);
for i=1:length_y
    diagB2i=ones(length_y+1-i,1);
    diagA2i=diag(A2,i-1);
    diagC2i=diag(C2,i-1);
    M=[diagB2i, diagA2i, diagC2i];
    Recon=qr(M,"econ");
    Rs2(1:size(Recon,1),:,i)=Recon;
end
Rs3=zeros(3,3,length_z);
for i=1:length_z
    diagB3i=ones(length_z+1-i,1);
    diagA3i=diag(A3,i-1);
    diagC3i=diag(C3,i-1);
    M=[diagB3i, diagA3i, diagC3i];
    Recon=qr(M,"econ");
    Rs3(1:size(Recon,1),:,i)=Recon;
end

squarednorm=0;
for i=1:length_x
    for j=1:length_y
        for k=1:length_z
            w=[Bsymbol(i,j,k)-Dsymbol(i,j,k);1;-1]; w=permute(w,[4 3 2 1]);
            array=sum(permute(Rs1(:,:,i),[1 3 4 2]).*w.*permute(Rs2(:,:,j),[3 1 4 2]).*permute(Rs3(:,:,k),[3 4 1 2]),4);
            %must collaspe the fourth index!
            squarednorm=squarednorm+norm(array,'fro')^2*...
                min(i,2)*min(j,2)*min(k,2);
            %the non-main-diagonals contribute double, as the
            %upper triangle is duplicated in the lower triangle
        end
    end
end
end