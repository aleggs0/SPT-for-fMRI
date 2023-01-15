function [reweighted_symbol,naive_symbol] = fasttoeplitzsymbol(Xtt)
%FASTTOEPLITZSYMBOL Toeplitz symbol of each frame
[length_x,length_y,length_z]=size(Xtt);

symbol=ifftn(abs(fftn(Xtt,[length_x*2-1,length_y*2-1,length_z*2-1])).^2)/(length_x*length_y*length_z);
naive_symbol=symbol(1:length_x,1:length_y,1:length_z);
reweighted_symbol=zeros(length_x,length_y,length_z);
for k=0:length_z-1
    for i=0:length_x-1
        for j=0:length_y-1
            reweighted_symbol(i+1,j+1,k+1)=naive_symbol(i+1,j+1,k+1)...
                /(length_x-i)/(length_y-j)/(length_z-k)*(length_x*length_y*length_z);
        end
    end
end
end