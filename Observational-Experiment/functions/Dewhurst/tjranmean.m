%Computes running mean of vector with NaNs

function [ybanded]=tjranmean(y,bands)

Nband=round(bands);
if bitget(Nband,1)==0 %bands is even
    if Nband>=bands %bands was rounded up
        Nband=Nband-1; 
    else
        Nband=Nband+1;
    end
    warning('Input "bands" changed to nearest odd')
end

startex=ceil(Nband/2);
endex=floor(length(y)-Nband/2);

ybanded=zeros(length(y),1);
for nn=startex:endex
    ybanded(nn)=tjnanmean(y(nn-(Nband-1)/2:nn+(Nband-1)/2));
end

ybanded(1:startex-1)=ybanded(startex);
ybanded(endex+1:length(y))=ybanded(endex);
