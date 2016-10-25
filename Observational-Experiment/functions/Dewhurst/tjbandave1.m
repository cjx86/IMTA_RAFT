%return band averaged series

function [xbanded ybanded]=tjbandave1(x,y,band)

Nband=round(band);
kk=1;
for nn=1:floor((length(y)-1)/band)

    xbanded(nn)=tjnanmean(x(kk:kk+Nband-1));
    ybanded(nn)=tjnanmean(y(kk:kk+Nband-1));
    kk=kk+Nband;

end