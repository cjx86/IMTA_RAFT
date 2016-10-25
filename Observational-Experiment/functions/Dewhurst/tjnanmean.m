%Computes mean of vector with NaNs
function [mean]=tjnanmean(ts)

summer=0;
kk=0;
for n=1:length(ts)
    if isnan(ts(n))==0
        summer=ts(n)+summer;
    kk=kk+1;
    end
end
mean=summer/kk;