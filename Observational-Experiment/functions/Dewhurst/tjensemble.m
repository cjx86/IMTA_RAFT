%Ensemble average (for spectra)

function [Sj2side ensSj1side freqs]=tjensemble(gn,delta,padyes,ensembles)

%Break-up, get Spectral Densities
Nens=floor(length(gn)/ensembles);
kk=1;
for nn=1:ensembles
    gnpart=gn(kk:kk+Nens-1);
    %Get FFT
    [Gfast(:,nn) freqs(:,nn)]=tjfft(gnpart,delta,padyes);
    Cj(:,nn)=sqrt((real(Gfast(:,nn))).^2+(imag(Gfast(:,nn))).^2);
    Sj2side(:,nn)=Nens*delta*Cj(:,nn).^2;
    kk=kk+Nens;
end

%Average ensembles
for nn=1:length(Sj2side)
    ensSj2side(1,nn)=tjnanmean(Sj2side(nn,:));
end

ensSj1side=2*ensSj2side(ceil(length(ensSj2side)/2)+1:length(ensSj2side));