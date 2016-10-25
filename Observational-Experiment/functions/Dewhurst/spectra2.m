%Calculates and returns 1-sided power spectral density of a random
%time series. Uses a hanning data window.
%Calls 

function [bandfreqs Sj1side conf bandwidth Sjarea gnvar]=spectra2(gn,freq,ensembles,bands,padyes,hannme,stds,passes,mean_bands,alph)

%gn %Time series
%freq %Sampling frequency
%ensembles %Number of ensembles to average
%bands %Number of bands over which to average
%padyes %Set to 1 to pad time series to take FFT. Recommended: padyes=0.
%hannme %Set to 1 to apply hanning window. Set to 0 otherwise.
%stds %Standard deviations to filter--high means don't filter
%passes %Number of times to filter data. 
%mean_bands. Bands over which to average when removing running mean. To
%remove single mean, set mean_bands=0

if mean_bands>=1
gn=gn-tjranmean(gn,mean_bands); %Remove running mean. May be unnecessary because mean of each ensemble removed
else
    gn=gn-tjnanmean(gn);
end

%Filter, get variance, etc.
[gn, gnmean, gnvar, ~, gnnumgood,~] =  stddev_filter(gn,stds,passes)

%Put in column vector for hanning.m,
[m n]=size(gn);
if n>m
    gn=gn';
end

delta=1/freq; %s
N=gnnumgood;
% %Ensemble
if hannme==1
    [Sj2side, Sj1side, freqs]=tjensemblehan2(gn,delta,padyes,ensembles);
%     %Boost
%     F=gnvar/gnvarhan;
%     Sj1side=Sj1side*F;
% tjensemblehan2 means that hanning happens to each ensemble. 
else
    [Sj2side, Sj1side, freqs]=tjensemble(gn,delta,padyes,ensembles);
end
%Band Average
% [bandfreqs Sj1side]=tjbandave1(freqs(ceil(length(freqs)/2:length(freqs))),Sj1side,bands);
[bandfreqs, Sj1side]=tjbandave1(freqs(ceil(length(freqs)/2)+1:end),Sj1side(1:end-1),bands);

%Calculate confidence
DOFs=2*bands*ensembles;
[chi2l, chi2h]=chi2_vals(alph,DOFs);

Sl=DOFs*Sj1side/chi2l; %Random S^ value
Sh=DOFs*Sj1side/chi2h;
conf=(Sh-Sl)'; %95 percent confidence interval

%Find variance by area
bandwidth=1/N/delta*bands*(ensembles+1/bands);

Sjarea=trapz(bandfreqs,Sj1side);
end
% [Sjarea gnvar]; %If analysis is correct, the area under the power spectrum should equal the variance of the original time series. 
