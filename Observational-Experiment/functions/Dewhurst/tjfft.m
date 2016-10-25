%FFT:
%Returns normalized, symmetrical Discrete Fourier Transform, using fft.m
%padyes=1 for zero padding.
function [Gfast freqs]=tjfft(gn,delta,padyes)
%gn=gn-tjnanmean(gn);
N=length(gn);

%padding

if padyes==0
    %Check even or odd
    if bitget(N,1)==1  %(TS is odd length)
        gnswap=[gn((N-1)/2+1:N);gn(1:(N-1)/2)];
        Gfastswap=fft(gnswap)/N;
        Gfast=[Gfastswap((N-1)/2+1:N);Gfastswap(1:(N-1)/2)];
        freqs=[(-1/2/delta-1/N/2/delta):1/N/delta:(1/2/delta-3/N/2/delta)];
    elseif bitget(N,1)==0  %(TS is even length)
        gnswap=[gn(N/2+1:N);gn(1:N/2)];
        Gfastswap=fft(gnswap)/N;
        Gfast=[Gfastswap((N)/2+1:N);Gfastswap(1:(N)/2)];
        freqs=[-1/2/delta+1/N/delta:1/N/delta:(1/2/delta)];
    end
elseif padyes==1
    NFFT = 2^nextpow2(N); % Next power of 2 from length of y
    padder=zeros(floor((NFFT-length(gn))/2),1);
    gnpad=[0;padder;gn;padder];
    N=NFFT;
    gnswap=[gnpad(N/2+1:N);gnpad(1:N/2)];
    Gfastswap=fft(gnswap)/N;
    Gfast=[Gfastswap((N)/2+1:N);Gfastswap(1:(N)/2)];
    freqs=[-1/2/delta+1/N/delta:1/N/delta:(1/2/delta)];
end