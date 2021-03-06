function [num, den, z, p] = butter(n, Wn, ftype)
%BUTTER	Butterworth digital filter design.
%	[B,A] = BUTTER(N,Wn) designs an N'th order lowpass digital
%	Butterworth filter and returns the filter coefficients in length
%	N+1 vectors B and A.  The cut-off frequency Wn must be
%	0.0 < Wn < 1.0, with 1.0 corresponding to half the sample rate.
%
%	If Wn is a two-element vector, Wn = [W1 W2], BUTTER returns an 
%	order 2N bandpass filter with passband  W1 < W < W2.
%	[B,A] = BUTTER(N,Wn,'high') designs a highpass filter.
%	[B,A] = BUTTER(N,Wn,'stop') is a bandstop filter if Wn = [W1 W2].
%	
%	When used with three left-hand arguments, as in
%	[Z,P,K] = BUTTER(...), the zeros and poles are returned in
%	length N column vectors Z and P, and the gain in scalar K. 
%
%	When used with four left-hand arguments, as in
%	[A,B,C,D] = BUTTER(...), state-space matrices are returned.
%
%	See also BUTTORD, CHEBY1, CHEBY2, FREQZ and FILTER.

%	J.N. Little 1-14-87
%	Revised 1-14-88 JNL, 4-29-88 LS
%	(c) Copyright 1987-88, by The MathWorks, Inc.

%	References:
%	  [1] T. W. Parks and C. S. Burrus, Digital Filter Design,
%	      John Wiley & Sons, 1987, chapter 7, section 7.3.3.

btype = 1;
if nargin == 3
	btype = 3;
end
if max(size(Wn)) == 2
	btype = btype + 1;
end

% step 1: get analog, pre-warped frequencies
fs = 2;
u = 2*fs*tan(pi*Wn/fs);

% step 2: convert to low-pass prototype estimate
if btype == 1	% lowpass
	Wn = u;
elseif btype == 2	% bandpass
	Bw = u(2) - u(1);
	Wn = sqrt(u(1)*u(2));	% center frequency
elseif btype == 3	% highpass
	Wn = u;
elseif btype == 4	% bandstop
	Bw = u(2) - u(1);
	Wn = sqrt(u(1)*u(2));	% center frequency
end

% step 3: Get N-th order Butterworth analog lowpass prototype
[z,p,k] = buttap(n);

% Transform to state-space
[a,b,c,d] = zp2ss(z,p,k);

% step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn
if btype == 1		% Lowpass
	[a,b,c,d] = lp2lp(a,b,c,d,Wn);

elseif btype == 2	% Bandpass
	[a,b,c,d] = lp2bp(a,b,c,d,Wn,Bw);

elseif btype == 3	% Highpass
	[a,b,c,d] = lp2hp(a,b,c,d,Wn);

elseif btype == 4	% Bandstop
	[a,b,c,d] = lp2bs(a,b,c,d,Wn,Bw);
end

% step 5: Use Bilinear transformation to find discrete equivalent:
[a,b,c,d] = bilinear(a,b,c,d,fs);

if nargout == 4
	num = a;
	den = b;
	z = c;
	p = d;
else	% nargout <= 3
% Transform to zero-pole-gain and polynomial forms:
	if nargout == 3
		[z,p,k] = ss2zp(a,b,c,d,1);
		num = z;
		den = p;
		z = k;
	else % nargout <= 2
		den = poly(a);
		num = poly(a-b*c)+(d-1)*den;
	end
end

