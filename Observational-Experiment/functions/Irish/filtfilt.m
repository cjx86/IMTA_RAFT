function y = filtfilt(b,a,x)
%FILTFILT Zero-phase forward and reverse digital filtering.
%	Y = FILTFILT(B, A, X) filters the data in vector X with the
%	filter described by vectors A and B to create the filtered
%	data Y.  The filter is described by the difference equation:
%
%	  y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%	                   - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%
%	After filtering in the forward direction, the filtered
%	sequence is then reversed and run back through the filter.
%	The resulting sequence has precisely zero-phase distortion
%	and double the filter order.  Care is taken to minimize
%	startup and ending transients by matching initial conditions.
%	See also FILTER.

%	L. Shure 5-17-88
%	(c) Copyright 1988, by The MathWorks, Inc.

[m,n] = size(x);
nb = length(b);
na = length(a);
nfilt = max(nb,na);
% set up initial conditions to remove dc offset problems at the beginning of
% the sequence
if nb == na
	zi = cumsum(b(nfilt:-1:1)-a(nfilt:-1:1)); zi = zi(nfilt-1:-1:1);
else
	if na < nfilt
		bb = b;
		aa = a;
		aa(nfilt) = 0;
	else	% nb < nfilt
		aa = a;
		bb = b;
		bb(nfilt) = 0;
	end
	zi = cumsum(bb(nfilt:-1:1)-aa(nfilt:-1:1)); zi = zi(nfilt-1:-1:1);
	clear bb aa 
end
% reflect beginning sequence of data so slopes match and filter this
% precursory signal so end effects are reduced
nfact = 3*(nfilt-1);
if m == 1	% row data
	y = [2*x(1)-x((nfact+1):-1:2) x];
else 	% column data
	y = [2*x(1)-x((nfact+1):-1:2);x];
end
y = filter(b,a,y,[zi*y(1)]);
% reverse data, prepend and filter again
y(:) = y(length(y):-1:1);
if m == 1	% row data
	y = [2*y(1)-y((nfact+1):-1:2) y];
else 	% column temp
	y = [2*y(1)-y((nfact+1):-1:2);y];
end
y = filter(b,a,y,[zi*y(1)]);
% remove pieces of y that were tacked on to eliminate end effects
y([1:nfact max(m,n)+nfact+(1:nfact)]) = [];
% reverse the final sequence to be in original direction
y = y(length(y):-1:1);
