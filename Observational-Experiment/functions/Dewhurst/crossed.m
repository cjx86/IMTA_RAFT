%Gives interpolated values of zero crossings in the form Os=[indices values]
%By Toby Dewhurst, 2013.
%Warning: Will not detect zero crossing in last two y values

%Returns:
%'Osdex'=index of zero crossing
%'Os'=x-value of zero crossing (from linear interpolation).
function [Osdex, Os]=crossed(x,y)

if min(y)>0
    error('No negative values found')
end

%Put x, y, in columns.
[m, n]=size(x);
if m<n
    x=x';
end
[m, n]=size(y);
if m<n
    y=y';
end

%create a shifted series
y2=[0;y];
y=[y;0];

signs=y2.*y; 

jj=1;
for kk=1:length(x)-1
    if signs(kk)<0
        Osdex(jj,1)=kk;
        Os(jj,1)=x(kk-1)-y(kk-1)/(y(kk)-y(kk-1))*(x(kk)-x(kk-1));
        jj=jj+1;
    elseif signs(kk)==0 && signs(kk+1)==0 %True iff y(kk)==0
        Osdex(jj,1)=kk;
        Os(jj,1)=y(kk);
        jj=jj+1;
    end
    %hold all;plot(x(kk),y(kk),'*') %I forget why I put this in here
end
