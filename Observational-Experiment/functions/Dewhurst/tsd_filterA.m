%# Standard Deviviation Filter for matrices
function [A]=tsd_filterA(A,Nstds,passes,dim,sub)
[m,n]=size(A);
for pass=1:passes+1 %Pass once for first stats
    %Find sum of all data points. (And set up for variance.)
        A(bsxfun(@gt,A,Nstds*nanstd(A,0,dim)))=sub; %Replace values outside threshold
end