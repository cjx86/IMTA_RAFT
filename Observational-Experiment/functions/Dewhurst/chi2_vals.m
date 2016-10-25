%Lookup chi-squared percentage points

function [chi2_l, chi2_h]=chi2_vals(alph,nDOF)
load chi2_vals.mat

%Testing
% alph=0.05;
% nDOF=12;

if nDOF<=max(n)
    
[X,Y]=meshgrid(p,n);

V=interp2(X,Y,c2_val,[1-alph alph], [nDOF nDOF]);
chi2_h=V(1);
chi2_l=V(2);

else
    chi2_h=0;
    chi2_l=0;
%     error('Sorry, that"s too much freedom!')
%    chi2_h=nDOF*(1-2/(9*nDOF)+alph*sqrt(2/(9*nDOF)))^2;
%    chi2_l=nDOF*(1-2/(9*nDOF)+(1-alph)*sqrt(2/(9*nDOF)))^2;
    
end