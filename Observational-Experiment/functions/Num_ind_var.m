function[N_st, Sa] =  Num_ind_var(minlag,maxlag,x,y)
%% Where minlag is the minimum lag in the x variable 
%% maxlag is the maximum lag in the x variable
%% x and y are the desired time series
%% del is the step in the x variable


%%

%Call Cross Correlation for the positive and negative long lags
corrxy = cross_corr(minlag,maxlag,x,y);
corryx = cross_corr(minlag,maxlag,y,x);

corrxy_sq = (corrxy).^2;
corryx_sq = fliplr(corryx).^2;

%Initialze variable for the sum of the squares of the x-corrs
Sa_p = 0;
Sa_n = 0;

M = (maxlag-minlag)+1;


for jj = 1:M
    
    Sa_p = Sa_p + corrxy_sq(jj);
    
    Sa_n = Sa_n + corryx_sq(jj);
    
end

% ii;
% Sa_p
% Sa_ptest = sum(corrxy_sq)
% Sa_n
% Sa_ntest = sum(corryx_sq)

Sa = (Sa_p+Sa_n)/2;

N_st = M/Sa;

end