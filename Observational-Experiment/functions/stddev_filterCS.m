function[y, tsf, y_mean, y_std, N_y, y_sq_mean] = stddev_filterCS(t,y,num_std,k)
%%Obtaining three main statistics of a given time series
%% N is the number of points y is the time series input
%% t is the time interval
%% num_std is the number of standard deviations outside of which to eliminate
%% k is the number of times to repeat the filter 


[y_mean, vari, y_std] = stats(y); 


for ii = 1:k
    
    m = 0;
    N_y = 0;
    sum = 0;
    sq_sum = 0;
    
    for n = 1:length(y)
    
        if y(n) > y_mean+num_std*y_std | y(n) < y_mean-num_std*y_std
        
            y(n) = nan;
        
        end
        
        if isnan(y(n)) ~= 1
            
            sum = sum+y(n);
            sq_sum = sq_sum + y(n)^2;
            N_y = N_y+1;
            m = m+1;
            tsf(m,2) = y(n);
            tsf(m,1) = t(n);
        
        end
    
    end
    
    N_y = length(tsf);
    
    %Calculating the mean of the sum
    y_mean = sum/N_y;
    y_sq_mean = sq_sum/N_y;
     
    %Calculating the variance
    varr = sq_sum/N_y - y_mean^2;

    %Calculate the standard deviation
    y_std = sqrt(varr);
    
end

end