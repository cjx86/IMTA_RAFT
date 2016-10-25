function[tsf, y_mean, y_var, y_std, N_y, y_sq_mean] = stddev_filter(y,num_std,pass)
%%Obtaining three main statistics of a given time series
%% N is the number of points y is the time series input

[y_mean, ~, y_std] = stats(y); 

for ii = 1:pass+1
    
    %m = 0;
    N_y = 0;
    sum = 0;
    sq_sum = 0;
    
    for n = 1:length(y)
    
        if y(n) > y_mean+num_std*y_std || y(n) < y_mean-num_std*y_std
        
            y(n) = nan;
                %If you want to know which pass and index that is cut
            %disp(['Pass ',num2str(ii),' Index ',num2str(n)])
            
            tsf(n) = NaN;
            
        end
        
        if isnan(y(n)) ~= 1
            
            sum = sum+y(n);
            sq_sum = sq_sum + y(n)^2;
            N_y = N_y+1;
            %m = m+1;
            tsf(n) = y(n);
            
        end
    
    end
    
    %N_y = length(tsf);
    
    %Calculating the mean of the sum
    y_mean = sum/N_y;
    y_sq_mean = sq_sum/N_y;
     
    %Calculating the variance
    y_var = sq_sum/N_y - y_mean^2;

    %Calculate the standard deviation
    y_std = sqrt(y_var);
    
    clear y
    y = tsf;

    
end

end