function[mean, var, std_dev] = stats(y)
%%Obtaining three main statistics of a given time series
%% N is the number of points y is the time series input

sum = 0;
sq_sum = 0;
N = length(y);

for n = 1:N
    
    sum = sum+y(n);
    sq_sum = sq_sum + y(n).^2;
    
end

%Calculating the mean of the sum
mean = sum/N;

%Calculating the variance
var = sq_sum/N - mean^2;

%Calculate the standard deviation
std_dev = sqrt(var);

end
