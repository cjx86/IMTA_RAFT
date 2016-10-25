      function [a, mr] = rate2angle20(r, f, a0)
% program originally written by Jason Gobat 
% messed up by Jim Irish in Jan 2001 to look at UNH wave data   
% The 20 indicates value used for cutoff freqency in the butterworth filter
%   inputs:
%      r = rate time series input - degrees per second
%      f = sample frequency = 10 Hz for motion package
%      a0 = integration constant in angle, usually zero
%   outputs:
%      a = angle
%
%   remove mean and trend from data
   mr=mean(r);
   r=r-mr;
   disp(['Mean rate removed = ',num2str(mr),' deg/s'])
%    clear mr
   r = detrend(r);
%
   dt = 1.0/f;
   disp('Creating high pass filter for angle')
   [bf,af] = butter(4,(1/20)/(f/2),'high');
   n = length(r);
   a = zeros(n,1);
%
% detrend and filter to remove low frequencies   
%   r = filtfilt(bf,af,detrend(r));
    r = filtfilt(bf,af,r);
% starting angle - a0
   disp('Integrating rate of rotation to angle')
   a(1) = a0;
% trapizoidal integration of velocity to position
   for i = 2:n
      a(i) = dt*(r(i) + r(i-1))/2 + r(i-1);
   end
% apply high pass Butterworth filter 
%    u = filtfilt(bf,af,detrend(u));
     a = filtfilt(bf,af,a);

