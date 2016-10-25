      function [u,v,ma] = accel2motion20(a, f, v0, u0)
% program originally written by Jason Gobat 
% messed up by Jim Irish in Jan 2001 to look at UNH wave data   
% The 20 indicates value used for cutoff freqency in the butterworth filter
%   inputs:
%      a = acceleration time series input
%      f = sample frequency = 10 Hz for waverider and motion package
%      v0 = integration constant in velocity, usually zero
%      u0 = integration constant in displacement, usually zero
%   outputs:
%      u = displacement
%      v = velocity
%
%   remove mean and trend from data
   ma=mean(a);
   a=a-ma;
	disp(['Mean removed = ',num2str(ma),' g''s'])
%    clear ma
   a = detrend(a);
% convert from g's to m/s²
   a=9.8.*a;
   dt = 1.0/f;
   disp('Creating high pass filter for displacement')
   [bf,af] = butter(4,(1/20)/(f/2),'high');
   n = length(a);
   v = zeros(n,1);
   u = zeros(n,1);
% starting velocity = v0   
   disp('Integrating to velocity')
   v(1) = v0;
% trapizoidal integration of acceleration to velocity  
   for i = 2:n
      v(i) = dt*(a(i) + a(i-1))/2 + v(i-1);
   end
% detrend and filter to remove low frequencies   
%   v = filtfilt(bf,af,detrend(v));
   v = filtfilt(bf,af,v);
% starting position - u0
   disp('Integrating to dispacement')
   u(1) = u0;
% trapizoidal integration of velocity to position
   for i = 2:n
      u(i) = dt*(v(i) + v(i-1))/2 + u(i-1);
   end
% apply Butterworth filter again
%    u = filtfilt(bf,af,detrend(u));
    u = filtfilt(bf,af,u);

