function [L,Lo] = find_L_dispersion(T,h)
% function num = f(k,sig,h)
% routine to determine the wave number
% this routine calls the function fminbnd 
% T is the wave period in seconds
% h is the water depth in m
 
pi=3.14159;
g=9.81;
om=2*pi./T;
Lo=g*T.^2/(2*pi); 
ko=2*pi./Lo;

k=fminbnd(@(k) dispers_func(k,om,h),ko/2,100*ko);
L=2*pi./k;

  
return

function num = dispers_func(k,om,h)
% function num = f(k,sig,h)
% routine to determine the dispersion relationship
% this routine is called by the fminbnd function
% k is the wave numer that will be determined
% sig is wave frequency (=2*pi/T)
% h is the water depth
 
g=9.81;
 
num =abs(om^2-g*k*tanh(k*h));
return