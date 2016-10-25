% Get amplitude and phase information from a time-series
% Inputs
% t - time
% y - independent var of same length as t
% 
function [y_amp]=amps_CS(t,y,good_stdex,good_endex,T)

dt=t(2)-t(1); % Time step;
[y_pks, y_negpks, y_pklocs, ~, ~, ~]= ... 
    OpkO_new(t(good_stdex:good_endex),y(good_stdex:good_endex)-mean(y(good_stdex:good_endex)));
%OpkO_new(t(good_stdex:good_endex),y(good_stdex:good_endex)-tjranmean(y(good_stdex:good_endex),round(T/dt)))
y_amp=mean(y_pks(1:min(length(y_pks),length(y_negpks)))-y_negpks(1:min(length(y_pks),length(y_negpks))))/2; %Truncte longer series
%y_phase=mean(t(good_stdex+y_pklocs)-(t(good_stdex)+T*(0:length(y_pklocs)-1)'))*2*pi/T;
      
end