function [Guu_ensemble,f_ensemble] = section_data(data,rate,num_section)
%% Section the data and ensemble average
%% David W. Fredriksson Fall 2001 - OE810


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sectioning the Data

n = length(data);				
L = n/num_section;		%length of each section


j =1;
k = 1;
for i=1:L:n
   zac_set(j,:)=data(1,i:j*L);
   j = j+1;
end

%size(zac_set);

%for i =1:(num_section+(num_section-1))
%   mean_set = mean(zac_set(i))
%   zac_set(i) = zac_set(i)-mean_set;
%end
%clear i

%zac_set(1,:)
%zac_set = zac_set;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Change over to frequency domain for section and 
%% Ensemble average

for i = 1:(num_section)
   [Guu_set(i,:),fset] = Guu_calc(zac_set(i,:),rate);
end

%Guu_set;
Guu_ensemble = mean(Guu_set);
f_ensemble = fset;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting

%figure(22)
%for i = 1:8
%subplot(211),plot(fset,Guu_set(1,:),'r')
%hold on,plot(fset,Guu_set(2,:),'k.')
%plot(fset,Guu_set(3,:),'b.')
%plot(fset,Guu_set(4,:),'g.')
%plot(fset,Guu_set(5,:),'m.')
%plot(fset,Guu_set(6,:),'y.')
%plot(fset,Guu_set(7,:),'c.')
%xlim([0 2])


%subplot(212)
%%plot(fuu,Guu_acc,'y')
%plot(fset,Guu_ensemble,'b')
%xlim([0 2])
%
%end