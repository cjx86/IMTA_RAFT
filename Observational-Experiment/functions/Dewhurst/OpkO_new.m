function [pks negpks pklocs pklocsneg Osdex Os]=OpkO_new(x,y) 

%Find positive and negative peaks between consecutive zero crossings.
%By Toby Dewhurst, 2013.
%Calls 'crossed.m'.

% pks: Value of upper peak
% negpks: Value of lower peak
% pkloc: Index of upper peak
% pklocneg: Index of lower peak
% Osdex: Indicies of zero crossings
% Os: x-values at zero crossings

% Example:  (Run from command line)
% X=(1:0.05:6*pi)';Y=-sin(X).*rand(length(X),1);
% [pks negpks pklocs pklocsneg Osdex Os]=OpkO(X,Y);
% amp=mean(pks-negpks)/2
% figure;
% plot(X,Y);
% hold all
% plot(X(Osdex),zeros(length(Osdex),1),'+');
% plot(X(pklocs),pks,'*');plot(X(pklocsneg),negpks,'*')
% ylim=[-1.2 1.2];
% legend('Data','Zero crossings','Positive peaks','Negative peaks')

%Put x, y, in columns.
[m, n]=size(x);
if m<n
    x=x';
end
[m, n]=size(y);
if m<n
    y=y';
end

[Osdex, Os]=crossed(x,y); %Find zero crossings [index value] %Returns index after crossing)

poskk=1;
negkk=1;

for kk=1:length(Osdex)-1
    
    if (Osdex(kk+1)-Osdex(kk))==1 %If consecutive crossings
        if y(Osdex(kk))>0 %I.e. if looking for positive peak
            pks(poskk,1)=y(Osdex(kk));
            pklocs(poskk,1)=0; %Osdex(kk) will be added later;
            poskk=poskk+1;
        else
            negpks(negkk,1)=y(Osdex(kk));
            pklocsneg(negkk,1)=Osdex(kk);
            negkk=negkk+1;
        end
    else %Typical
        if y(Osdex(kk)+1)>0 %I.e. if looking for positive peak
            [pks(poskk,1),pklocs(poskk,1)]=max(y(Osdex(kk):Osdex(kk+1)));
            poskk=poskk+1;
        else
            [negpks(negkk,1),pklocsneg(negkk,1)]=min(y(Osdex(kk):Osdex(kk+1)));
            negkk=negkk+1;
        end
    end
end

%Adjust for the time series being cut
% if y(Osdex(1)+1)>0; %If first peak is positive
%     pklocs=pklocs+Osdex(1:2:end-1)-1;
%     pklocsneg=pklocsneg+Osdex(2:2:end-1)-1;
% else 
%     pklocs=pklocs+Osdex(2:2:end-1)-1;
%     pklocsneg=pklocsneg+Osdex(1:2:end-1)-1;
% end

then=1
