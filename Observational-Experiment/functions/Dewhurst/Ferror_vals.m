%Copied from myFerror, a temporaray F-distributed error bar calculator and plotter

function [Flow,Fhigh]=Ferror_vals(alph,df1,df2)
load('F_dist_vals.mat')
% %Trim values to eliminate Inf values.
% df1s=df1s(1:end-1);
% df2s=df2s(1:end-1);
% Fs=Fs(1:end-1,1:end-1,:);
%Replace Inf with large number:
df1s(end)=10E9;
df2s(end)=10E9;


akk=find(alphas==alph,1);


if isempty(akk)
    error('Number of DOFs not handled')
end

[X,Y]=meshgrid(df1s,df2s);

Fhigh=interp2(X,Y,Fs(:,:,akk),df1,df2);
Flow=1/interp2(X,Y,Fs(:,:,akk),df2,df1);


% %Application for RAOs should look like:
% Lerr=(Y/sqrt(Fhigh));
% Uerr=(Y/sqrt(Flow));
% %Because RAO=sqrt(Sj/Sw)

