%IMU_meta_CS
% y
% y
% y
%  x x x
% Corey Sullivan, 2016
    % based on code from Toby Dewhurst, 2015.
clear;clc;close all;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))
addpath(genpath('../'));
addpath(genpath('../functions'));

% %Load filenames
%  dep='First';
% pathdir=['..\Data\' dep '\'];
% 
% dirData = dir(pathdir);      %# Get the data for the current directory
% dirIndex = [dirData.isdir];  %# Find the index for directories
% fileList = {dirData(~dirIndex).name}';  %# Get a list of the files
% 
% if strcmp(dep,'First')
% 
%     st_id=3; %First file of interest
%     ids=st_id:length(fileList);
%     
% end
% 
% %Pre-allocate
% Disp_std=zeros(length(ids),3);
% Ang_std=zeros(length(ids),3);
% a_means=zeros(length(ids),3);
% rot_means=zeros(length(ids),3);
% batts=zeros(length(ids),1);
% 
% Sj1side20=zeros(length(ids),65,6); %Will vary with ensembles in mpnorm_fun
% Sjarea20=zeros(length(ids),1,6);
% tsvar20=zeros(length(ids),6);
% 
% %Spectrating options
% padyes=0;
% hannme=1; %Variance error comes from hanning=1...?; Fixed.
% ensembles=1;
% alph=0.05; %Confidence margin. I.e. alph=0.10 for 90% confidence. (See Bendat and Piersol, p. 90).
% bands=1;
% stds=100; %Standard deviations to filter--high # means don't filter
% passes=1;
% %Set number of indices over which to find a running mean. To just demean: mean_bands=0
% mean_bands(1)=0; 


%  for id=ids;
%      kk=id+1-st_id;
%      [date, time, con95, Disp, Ang, Disp_nf, Ang_nf, t_plot, batt, a_mean,a_mean_nf, rot_mean, bandfreqs, Sj1side20(kk,:,:), Sj1side30(kk,:,:), Sj1side_nf(kk,:,:), bands_IMU, bandwidth, Sjarea20(kk,:,:), tsvar20(kk,:), Sjarea30(kk,:,:), tsvar30(kk,:), Sjarea_nf(kk,:,:), tsvar_nf(kk,:)]=mpnorm_fun_CS([pathdir '\' fileList{id}],ensembles,padyes,hannme,stds,passes,mean_bands,alph);
%          %Get information from raw data
%          
%      ts_plot(kk,:) = t_plot;    %time vector to plot the data
%      Displ(kk,:,:) = Disp;   %Displacement 3, 2D matrices x,y,z respectively
%      Angl(kk,:,:) = Ang;  %Rotation 3, 2D about: x,y,z respectively
%      Displ_nf(kk,:,:) = Disp_nf;  
%      Angl_nf(kk,:,:) = Ang_nf; 
%      con95(:,:,kk)=con95;
%      %bfreq(kk,:) = bandfreqs;
%          
%     
%      Disp_std(kk,:)=std(Disp,1);     %Std dev of each displacement
%      Ang_std(kk,:)=std(Ang,1);   %Std dev of each rotation
%      clear('Disp','Ang');
%      
%      
%      dates(kk,:)=date;   %Fill a matrix of each date for t.s. plotting
%      times(kk,:)=time;   %Fill a matrix of times for t.s. plotting
%      a_means(kk,:)=a_mean; %Array of mean accelerations
%      a_means_nf(kk,:) = a_mean_nf;
%      rot_means(kk,:)=rot_mean;   %Array of mean rotations
%      batts(kk)=str2num(batt);    %Battery Output
%  
%  end
%  save('data_read_from_IMU.mat')

load('data_read_from_IMU.mat')

close all

    %Produce the T.S. from each 18  min sample
ts_plot = ts_plot(1,:);   

displ_x = Displ([1:28 30:end],:,1);
displ_y = Displ([1:28 30:end],:,2);
displ_z = Displ([1:28 30:end],:,3);

%Mean angle
mean_roll=180/pi*(pi/2+atan(a_means([1:28 30:end],3)./a_means([1:28 30:end],1))); %Mean roll angle, deg
mean_roll(mean_roll>90)=mean_roll(mean_roll>90)-180;

mean_roll_nf = 180/pi*(pi/2+atan(a_means_nf([1:28 30:end],3)./a_means_nf([1:28 30:end],1)));
mean_roll_nf(mean_roll_nf>90)=mean_roll_nf(mean_roll_nf>90)-180;


mean_pitch=180/pi*(pi/2+atan(a_means([1:28 30:end],3)./a_means([1:28 30:end],2))); %Mean pitch angle, deg
mean_pitch(mean_pitch>90)=mean_pitch(mean_pitch>90)-180;

mean_pitch_nf=180/pi*(pi/2+atan(a_means_nf([1:28 30:end],3)./a_means_nf([1:28 30:end],2)));
mean_pitch_nf(mean_pitch_nf>90)=mean_pitch_nf(mean_pitch_nf>90)-180;

%% Correct z-displ to COM
    
    %Distances IMU was from COM
L_x = 76*.0254; %Distance IMU is from COM in short direction
L_y = 40*.0254; %Distance IMU is from COM in long direction

    %Build a Low Pass Filter for the Pitch Spectra for Correcting Data
[bfa,afa] = butter(1,2/(10/2)); 

    %Accounting for displ at COM via pitch and roll angles with unfiltered
zdisp_nf = Displ_nf([1:28 30:end],:,3) + L_x.*sind(Angl_nf([1:28 30:end],:,2)) - L_y.*sind(Angl_nf([1:28 30:end],:,1));

IMU_freq = 10;

Ffreq=IMU_freq/(floor(length(displ_z)/ensembles)); %Fourier frequency spacing

bands2=round(0.5/Ffreq/64);

    %Just correcting back to C.O.M. with filtered angles and displacements
zdisp = displ_z + L_x.*sind(Angl([1:28 30:end],:,2)) - L_y.*sind(Angl([1:28 30:end],:,1));

    %Applying the 20 butterworth filter
[bf,af] = butter(4,(1/20)/(10/2),'high');

for ii = 1:length(ids)-1
    
        %Butter Filtering corrected displ - replacing this with mpnorm_fun
    disp_z_corr(ii,:) = filtfilt(bf,af,zdisp_nf(ii,:));
    
        %Spectra of data corrected then filtered
    [bfreq, Sj1side_cor(ii,:), ~, ~, Sjarea_cor(ii,:), tsvar_cor(ii,:)]=spectra2(disp_z_corr(ii,:),IMU_freq,ensembles,bands2,padyes,hannme,stds,passes,mean_bands,alph);
          
        %Spectra, filtered then corrected
    [bfreq2, Sj1side_cor2(ii,:), ~, ~, Sjarea_cor2(ii,:), tsvar_cor2(ii,:)]=spectra2(zdisp(ii,:),IMU_freq,ensembles,bands2,padyes,hannme,stds,passes,mean_bands,alph);
    
       
end

    %Truncate to 0.5 Hz
cut1 = find(bfreq>0.5,1);
bfreqc = bfreq(1:cut1);
Sj1side_cor_c=Sj1side_cor(:,1:cut1);

cut2 = find(bfreq2>0.5,1);
bfreqc2 = bfreq2(1:cut2);
Sj1side_cor_c2=Sj1side_cor2(:,1:cut2);


%% Check on Correction

%Variance
figure,plot(tsvar20([1:28 30:end],3)),hold on,plot(tsvar_cor),plot(tsvar_cor2)
title('Variance')
legend('Original','Corrected then filtered','Filtered then corrected')
hold off

%compare corrected after filtering and corrected before filetering
figure,plot(displ_z(23,:)),hold on,plot(zdisp(23,:),'g'),plot(disp_z_corr(23,:),'r--')
title('Displacement Comparison and corresponding error')
legend('Original','Filtered then corrected','Corrected then filtered')
hold off

%Compare Spectra
figure,semilogy(bandfreqs(2:end),Sj1side20(23,2:end,3),bfreqc(2:end),Sj1side_cor_c(23,2:end))
title('Spectra of Heave Motion from IMU')
xlabel('Frequency (Hz)')
ylabel('Spectral Energy Density, m^2/Hz')
legend('Original Z-Displ Spectrum','Corrected to C.O.M. Spectrum')

    %Replace Sj1side20 with the corrected back to COM
Sj_zcor = Sj1side_cor_c;

%% Testing a general Low Pass Filter for z displacement (3/22/16)

[bfa,afa] = butter(1,.2); 

 for ii = 1:length(ids)-1   

    test_z_filt(ii,:) = filtfilt(bfa,afa,disp_z_corr(ii,:));
    [testfreq, Sj1side_cortest(ii,:), ~, ~, Sjarea_cortest(ii,:), tsvar_cortest(ii,:)]=spectra2(test_z_filt(ii,:),IMU_freq,ensembles,bands2,padyes,hannme,stds,passes,mean_bands,alph);
    
 end
 
cut3 = find(testfreq>0.5,1);
bfreqtest = testfreq(1:cut3);
Sj1side_cortest_c=Sj1side_cortest(:,1:cut3);

figure,plot(displ_z(23,:)),hold on,plot(disp_z_corr(23,:),'g'),plot(test_z_filt(23,:),'r--')
title('Displacement correction filter check')
legend('Orig','Noisy','Filtered')
hold off
% 
figure,semilogy(bandfreqs(2:end),nanmean(Sj_zcor(:,2:end),1),bfreqtest(2:end),nanmean(Sj1side_cortest_c(:,2:end),1))
title('Spectra of Heave Motion from IMU')
xlabel('Frequency (Hz)')
ylabel('Spectral Energy Density, m^2/Hz')
legend('Noisy','Filtered')

%% Create overall matrices of data

roll = Angl([1:28 30:end],:,2); %Angle of Short Side
pitch = Angl([1:28 30:end],:,1); %Angle of Long Side
yaw = Angl([1:28 30:end],:,3);

Disp_totz = Displ([1:28 30:end],:,3);
%t_tot_plot = ts_plot(1,([1:28 30:end]));

for ii = 2:length(ts_plot(:,1))
    
    t_tot_plot = [t_tot_plot ts_plot(1,:).*ii];
    Disp_totz = [Disp_totz Displ(ii,:,3)];
    
end

%% Plotting specific days UNCOMMENT TO VIEW

% %Days of interest(mm/dd HH): 01/16 ~08  01/16 ~20
%                       % samples 13-16      25-28
figure,subplot(2,1,1),plot(t_plot/60,displ_z(29,:)),ylim([-0.4 0.4])
ylabel('Vertical Displacement [m]')
title('Displacement from 9 PM sample')
xlim([0 1050/60])
subplot(2,1,2),plot(t_plot/60,displ_z(30,:)),ylim([-0.4 0.4])
xlim([0 1050/60])
ylabel('Vertical Displacement [m]')
title('Displacement from 10 PM sample')
xlabel('Time [mins]')
suptitle('Heave Motion Mid Morning 01/16');

    %Plotting the Pitch
    
figure,subplot(4,1,1),plot(t_plot,pitch(30,:),t_plot,roll(30,:)),ylim([-5 5])
subplot(4,1,2),plot(t_plot,pitch(31,:),t_plot,roll(31,:)),ylim([-5 5])
subplot(4,1,3),plot(t_plot,pitch(32,:),t_plot,roll(32,:))
subplot(4,1,4),plot(t_plot,pitch(33,:),t_plot,roll(33,:))
suptitle('Pitch Angle Morning 01/16');    
legend('Pitch','Roll')

figure,
l = 30;
subplot(3,1,1),plot(t_plot(:)./60,pitch(l,:)),xlim([0 17])
title('Pitch Motion Time Series','FontSize',24)
ylabel('Angle of Rotation [deg]')
set(gca,'FontSize',20)
ylim([-4 4])
subplot(3,1,2),plot(t_plot(:)./60,roll(l,:)),xlim([0 17])
title('Roll Motion Time Series','FontSize',24)
ylabel('Angle of Rotation [deg]')
set(gca,'FontSize',20)
ylim([-4 4])
subplot(3,1,3),plot(t_plot(:)./60,yaw(l,:)),xlim([0 17])
title('Yaw Motion Time Series','FontSize',24)
ylabel('Angle of Rotation [deg]')
xlabel('Time [min]')
set(gca,'FontSize',20)
ylim([-4 4])


%Pitch and roll Spectra
figure,
plot(bandfreqs,nanmean(Sj1side20([1:28 30:end],:,4),1),'--','LineWidth',1.5),hold on
plot(bandfreqs,nanmean(Sj1side20([1:28 30:end],:,5),1),'-','LineWidth',1.5)
title('Averaged Spectra of Pitch and Roll from the IMU','FontSize',22)
xlabel('Frequency [Hz]')
ylabel('Spectral Energy Density [m^2/Hz]')
set(gca,'FontSize',18)
legend(' Pitch',' Roll')

%figure,
%for l = 1:length(pitch(:,1))
    plot(t_plot./60,pitch(40,:),'LineWidth',1.25),ylim([-2 2])
title('Pitch Motion Time Series','FontSize',22)
ylabel('Angle of Rotation [deg]')
xlabel('Time [min]')
set(gca,'FontSize',18)
%xlim([15 45])
%set(gca,'XTick',(15:5:45))
%pause
%end

for id =1:length(disp_z_corr(:,1))
    
[amp_heave(id)] = amps_CS(t_plot(:),disp_z_corr(id,:),1,length(disp_z_corr(id,:)));

end

%% Stats Calculations etc.

%SOMETIMES YOU HAVE TO F9 from here to the end to finish it

%Assemble datenum vector
yrstr=cellstr(dates([1:28 30:end],1:2));
yr=2000+str2double(yrstr);
mostr=cellstr(dates([1:28 30:end],4:5));
dastr=cellstr(dates([1:28 30:end],7:8));

hrstr=(times([1:28 30:end],1:2));
mnstr=(times([1:28 30:end],4:5));
scstr=(times([1:28 30:end],7:8));

mo=str2double(mostr);
da=str2double(dastr);
hr=str2double(cellstr(hrstr));
mn=str2double(cellstr(mnstr));
sc=str2double(cellstr(scstr));

datenums=datenum(yr,mo,da,hr,mn,sc);
chk_dat = datestr(datenums);


xlimit=([datenum(2016,1,15,21,0,0) datenums(end)]);
stkk=1; %Index of first good data point. 15-Jan-2016 21:00 (16:00 PM EST)
endkk=length(datenums)+1;
% savename='..\..\Combined\Data\First\IMU_meta.mat';
   
if length(ids)<=2
   stkk=1;
   endkk=length(ids);
end

dropouts=isnan(datenums); %Check for nans
datenums=inpaint_nans(datenums); %There's no nans this is just a check




%% Plot Surge sWAY TS

figure,
plot(t_plot/60,displ_x(47,:),'LineWidth',1)
xlabel('Time [min]')
ylabel('Displacement [m]')
set(gca,'FontSize',20)
title('Sway - X Axis Motion','FontSize',24)
ylim([-0.153 0.15])
xlim([0 17.25])

figure,
plot(t_plot/60,displ_y(47,:),'LineWidth',1)
xlabel('Time [min]')
ylabel('Displacement [m]')
set(gca,'FontSize',20)
title('Surge - Y Axis Motion','FontSize',24)
ylim([-0.12 0.12])
xlim([0 17.25])

figure,
plot(t_plot/60,disp_z_corr(47,:),'LineWidth',1)
xlabel('Time [min]')
ylabel('Displacement [m]')
set(gca,'FontSize',20)
title('Heave - Z Axis Motion','FontSize',24)
ylim([-0.3 0.3])
xlim([0 17.25])

%%
%Plot STD

figure
plot(datenums,Disp_std([1:28 30:end],:))
ylabel('STD of Displacement, m')
%legend('Sway (In/Out of Channel)', 'Surge( Push/Pull on short side)', 'Heave')
datetick('x','mm/dd HH','keepticks')
xlim(xlimit)
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)

figure
plot(datenums,Ang_std([1:28 30:end],:))
ylabel('STD of angle, deg')
legend('Roll (About y)', 'Pitch (about x)', 'Yaw')
datetick('x','mm/dd HH','keepticks')
xlim(xlimit);
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)

%Mean accelerations
figure
plot(datenums,a_means([1:28 30:end],:))
ylabel('Mean Acceleration, g')
legend('x', 'y', 'z')
datetick('x','mm/dd HH','keepticks')
xlim(xlimit);
%% Trim to align with other instruments

dropouts_IMU=dropouts;%([stkk:28 30:endkk]);
datenums_IMU=datenums;

Disp_std=Disp_std([stkk:28 30:endkk],:);
Ang_std=Ang_std(([stkk:28 30:endkk]),:);

a_means=a_means(([stkk:28 30:endkk]),:);
a_means_nf=a_means_nf(([stkk:28 30:endkk]),:);

%mean_pitch=mean_pitch(([stkk:28 30:endkk]));
%mean_roll=mean_roll(([stkk:28 30:endkk]));

bandfreqs_IMU=bandfreqs;

dump = Sj1side20(([stkk:28 30:endkk]),:,:);
clear Sj1side20
Sj1side20=dump;
clear dump

Sj1side30=Sj1side30(([stkk:28 30:endkk]),:,:);
Sj1side_nf=Sj1side_nf(([stkk:28 30:endkk]),:,:);
bandwidth=bandwidth(1);

dump = Sjarea20(([stkk:28 30:endkk]),:);
clear Sjarea20
Sjarea20 = squeeze(dump);
clear dump
dump = tsvar20(([stkk:28 30:endkk]),:);
tsvar20 = squeeze(dump);
clear dump

 
Sjarea30=Sjarea30(([stkk:28 30:endkk]),:);
tsvar30=tsvar30(([stkk:28 30:endkk]),:);

dump = Sjarea_nf(([stkk:28 30:endkk]),:);
clear Sjarea_nf
Sjarea_nf= dump;
clear dump

dump = tsvar_nf(([stkk:28 30:endkk]),:);
tsvar_nf = dump;
clear dump

figure,semilogy(bfreqc(2:end),nanmean(Sj1side_cor_c(:,2:end)),'LineWidth',1.5)
xlabel('Frequency [Hz]')
ylabel('Spectral Energy Density, [m^2/Hz]')
set(gca,'XLim',[0.05 0.5])
%set(gca,'YLim',[10^-10 10^-1])
set(gca, 'YTick', [10^-10 10.^-9 10^-8 10.^-7 10^-6 10.^-5 10^-4 10.^-3 10^-2 10.^-1])
%set(gca, 'YTick', [10^-10 10^-8 10^-6 10^-4 10^-2 10^0])
set(gca,'FontSize',22)
title('Averaged Spectra of Heave Motion from IMU','FontSize',28)

%% % IMU Error calcs. (Commented out because you can do this later)
alph=0.05;
DOFs=2*bands_IMU*ensembles;
[chi2l, chi2h]=chi2_vals(alph,DOFs);
Sl=DOFs/chi2l*Sj1side20; %Needs to be a function of Sj magnitude if not on log plot
Sh=DOFs/chi2h*Sj1side20; %Ends of error bars
errL=Sj1side20-Sl; %Length of error bars
errH=Sh-Sj1side20; %Length of error bars

SlC=DOFs/chi2l*Sj_zcor; %Needs to be a function of Sj magnitude if not on log plot
ShC=DOFs/chi2h*Sj_zcor; %Ends of error bars
CerrL=Sj_zcor-SlC; %Length of error bars
CerrH=ShC-Sj_zcor; %Length of error bars

figure,
%for n = 1:length(Sj1side20(:,1,3))   
semilogy(bandfreqs_IMU(2:end),Sj1side20(23,2:end,3)),hold on
errorbar(bandfreqs_IMU(2:end),Sj1side20(23,2:end,3),errL(23,2:end,3),errH(23,2:end,3),'.');
semilogy(bandfreqs_IMU(2:end),Sj_zcor(23,2:end))
errorbar(bandfreqs_IMU(2:end),Sj_zcor(23,2:end),CerrL(23,2:end),CerrH(23,2:end),'.');
title('Checking Error Bars on Sj from IMU')
hold off
%pause
%end

% save(savename,'datenums_IMU','dropouts_IMU','Disp_std','Ang_std',...
%     'a_means','amp_heave','mean_pitch','mean_roll','bandfreqs_IMU','Sj1side20','Sj_zcor','Sj1side30', 'Sj1side_nf',...
%     'bandwidth','Sjarea20','tsvar20','Sjarea30','tsvar30','Sjarea_nf','tsvar_nf',...
%     'ensembles','bands_IMU', 'padyes','hannme','stds','passes','mean_bands');
% %toc

%Finding noise floor
figure;loglog(bandfreqs(2:end),squeeze(Sj1side20(47,2:end,3)))
hold on
%legend('X','Y','Z')
xlabel('Frequency, Hz')
ylabel('Spectral Energy Density, m^2/Hz')
title(['Spectrum from stationary IMU, ',fileList(ids(23))])
hold off
