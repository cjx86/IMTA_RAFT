%Read, IMU, wave data
%Calculate RAOs, etc.

% Toby Dewhurst, 2015

% Updated by Corey Sullivan 2016, for his deployment

clear;clc;close all;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))
addpath(genpath('../'));
addpath(genpath('../functions'));

%% Load metadata files

% %Load filenames
showmeans=1; %1 for showing mean RAOs. 0 for showing specific instance. %Obsolete...?

%LOADING WAVES DATA AND IMU DATA
load(['..\Observational Experiment\Waves_meta.mat']);
load(['..\Observational Experiment\IMU_meta.mat']);

%Time to show: IMU is the limiting sensor converted to EST
start_datestr=datestr(datenums_IMU(1)); %All
end_datestr=datestr(datenums_IMU(end)); %All

start_datestr_wv=datestr(datenums_IMU(1)-datenum(0,0,0,5,0,0)); %All
end_datestr_wv=datestr(datenums_IMU(end)-datenum(0,0,0,5,0,0)); %All

%LOADING LOADCELL DATA
LC_meta=load('..\Observational Experiment\LC2_meta.mat');

    Sj1side_lc=LC_meta.Sj1_lc;
    bandfreqs_lc=LC_meta.bandf_lc;

    
%Setting the dates to cut the data to based on IMU battery life    
start_datenum=datenums_IMU(1);
end_datenum=datenums_IMU(end);

ttol=datenum(0,0,0,0,1,0); %Within a minute

startkk=find(abs(datenums_IMU-start_datenum)<ttol,1); 
endkk=find(abs(datenums_IMU-end_datenum)<ttol,1);


%% Constants
 
g=9.81;
decl=-15.31667; %Magnetic declination in degrees. (Apply to all measurements)

%% Align and compute/plot RAOs

date_chk_wv = datestr(datenums_Waves);
date_chk_IMU = datestr(datenums_IMU);

%Correct heading/tailing
veldir_f=veldir_f+decl;
Dp=Dp+decl;

%Line up frequencies
    %freq step from the Waves spectra
f_step = .0078125; %Spacing of freqencies from Teledynes Spectra
    %Create an array of these frequencies
freqs_Waves = f_step:f_step:f_step*64;
    %Convert from mm to m. 
    %Square to get from m/sqrt(Hz) to energy spectrum (m^2/Hz).
    %See Hs_area_comp.
Waves_Spectra=1/10^6*Waves_Spectra(:,:).^2; 

Waves_Spectra_all=Waves_Spectra;

%MAYBE USE INTERP HERE 
    %tide1 = interp1(data.tide(:,1),data.tide(:,2),data.t);
    
clear i, clear Waves_Spectra
for i = [1:28 30:ceil(length(Waves_Spectra_all(:,1))/3)]
       
    if i ~= 29 && i<=28
        m = i;
    Waves_Spectra(m,:) = Waves_Spectra_all(3*i-2,:); %This is 0 min output
      
    date_chk_wav(m,:) = date_chk_wv(3*i-2,:);
    
    datenum_Waves(m,:) = datenums_Waves(3*i-2,:);
    
    Hs_t(m,:) = Hs(3*i-2,:);
    
    H13_t(m,:) = H13(3*i-2,:);
    
    Dp_t(m,:) = Dp(3*i-2,:);  
    
    Tp_t(m,:) = Tp(3*i-2,:); 
    
    T13_t(m,:) = T13(3*i-2,:); 
    
    elseif i ~= 29 && i>=30
    m = i-1;
    Waves_Spectra(m,:) = Waves_Spectra_all(3*i-2,:); %This is 0 min output
      
    date_chk_wav(m,:) = date_chk_wv(3*i-2,:);
    
    datenum_Waves(m,:) = datenums_Waves(3*i-2,:);
    
    Hs_t(m,:) = Hs(3*i-2,:);
    
    H13_t(m,:) = H13(3*i-2,:);
    
    Dp_t(m,:) = Dp(3*i-2,:);  
    
    Tp_t(m,:) = Tp(3*i-2,:); 
    
    T13_t(m,:) = T13(3*i-2,:); 
    
    else
    
    end
    
end

clear date_chk_wv, clear datenums_Waves
clear Hs, clear H13, clear Dp, clear Tp,clear T13
datenums_Waves = datenum_Waves;
Hs = Hs_t; clear Hs_t
H13 = H13_t; clear H13_t
Dp = Dp_t; clear Dp_t
Tp = Tp_t; clear Tp_t
T13 = T13_t; clear T13_t
clear datenum_Waves


wvs_spec_check = Waves_Spectra;

    %Interpolate values of IMU frequencies from Actual Frequencies
Waves_Spectra=interp1(freqs_Waves,Waves_Spectra',bandfreqs_IMU);

Waves_Spectra=Waves_Spectra';

wav_spec_interp_chk = Waves_Spectra;

freqs=bandfreqs_IMU; %Generalizing name

%Compute derived wave spectra dimensions
T=1./freqs; %s, Periods relating to banded frequencies
deep=0; %Not in deep water yet (if staff_T decreases as kk increases)
L=zeros(length(T),1);
L2=zeros(length(T),1);
k=zeros(length(T),1);
k2=zeros(length(T),1);
h=8.9579; %m. Depth. Averaged from ADCP Depth Output 

for kk=1:length(T)
    
    if deep==0; %If not in deep water
       [L(kk,1),~] = find_L_dispersion(T(kk),h);
       k(kk,1) = 2*pi./L(kk,1);
        
        if L(kk)<2*h %Check for deep water
            deep=1;
        end
                
    else
        k(kk,1)=(2*pi./T(kk)).^2/g; %Deep water dispersion relation
    end
    
end

%Compute pitch spectrum from heave spectrum
%1--Surge. 2--Sway. 3--Heave. 4--Pitch. 5--Roll. 6--Yaw.

Waves_Spectra=inpaint_nans(Waves_Spectra);%But WATCH OUT for any weird behavior

Waves_Spectra=repmat(Waves_Spectra,1,1,6);  

Wv_Spec=Waves_Spectra;  
%clear Waves_Spectra
%Waves_Spectra = Wv_Spec([1:28 30:end],:,:);


    %Using LWT to convert vertical displ spectra, to surge spectra
Waves_Spectra(:,:,1)=bsxfun(@times,1./tanh(k*h)'.^2,Waves_Spectra(:,:,3));
    %Using LWT to convert vertical displ spectra, to sway spectra
Waves_Spectra(:,:,2)=Waves_Spectra(:,:,1);
    %Using LWT to convert vertical displ spectra, to pitch spectra
Waves_Spectra(:,:,5)=bsxfun(@times,k.^2',Waves_Spectra(:,:,3)); 
    %Using LWT to convert vertical displ spectra, roll spectra
Waves_Spectra(:,:,4)=Waves_Spectra(:,:,5);
    %Using LWT to convert vertical displ spectra, to yaw spectra
Waves_Spectra(:,:,6)=Waves_Spectra(:,:,5);

clrkk=1; %counter to assign colors to plot
leges={'Surge','Sway','Heave','Pitch','Roll','Yaw'};

%% IMU Error bars/Confidence

%IMU confidence
alph=0.05;  % Setting desired confidence
DOFs_IMU=2*bands_IMU*ensembles; % Setting DOF's
[chi2l, chi2h]=chi2_vals(alph,DOFs_IMU);    %Automated way to find CHI^2
    %Error bar creation for entire plot
Sl_IMU=DOFs_IMU*Sj1side20/chi2l; %Needs to be a function of Sj magnitude if not on log plot
Sh_IMU=DOFs_IMU*Sj1side20/chi2h; %Ends of error bars
    
    %Actual sizing of error bars to plot
errL=Sj1side20-Sl_IMU; %Length of error bars
errH=Sh_IMU-Sj1side20; %Length of error bars

%% Waves Confidence: THIS IS EFFECTIVELY GUESSING AT WHAT TELEDYNE DID

ensembles_wv=1;

N_frr = 1200; %Number of samples, I know that I have 1200 pings per sample 
del_ADCP = 1; %seconds, time gap between samples, i.e. ping rate
f_frr=1/(2*N_frr*del_ADCP); %Spacing of the fourier freq produced by F.T.

bands1_wv=floor(0.5/64/f_frr); %Number of bands already averaged to get 64 freqs from 0 to 0.5 Hz.

DOFs_wv=2*bands1_wv*ensembles_wv;

DOFs_wv_noise=2*bands1_wv*ensembles_wv*(length(Waves_Spectra)-1);%*meaning;
[chi2l_wv, chi2h_wv]=chi2_vals(alph,DOFs_wv);

Sl_wv=DOFs_wv*Waves_Spectra(:,2:end-1,3)/chi2l_wv; %Needs to be a function of Sj magnitude if not on log plot
Sh_wv=DOFs_wv*Waves_Spectra(:,2:end-1,3)/chi2h_wv; %Ends of error bars
errsL_wv=Waves_Spectra(:,2:end-1,3)-Sl_wv; %Length of error bars (matrix of all )
errsH_wv=Sh_wv-Waves_Spectra(:,2:end-1,3); %Length of error bars (matrix)

%Plot Pitch spectra of IMU 
figure, 
plot(freqs(2:end-1),mean(Sj1side20(:,2:end-1,4),1),'--','Color',[0 .45 .74],'LineWidth',1.5),hold on
plot(freqs(2:end-1),mean(Sj1side20(:,2:end-1,5),1),'Color',[.85 .33 .1],'LineWidth',1.5);
xlabel('Frequency, Hz') 
ylabel('Spectral Energy Density, [rad^2/Hz]')
set(gca,'FontSize',22)
xlim([0.05 0.5]) 
legend(' Pitch',' Roll')
title('Average of all IMU Pitch and Roll Spectra and Wave Pitch Spectra','FontSize',28)
hold off

% Plot spectra of IMU and Wave to Find noise floors
figure; 
h_1 = semilogy(freqs(2:end-1),mean(Waves_Spectra(:,2:end-1,3),1),'LineWidth',2);
hold on
h2 = semilogy(freqs(2:end-1),mean(Sj1side20(:,2:end-1,3),1));
hh = semilogy(freqs(2:end-1),nanmean(Sj_zcor(:,2:end-1),1),'LineWidth',2);
h3 = errorbar(freqs(2:end-1),mean(Waves_Spectra(:,2:end-1,3),1),mean(errsL_wv,1),mean(errsH_wv,1));
xlabel('Frequency, [Hz]')
ylabel('Spectral Energy Density, [m^2/Hz]')
xlim([0.05 0.5]) 
legend('Wave Spectra','Z Displacement Spectra','Corrected Z Displ Spectra','Wave Error Bars')
uistack(h_1,'top')
set(gca,'FontSize',22)
title('Average of Heave Spectra and Wave Vertical Spectra','FontSize',28)
hold off

figure,
semilogy(freqs(2:end-1),mean(Waves_Spectra(:,2:end-1,3)),'o','LineWidth',1.5);
hold on
semilogy(freqs(2:end-1),mean(Sj1side20(:,2:end-1,3),1),'--','LineWidth',1.5);
semilogy(freqs(2:end-1),nanmean(Sj_zcor(:,2:end-1),1),'LineWidth',1.5);
title('Average of all Z-Displ Spectra and Wave Vertical Spectra','FontSize',22)
xlabel('Frequency, [Hz]')
ylabel('Spectral Energy Density, [m^2/Hz]')
set(gca,'FontSize',18)
legend('Wave Spectra','IMU Vertical Motion Spectra','Heave Motion Spectra (corrected IMU)')
hold off

%% Waves Noise Floor: Visually find and plot the apparent ADCP noise floor 

ns_flr_wv=2*10^-3; 

ns_flrkk=find(freqs>0.43,1);%; index at which noise floor begins, roughly. (Based on visual observations).
    %Find mean value of spectra 
ns_flr_smp=mean(Waves_Spectra(:,ns_flrkk:end,3),2);
ns_flr_smp_max=max(Sh_wv(:,ns_flrkk:end),[],2); %based on peak of noise+error bar

delta = 1; %seconds between samples
nyq=1/(2*delta); %Hz, Nyquist frequency=1/(2*delta)
ns_var=ns_flr_wv*nyq; %m^2. Variance of noise signal (integrating everything below noise floor up to nyquist frequency)
ns_std=sqrt(ns_var); %standard dev of noise signal
ns_pk2pk=2*sqrt(2)*ns_std/2; %Peak to Peak ?

%OR
% ns_pk2pk=8*ns_std; %http://www.dspguide.com/ch2/2.htm
ns_pk2pk=6.6*ns_std; %https://www.youtube.com/watch?v=-KcODSYXiZA. "We want to convert from RMS noise into peak to peak noise. Because noise has a Gaussian distribution, to span 99.9% of the waveform (approximately peak to peak), we need 6.6 standard deviations. Because the average value for the noise is 0, RMS equals the standard deviation, so we can just multiply the RMS noise by 6.6."

DOFs_wvns=2*bands1_wv*ensembles_wv;
[chi2l_wvns, chi2h_wvns]=chi2_vals(alph,DOFs_wvns);

Sl_wvns=DOFs_wv*mean(Waves_Spectra(:,:,3),1)/chi2l_wvns; %Needs to be a function of Sj magnitude if not on log plot
Sh_wvns=DOFs_wvns*mean(Waves_Spectra(:,:,3),1)/chi2h_wvns; %Ends of error bars

ns_flr_errh=DOFs_wvns*ns_flr_wv/chi2h_wvns;
line([0 nyq],[ns_flr_errh ns_flr_errh],'LineStyle','--','Color','g')

%Ignore values whose error bars dip below noise floor
Waves_Spectra_tmp=squeeze(Waves_Spectra(:,:,3));
% thresh=ns_flr_smp_max
thresh=DOFs_wv*ones(size(ns_flr_smp))*ns_flr_wv./chi2h_wv; %Ends of error bars
% Waves_Spectra_tmp(bsxfun(@le,Sl_wv,thresh))=NaN;
thresh1 = .002.*ones(length(thresh),1);
Waves_mask=bsxfun(@le,Sl_wv,thresh); %Order is the 1st array, fcn then the 2nd
Waves_mask1=bsxfun(@le,Sl_wv,thresh1); %Returns 0 if Sl_wv is greater than 0.01
% Waves_Spectra(:,:,3)=Waves_Spectra_tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RAOs: Calculating the RAOs from the IMU and the Waves

RAOs_uncorr=sqrt(abs(Sj1side20(:,2:end-1,:))./abs(Waves_Spectra(:,2:end-1,:)));

COM_corrected = 1;

if COM_corrected
   
    Sj_z_uncorr = Sj1side20(:,:,3);
    Sj1side20(:,:,3) = Sj_zcor;
    disp('The Z disp is corrected to the Center of Mass')
    pause(0.5)
    
else
    
    disp('The Z disp IS NOT corrected to the Center of Mass')
    pause(1)
    
end



    %General RAOs from IMU data 
RAOs_raw=sqrt(abs(Sj1side20(:,2:end-1,:))./abs(Waves_Spectra(:,2:end-1,:)));
    %Checking for non-real data
RAOs_raw(real(RAOs_raw)==Inf)=NaN;
RAOs_uncorr(real(RAOs_uncorr)==Inf)=NaN;

    %Simplifying notation

RAOs=RAOs_raw;
RAOs1 = RAOs_raw;
%look1 = RAOs(:,:,3); Used to check what the Waves_Mask does to the RAOs

    %Ignoring RAO if waves data drops below noise floor
RAOs(repmat(Waves_mask,1,1,6))= NaN; %
RAOs_uncorr(repmat(Waves_mask,1,1,6))= NaN; %
    %Playing with the noise floor manually.
RAOs1(repmat(Waves_mask1,1,1,6))= NaN;

chk1 = RAOs(:,:,3);
chk2 = RAOs1(:,:,3);

%Compute mean RAOs for certain conditions
cdn_incld=ones(length(Dp),1);
%cdn_incld(excldkk)=0;

degwin=10; %Allowable margin for direction. degwin/2 on each side of axis

%pause
disp('YOU NEED TO FIND ACTUAL MOORING HEADING THIS IS A GUESS')
moorax=-30+decl; %deg. Mooring axis (from compass readings)

headonly=0; %0 includes following seas in mean 

if headonly %This is a logical, if headonly = 1 the if staement executes
    cdn_head=(Dp-veldir_f')<degwin %Condition for alignment of waves and tides
else
    cdn_head=ones(length(Dp),1);
end

    %Max and min possible direction?
Dp_min=0;%moorax-degwin/2;
Dp_max=360;%moorax+degwin/2;

    
    %Max and min H13
H13_min=0.0;
H13_max=0.6; %Inf

% Addressing any negative angles in Dp 
for ii = 1:length(Dp)
    if Dp(ii) < 0
        Dp(ii) = Dp(ii)+360;
    end
end
    
    %Allowing waves from any direction
cdn_Dp=(Dp>=Dp_min&Dp<=Dp_max)|(Dp>=Dp_min+180&Dp<=Dp_max+180); %Allowing either direction (but heading condition)
% cdn_CD=(veldir_f>=CD_min&veldir_f<=CD_max)|(veldir_f>=CD_min+180&veldir_f<=CD_max+180); %Allowing either direction (but heading condition)
% cdn_CD=(veldir_f>=CD_min&veldir_f<=CD_max); %Allowing only flooding tides


velmag_chk_sum(1) = sum(velmag_f(1:19));
veldir_chk_sum(1) = sum(veldir_f(1:19));

for jj = 1:length(datenums_Waves)
    
    if jj ~= 29 && jj <= 28
        n = jj;
    mag_ind(n) = 60*jj+1;
    mag_ind2(n) = 60*jj+19;
    velmag_chk_sum(n+1) = sum(velmag_f(mag_ind(n):mag_ind2(n)))';
    veldir_chk_sum(n+1) = sum(veldir_f(mag_ind(n):mag_ind2(n)))';
    
    elseif jj~= 29 && jj >= 30
         n = jj-1;
    mag_ind(n) = 60*jj+1;
    mag_ind2(n) = 60*jj+19;
    velmag_chk_sum(n+1) = sum(velmag_f(mag_ind(n):mag_ind2(n)))';
    veldir_chk_sum(n+1) = sum(veldir_f(mag_ind(n):mag_ind2(n)))';
    
    else
        
    end
    
end

velmag_chk = velmag_chk_sum'./length(velmag_f(mag_ind(n):mag_ind2(n)));
veldir_chk = veldir_chk_sum'./length(veldir_f(mag_ind(n):mag_ind2(n)));

for ii = 1:length(veldir_chk)
    if veldir_chk(ii) < 0
        veldir_chk(ii) = veldir_chk(ii)+360;
    end
end

%Setting max and min current directions relative to raft that sets the
%condition that the current and waves are aligned
    
CD_min=moorax-degwin/2; %CD is direction current is headed(as per load cells)
CD_max=moorax+degwin/2;


    %Allowing for either flow direction
%cdn_CD1=(veldir_f>=CD_min & veldir_f<=CD_max); 
%cdn_CD_eb1=(veldir_f>=(CD_min-180) & veldir_f<=(CD_max-180)); %Allowing ebb
cdn_CD1=(veldir_chk>=CD_min & veldir_chk<=CD_max); 
%360 added to keep the ebb flow direction a positive degree
cdn_CD_eb1=(veldir_chk>=(360+CD_min-180) & veldir_chk<=(360+CD_max-180));

%Done to make it run until I figure out wtf is going on
cdn_CD = ones(length(cdn_Dp),1);
cdn_CD_eb = ones(length(cdn_Dp),1);

    %Ensuring H13 is between min and max
cdn_H13=(H13>=H13_min&H13<=H13_max);%H13>=Hs_min;
    
    %Setting a counter of all the data from the IMU
kks=(1:length(datenums_IMU))';
    % Making sure that IMU data is between known start and end dates
cdn_date_RAO=and(kks>=startkk,kks<=endkk);
    %Setting a counter of all the data from the waves from ADCP
kks_wv=(1:length(datenums_Waves))';
    % Making sure that IMU data is between known start and end dates
cdn_date_wv=and(kks_wv>=startkk,kks_wv<=endkk);

    %IF the starting date is the ending date set something weird in the
    %logicals
if startkk==endkk
    
    cdn_all_RAO=zeros(length(Disp_std),1);
    cdn_all_RAO(startkk)=1;
    
    cdn_all_wv=zeros(length(H13),1);
    cdn_all_wv(startkk)=1;
    
else 
        %Basically saying all of the condidions need to be met for
        %cdn_all_RAO to be relevent
    %removed so it would run until I figure it out ".*cdn_CD(1:length(cdn_date_RAO))"
    cdn_all_RAO=cdn_head(1:length(cdn_date_RAO)).*cdn_Dp(1:length(cdn_date_RAO)).*cdn_date_RAO.*cdn_H13(1:length(cdn_date_RAO)).*cdn_incld(1:length(cdn_date_RAO)); %Triple AND function
    %".*(cdn_CD+cdn_CD_eb)"
    cdn_all_wv=cdn_Dp.*cdn_date_wv.*cdn_H13.*cdn_incld; %Quadruple AND function   
end

RAOs_filt = tsd_filterA(RAOs,5,1,1,NaN);

    %applying the mean using the mask to each spectra
Waves_Spectra_mean=nanmean(Waves_Spectra(logical(cdn_all_RAO),:,:),1);
Sj1side_mean=nanmean(Sj1side20(logical(cdn_all_RAO),:,:),1);
Sj1side_mean_uncor=nanmean(Sj_z_uncorr(logical(cdn_all_RAO),:),1);
Sj1side_mean_nf=nanmean(Sj1side_nf(logical(cdn_all_RAO),:,:),1);

RAOs_mean=nanmean(RAOs(logical(cdn_all_RAO),:,:),1);
RAOs_mean2 = nanmean(RAOs_filt(logical(cdn_all_RAO),:,:),1);
meaning=sum(cdn_all_RAO);

    %Calculate RAO confidence
[Flow,Fhigh]=Ferror_vals(alph,DOFs_IMU,DOFs_wv);
Lerr=(RAOs_filt/sqrt(Fhigh));
Uerr=(RAOs_filt/sqrt(Flow));

[Flow_mn,Fhigh_mn]=Ferror_vals(alph,DOFs_IMU*meaning,DOFs_wv*meaning);
Lerr_mn=(RAOs_mean./sqrt(Fhigh_mn));
Uerr_mn=(RAOs_mean./sqrt(Flow_mn));
Lerr_mn=Lerr_mn.*RAOs_mean./RAOs_mean;
Uerr_mn=Uerr_mn.*RAOs_mean./RAOs_mean;

%% Plotting HEAVE Spectra and RAOs

%choose a day
des_day = 17;
des_Hr = 23; %Hour 3 in EST
IMUdate_desir = datenum(2016,01,des_day,des_Hr,00,00);

    %Find closest this date in the datenums array (min,sec could be off)
[~, ind] = min(abs(datenums_IMU-IMUdate_desir));

figure,
subplot(2,1,1); 
semilogy(freqs(2:end-1),Waves_Spectra(ind,2:end-1,3));
hold on
semilogy(freqs(2:end-1),Sj_z_uncorr(ind,2:end-1));
semilogy(freqs(2:end-1),Sj1side20(ind,2:end-1,3));
%errorbar(freqs(2:end-1),Waves_Spectra(:,2:end-1,3))
title(['Z-Displ Spectra and Wave Vertical Spectra ', datestr(IMUdate_desir,6),' ',num2str(des_Hr)])
xlabel('Frequency, Hz')
ylabel('Spectral Energy Density, m^2/Hz')
legend('Wave Spectra','IMU Spectra','Corrected IMU')
hold off

subplot(2,1,2)
plot(bandfreqs_IMU(2:end-1),RAOs(ind,:,3),bandfreqs_IMU(2:end-1),RAOs1(ind,:,3),'--')
title(['Heave RAO for ', datestr(IMUdate_desir,6),' ',num2str(des_Hr),' (not corrected to C.O.M. yet)'])
xlim([0.05 0.5])
xlabel('Frequency, Hz')
linkaxes(findall(gcf,'type','axes'), 'x');
ylabel('RAO')
legend('Calculated threshold','Visually Applied Threshold')

%%Plot the average RAO over the course of specific dates
    %These ones are during the storm event
avg_RAO_st = 15;
RAO_hr_st = 21;
avg_RAO_end = 19;
RAO_hr_end = 14;
date2avg1 = datenum(2016,01,avg_RAO_st,RAO_hr_st,00,00);
date2avg2 = datenum(2016,01,avg_RAO_end,RAO_hr_end,00,00);


    %Find closest this date in the datenums array (min,sec could be off)
[~, inde] = min(abs(datenums_IMU-date2avg1));
[~, inde2] = min(abs(datenums_IMU-date2avg2));

figure,plot(bandfreqs_IMU(2:end-1),RAOs_mean(:,:,3),'LineWidth',1.5),hold on

xlabel('Frequency, [Hz]')
ylabel('Response Amplitude')
xlim([0.05 0.5])
ylim([0.3 1.3])
set(gca,'YTick',[0:0.1:1.3])
set(gca,'FontSize',22)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca, 'FontSize', 22)
title('Date Averaged, Deployment Heave RAOs','FontSize',28)

figure,plot(bandfreqs_IMU(2:end-1),nanmean(RAOs(inde:inde2,:,4),1),'-','LineWidth',1.5)
hold on,plot(bandfreqs_IMU(2:end-1),nanmean(RAOs(inde:inde2,:,5),1),'--','LineWidth',1.5)
xlabel('Frequency, [Hz]')
ylabel('Response Amplitude')
xlim([0.05 0.5])
set(gca,'YTick',[0:0.5:3])
legend(' Pitch',' Roll')
set(gca,'FontSize',22)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca, 'FontSize', 22)
title('Date Averaged, Deployment Pitch and Roll RAOs','FontSize',28)

dim_pit_RAO = sqrt(abs(Sj1side20(inde:inde2,2:end-1,4))./abs(Waves_Spectra(inde:inde2,2:end-1,3)));
dim_pit_RAO(real(dim_pit_RAO)==Inf)=NaN;
dim_roll_RAO = sqrt(abs(Sj1side20(inde:inde2,2:end-1,5))./abs(Waves_Spectra(inde:inde2,2:end-1,3)));
dim_roll_RAO(real(dim_roll_RAO)==Inf)=NaN;

figure,plot(bandfreqs_IMU(2:end-1),nanmean(dim_pit_RAO,1),'-','LineWidth',1.5)
hold on,plot(bandfreqs_IMU(2:end-1),nanmean(dim_roll_RAO,1),'--','LineWidth',1.5)
xlabel('Frequency, [Hz]')
ylabel('Response Amplitude [rad/m]')
xlim([0.05 0.5])
legend(' Pitch',' Roll')
set(gca,'YTick',[0:0.1:0.6])
set(gca,'FontSize',22)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca, 'FontSize', 22)
title('Date Averaged, Deployment Dimensional Pitch and Roll RAOs','FontSize',28)

%PEAKS ARE SEEN @ 0.16 AND 0.224 Hz in single and averages
    %@ 3.8 and 4.5 sec
    %Printing out the wavelengths seen at mean water @ the site
L_peak1 = find_L_dispersion(1/.224,h);
Length1 = ([num2str(L_peak1./.3048),' Feet'])
L_peak2 = find_L_dispersion(1/.2632,h);
Length2 = ([num2str(L_peak2./.3048),' Feet'])

%% Significant amplitudes of low-frequency motions

As_lf=sqrt(2*trapz(freqs(4:6),nanmean(Sj1side20(logical(cdn_all_RAO),4:6,:),1),2));

%% Current Rose: Showing the magnitude weighted direction of the current

lng(kk)=length(H13);
veldir_f(veldir_f>360)=NaN;

[N_CD, edges_CD, bin_CD]=histcounts(veldir_f,36);

for CDkk=1:length(N_CD)
    
    CM_bin(CDkk)=nansum(velmag_f(bin_CD==CDkk));
    
end

figure;title('Current Direction Probability weighted by magnitude.')
[Xpol,Ypol]=pol2cart(((edges_CD(2:end))*pi./180),CM_bin/sum(CM_bin));
compass(Xpol,Ypol)

%Show raft orientation: 
graft = hgtransform;
raft = rectangle('position',0.05*[-1,-1,2,2],'Parent',graft);
graft.Matrix = makehgtform('zrotate',moorax*pi/180);
view([90 -90])

dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@polarangs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot RAOs (all together)

mkks = 1:1:6;

%Figures;
%Set figure preferences
mpos=get(0,'MonitorPositions');
if sum(size(mpos))>4 %More than one monitor
    figure('units','normalized','Position',[1.2 0.2 0.6 0.4]);
else
    figure('units','normalized','Position',[0.2 0.2 0.6 0.4])
end

    %Finding which date to plot the RAOs of
viewdatestr = datenum(IMUdate_desir);
viewdatenum=datenum(viewdatestr);
viewkk=find(abs(datenums_IMU-viewdatenum)<ttol,1); 

    %Plotting the RAOs for each component of the motion
for mkk=mkks
    
    hold on
    if showmeans
        plot(freqs(2:end-1),RAOs_mean(:,:,mkk),'DisplayName',leges{mkk})
    else
        plot(freqs(2:end-1),RAOs(viewkk,:,mkk),'DisplayName',leges{mkk})
    end
    % clrkk=clrkk+1; %counter to assign colors to plot
    
end

set(gca, 'ColorOrderIndex', 1)
    
startdatestr_p = datestr(datenums_IMU(1)-datenum(0,0,0,5,0,0));
enddatestr_p = datestr(datenums_IMU(end)-datenum(0,0,0,5,0,0));


    %Setting the title of the RAO plot and checking whether its means or
    %one instance
if showmeans
%     title(['Mean RAOs'])% for ' num2str(Dp_min) '<Dp<' num2str(Dp_max) ', ' num2str(CD_min) '<veldir_f<' num2str(CD_max)])
    if headonly
        titlestr=strvcat(['Mean RAOs for Wave Direction and Current Dir. aligned and within ' num2str(degwin) ' deg. of mooring axis. ' dep],[startdatestr_p ' to ' enddatestr_p '. ' num2str(meaning) ' samples'])
    else
        titlestr=strvcat(['Mean RAOs for Wave Direction and Current Mag.'],[startdatestr_p ' to ' enddatestr_p '. ' num2str(meaning) ' samples'])
    end
    title(titlestr)

else %If the means of the RAOs are not what is desired this plots one instance 
    title(viewdatestr)
    plot(freqs(2:end-1),Waves_Spectra(startkk:endkk,2:end-1,1),'--k','DisplayName','Wave (Heave)')%,'Color',pltcs{clrkk}) %Multiply freqs by 2pi to plot over rad/s
   
    for mkk=mkks
        
        if mkk==3
            plot(freqs(2:end),Sj_zcor(startkk:endkk,2:end-1),':','DisplayName',['IMU.' leges{mkk}])%,'Color',pltcs{clrkk}) %Multiply freqs by 2pi to plot over rad/s
        else
            plot(freqs(2:end),Sj1side20(startkk:endkk,2:end-1,mkk),':','DisplayName',['IMU.' leges{mkk}])%,'Color',pltcs{clrkk}) %Multiply freqs by 2pi to plot over rad/s
            % clrkk=clrkk+1; %counter to assign colors to plot
        end
        
    end
    
    %Or, if you want to see what you're averaging:
    figure;plot(freqs(2:end-1),RAOs(logical(cdn_all_RAO),:,3),':','DisplayName',['IMU.' leges{mkk}])%,'Color',pltcs{clrkk}) %Multiply freqs by 2pi to plot over rad/s
    
end

legend(gca,'show')
xlabel('Frequency, Hz');
ylabel('Response Amplitude Operator')
xlim([0 0.5])
ylim([0 1.2])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Show Wave Spectrum and Heave Response

figure;
loglog(freqs(2:end-1),Waves_Spectra_mean(1,2:end-1,3),'--','DisplayName','Wave')
    hold on
    loglog(freqs(2:end-1),Sj1side_mean(1,2:end-1,3),'DisplayName','Raft, IMU Z-Displ')
    loglog(freqs(2:end-1),Sj1side_mean_uncor(1,2:end-1),'DisplayName','Raft, Heave(Corrected IMU)')
    %loglog(bandfreqs_lc,mean(Sj1side_lc(logical(cdn_all_RAO),:),1),'DisplayName','Load Cell 2')%,'Color',pltcs{clrkk}) %Multiply freqs by 2pi to plot over rad/s
%     errorbar(freqs(2:end),Sj1side20(startkk,2:end,3),errL_IMU(startkk,2:end,3),errH_IMU(startkk,2:end,3))
legend(gca,'show')
title('Mean Spectral Energy')
xlabel('Frequency, Hz');
ylabel('Spectral Energy Density, m^2/Hz')
xlim([0 1])
hold off

%% Show effect of current magnitude: Shows how RAOs change with current mag

CM_edges=0:0.1010:max(velmag_f);

figure('units','normalized','Position',[0.2 0.2 0.6 0.4])

velmag_chk_sum(1) = sum(velmag_f(1:19));

for jj = 1:length(datenums_Waves)-1
    
    mag_ind(jj) = 60*jj+1;
    mag_ind2(jj) = 60*jj+19;
    velmag_chk_sum(jj+1) = sum(velmag_f(mag_ind(jj):mag_ind2(jj)))';
    
end

velmag_chk = velmag_chk_sum'./length(velmag_f(mag_ind(jj):mag_ind2(jj)));

cmm = 1;
for CMkk=2:length(CM_edges)
    
    cdn_CM=and(velmag_chk >= CM_edges(CMkk-1),velmag_chk <= CM_edges(CMkk));
    sum(cdn_CM);
    RAOs_CMkk=nanmean(RAOs(logical(cdn_all_RAO.*cdn_CM),:,:),1);
    
    hey1 = subplot(2,1,1);
    for mkk = 1
        hold on
        plot(freqs(2:end-1),RAOs_CMkk(:,:,mkk));%,'Color',pltcs{clrkk}) %Multiply freqs by 2pi to plot over rad/s
        clrkk=clrkk+1; %counter to assign colors to plot
    end
    
    title('Mean Surge (X-Axis Motion) RAOs based on different Velocities')
    xlabel('Frequency (Hz)')
    ylabel('RAO')
    
    subplot(2,1,2)
    for mkk = 2
        hold on
        plot(freqs(2:end-1),RAOs_CMkk(:,:,mkk),'DisplayName',[num2str(CM_edges(CMkk-1)), ' < Cur_{mag} < ',num2str(CM_edges(CMkk))])%,'Color',pltcs{clrkk}) %Multiply freqs by 2pi to plot over rad/s
        clrkk=clrkk+1; %counter to assign colors to plot
    end
cmm = cmm+1;
end

legend(hey1,'    0 < Cur_{mag} < 0.101','0.101 < Cur_{mag}  <0.202','0.202 < Cur_{mag} < 0.303','0303 < Cur_{mag} < 0.404','0.404 < Cur_{mag} < 0.505')
legend(gca,'show')
title('Mean Sway (Y-Axis Motion) RAOs based on different Velocities')
xlabel('Frequency (Hz)')
ylabel('RAO')
set(gca, 'ColorOrderIndex', 1)
hold off

%% Show effect of H13: Shows how RAOs change with H_sig seen during the deployment


Hs_edges=0:max(H13)/6:max(H13);

figure('units','normalized','Position',[0.2 0.2 0.6 0.4])
m = 1;

for Hskk=2:length(Hs_edges);
    
    cdn_Hs=and(Hs>=Hs_edges(Hskk-1),Hs<=Hs_edges(Hskk));
    
    cdn_H13=and(H13>=Hs_edges(Hskk-1),H13<=Hs_edges(Hskk));
    sum(cdn_H13)
    
    RAOs_H13kk=nanmean(RAOs(logical(cdn_all_RAO.*cdn_H13),:,:),1);
    
    RAOsUC_H13kk=nanmean(RAOs_uncorr(logical(cdn_all_RAO.*cdn_H13),:,:),1);
    
    out_heave(Hskk,:,:) = RAOs_H13kk(:,:,3);
    
    out_heaveUC(Hskk,:,:) = RAOsUC_H13kk(:,:,3);
    
    for mkk=3%mkks
        hold on
            plot(freqs(2:end-1),RAOs_H13kk(:,:,mkk),'DisplayName',[sprintf('%0.3f',Hs_edges(Hskk-1)) ' < H_{1/3} < ' sprintf('%0.3f',(Hs_edges(Hskk))) ' m'])%,'Color',pltcs{clrkk}) %Multiply freqs by 2pi to plot over rad/s
        % clrkk=clrkk+1; %counter to assign colors to plot
    end
    
    
%     pause
m = m+1;
end


save('RAO_dep','freqs','out_heave','out_heaveUC','H13','RAOs','RAOs_uncorr','Sj1side20','Sj_z_uncorr','Waves_Spectra')

%ylim([0 1.1])
title(['Mean ' leges(:,mkk) 'RAOs for different Significant Wave Heights'])
ylabel('RAOs')
legend(gca,'show','Location','NorthWest')
set(gca, 'ColorOrderIndex', 1)
hold off



%% Plot Individual Mean RAOs and Spectra on same plot

% % Checking units of Waves_Spectra
%Hs_area_comp=[4*sqrt(trapz(freqs,Waves_Spectra(viewkk,:,3))) Hs(viewkk)]; 
    % If they are almost equal shows Waves_Spectra is in correct units
 
%RAOs-separate
figure
%Title
if headonly
    titlestr=strvcat(['Mean RAOs for Dp, veldir_f aligned and within ' num2str(degwin) ' deg. of mooring axis. '],[start_datestr ' to ' end_datestr '. ' num2str(meaning) ' samples'])
else
    titlestr=strvcat(['Mean RAOs for Dp, veldir_f within ' num2str(degwin) ' deg. of mooring axis. '],[start_datestr ' to ' end_datestr '. ' num2str(meaning) ' samples'])
end

% mkks=[3 5 1 2 4]; %Sets order of graphs/DOFs
ylims=[.5 .375 1.5 1.25 1.5 .75];%Should be in order of [1-6], i.e. surge-yaw.
y2lims =[0.25 .25 0.05 2E-3 2E-3 2e-3];
Slabels={'S_\zeta, m^2/Hz','S_{eta_3}, m^2/Hz','S_\eta, m^2/Hz','S_{eta_3}, rad^2/Hz','S_\theta, rad^2/Hz','S_{eta_3}, rad^/Hz'};
skip=2; %For error bars

%Waves_Spectra_avg=mean([Waves_Spectra(cdn_date_RAO,:,:);Waves_Spectra(cdn_date_RAO,:,:)],1); %For showing average spectrum

    %Plotting the individual RAOs
for kk=1:length(mkks)
    
    subplot(length(mkks)*100+10+kk)
    
    %h1 is the RAO of each part    %hf is the corresp. wave spectra
    [axs,h1,hf]=plotyy(freqs(2:end-1),RAOs_mean(:,:,mkks(kk)), ... % Data for first plot
    [freqs(2) freqs(2:end-1) freqs(end)],[0 Waves_Spectra_mean(:,2:end-1,mkks(kk)) 0], ... % Data for second plot
    @(X, Y) plot(X, Y), ... % Function to use for first plot
    @(X, Y) fill(X, Y, [0.9 0.9 1],'EdgeColor',[0.7 0.7 0.7])); % Function to use for second plot    
    
    hold on
    %Without Wave Spectrum 
    %     plot(freqs(2:end-1),RAOs_mean(:,2:end-1,mkks(kk)),'DisplayName','Surfaced'),grid
    
    ylim([0 ylims(mkks(kk))])
    eb1=errorbar(freqs(2:skip+1:end-1),RAOs_mean(:,1:skip+1:end,mkks(kk)),RAOs_mean(:,1:skip+1:end,mkks(kk))-Lerr_mn(:,1:skip+1:end,mkks(kk)),Uerr_mn(:,1:skip+1:end,mkks(kk))-RAOs_mean(:,1:skip+1:end,mkks(kk)),'.');
   
    
ax1=axs(1);
ax1.YTick=0:ylims(mkks(kk))/2:ylims(mkks(kk));
    ylabel(leges(mkks(kk)))
    hold all

            if kk==1

        axtop=LinkTopAxisData(freqs(2:3:end),round(1./freqs(2:3:end)*10)/10);
       legend([h1 hf],{'Mean RAO','Waves Spectra'})
            else
        axtop=LinkTopAxisData([0],[0]); %Hack. (To bring axs(1) to top

            end

     ax2=axs(2);
     
     set(axs(1),'YLim',[0 ylims(mkks(kk))])
     set(axs(2),'YLim',[0 y2lims(mkks(kk))])
     
     ax2.YTick=0:y2lims(mkks(kk))/2:y2lims(mkks(kk));
     ylb=ylabel(Slabels{mkks(kk)});
     set(ylb, 'Units', 'Normalized', 'Position', [1.07, 0.5, 0]);
     
end

axes(ax1);
xlabel('Frequency, Hz')
linkaxes(findall(gcf,'type','axes'), 'x');
xlim([1/20 freqs(ns_flrkk)])
set(gcf,'Color','w')




%% Check Data Alignment: Check LC, ADCP and IMU alignment and plot

LCind1 = 28; %Matches 01/15/16 hour: 21 UTC 16 EST
LCind2 = 117; %Matches 01/19/16 hour 14 UTC 09 EST

    %Convert the LC dates to EST
datenums_LC = LC_meta.yd(LCind1:LCind2)-datenum(0,0,0,5,0,0);
%date_chk_LC_UTC = datestr(LC_meta.yd(LCind1:LCind2));
date_chk_LC_EST = datestr(datenums_LC);

     %Filtering LC data
[LC_filt,LC_filt_mean,~,~,numgg,~]=stddev_filter(LC_meta.mn_lb(LCind1:LCind2),6,3);

%Shorten LC std_dev
LC_std = LC_meta.std_lb(LCind1:LCind2);

figure; %Plot current effects  %FIGURE (4)
subplot(311),plot(datenums_cur,velmag_f,'DisplayName','Current Magnitude'),datetickzoom,legend('Show')
ylabel('Velocity Mag (m/s)')
subplot(312),plot(datenums_cur,veldir_f,'r','DisplayName','Current Direction'),datetickzoom,legend('Show')
ylabel('Velocity Dir (degrees)')
subplot(313),plot(datenums_LC,LC_filt,'g','DisplayName','LC2 Mean')
hold on
% plot(Disp_std(:,3)/mean(Disp_std(:,3)),'DisplayName','IMU SD')
datetickzoom
%plot(datenums_LC(logical(cdn_all_wv.*cdn_CD)),LC_filt(logical(cdn_all_wv.*cdn_CD)),'og')
%plot(datenums_LC(logical(cdn_all_wv.*cdn_CD_eb)),LC_filt(logical(cdn_all_wv.*cdn_CD_eb)),'.g')
legend('Show')
xlabel('Date')
ylabel('Mean Tension, lb')
linkaxes(findall(gcf,'type','axes'), 'x');
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
suptitle('Current Effects on Tension')
hold off

%% Check How Different Factors Affect Tensions

%Plot of water level and current direction
figure,subplot(2,1,1),plot(datenums_cur,velmag_f,'LineWidth',1.5),hold on
datetickzoom('x','mm/dd','keepticks')
set(gca,'FontSize',18)
ylim([0 0.7])
ylabel('Velocity [m/s]')
title('Current Magnitude')
plot(datenums_LC(12).*ones(length(0:.1:.8)),[0:.1:.8],'k--','Linewidth',1.5)
plot(datenums_LC(36).*ones(length(0:.1:.8)),[0:.1:.8],'k--','Linewidth',1.5)
plot(datenums_LC(62).*ones(length(0:.1:.8)),[0:.1:.8],'k--','Linewidth',1.5)
plot(datenums_LC(87).*ones(length(0:.1:.8)),[0:.1:.8],'k--','Linewidth',1.5),hold off

subplot(2,1,2),plot(datenums_LC,LC_filt,'Color',[.85 .33 .1],'LineWidth',1.5),hold on
plot(datenums_LC(12).*ones(length(200:10:300)),[200:10:300],'k--','Linewidth',1.5)
plot(datenums_LC(36).*ones(length(200:10:300)),[200:10:300],'k--','Linewidth',1.5)
plot(datenums_LC(62).*ones(length(200:10:300)),[200:10:300],'k--','Linewidth',1.5)
plot(datenums_LC(87).*ones(length(200:10:300)),[200:10:300],'k--','Linewidth',1.5),hold off
datetickzoom('x','mm/dd','keepticks')
ylim([210 280])
ylabel('Tension [lbs]')
xlabel('Date [mm/dd]')
set(gca,'FontSize',18)
linkaxes(findall(gcf,'type','axes'), 'x');
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
title('Mean Load Cell Tension')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LOW WAVE EVENT INFORMATION OVERLAYED

low_wavst = 18;
low_wavst_hr = 03;
hi_wavend = 19;
low_wavend_hr = 09;
dat1 = datenum(2016,01,low_wavst,low_wavst_hr,00,00);
dat2 = datenum(2016,01,hi_wavend,low_wavend_hr,00,00);

    %Find closest this date in the datenums array (min,sec could be off)
[~, in] = min(abs(datenums_cur-dat1));
[~, in2] = min(abs(datenums_cur-dat2));

[~, inwv] = min(abs(datenums_Waves-dat1));
[~, in2wv] = min(abs(datenums_Waves-dat2));

[~, inLC] = min(abs(datenums_LC-dat1));
[~, in2LC] = min(abs(datenums_LC-dat2));

% figure; %Plot current effects during LOW WAVE EVENTS 1/18 09 to 1/19 09 
% subplot(311),[axx,~,~] = plotyy(datenums_Waves(inwv:in2wv),H13(inwv:in2wv),datenums_cur(in:in2),velmag_f(in:in2));
% set(axx(2),'YLim',[0 0.6])
% ylabel(axx(1),'Sig Wave Ht (m)')
% ylabel(axx(2),'Current Magnitude (m/s)')
% datetickzoom
% legend('Sig Wave Ht','Current Magnitude','Location','NorthEast')
% subplot(312),[a_x,~,~] = plotyy(datenums_Waves(inwv:in2wv),Dp(inwv:in2wv),datenums_cur(in:in2),veldir_f(in:in2));
% datetickzoom
% ylabel(a_x(1),'Wave Dir (deg)')
% ylabel(a_x(2),'Current Dir (deg)')
% legend('Wave Direction','Current Magnitude','Location','NorthWest')
% subplot(313),[a_x2,~,~] = plotyy(datenums_LC(inLC:in2LC),LC_filt(inLC:in2LC),datenums_cur(in:in2),w_lev(in:in2));
% set(a_x2(1),'YLim',[200 250])
% set(a_x2(2),'YLim',[7 11])
% legend('LC2 Mean','Water Level','Location','SouthWest')
% hold on
% datetickzoom
% legend('Show')
% xlabel('Date')
% ylabel(a_x2(1),'Mean Tension (lb)')
% ylabel(a_x2(2),'Water Level (m)')
% linkaxes(findall(gcf,'type','axes'), 'x');
% dcm_obj = datacursormode(figure(gcf));
% set(dcm_obj,'UpdateFcn',@datestrcurs)
% suptitle('Overlay Effects from Low Wave Event on Tension 1/18 - 1/19')
% hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %HIGH WAVE EVENT INFORMATION OVERLAYED
hi_wavst = 15;
hi_wavst_hr = 16;
hi_wavend = 17;
hi_wavend_hr = 11;
dat3 = datenum(2016,01,hi_wavst,hi_wavst_hr,00,00);
dat4 = datenum(2016,01,hi_wavend,hi_wavend_hr,00,00);

    %Find closest this date in the datenums array (min,sec could be off)
[~, inh] = min(abs(datenums_cur-dat3));
[~, inh2] = min(abs(datenums_cur-dat4));

[~, inhwv] = min(abs(datenums_Waves-dat3));
[~, inhwv2] = min(abs(datenums_Waves-dat4));

[~, inhLC] = min(abs(datenums_LC-dat3));
[~, inhLC2] = min(abs(datenums_LC-dat4));

%% Look at Tension Changes for different Sig Have Heights
% 
% figure;
% ckk = 1;
% for Hskk=2:length(Hs_edges);
% 
%     
%     LC_ckk = LC_filt(:,logical(cdn_Hs(:,ckk).*cdn_all_RAO));
%     
%     for mkk=3%mkks
%         hold on
%             plot(datenums_LC(logical(cdn_Hs(:,ckk).*cdn_all_RAO)),LC_ckk,'DisplayName',[num2str(Hs_edges(Hskk-1)) ' < H13 < ' num2str(Hs_edges(Hskk)) ' m'])%,'Color',pltcs{clrkk}) %Multiply freqs by 2pi to plot over rad/s
%         % clrkk=clrkk+1; %counter to assign colors to plot
%     end
% %     pause
% ckk = ckk+1;
% end
% datetickzoom
% %ylim([0 1.1])
% title([leges(mkk) 'Mean Tensions for different Significant Wave Heights'])
% ylabel('Tension (lbs)')
% legend(gca,'show')
% set(gca, 'ColorOrderIndex', 1)

%% RAO of LC Tension and Current Magnitude

ttol=datenum(0,0,0,0,0,10); %Within 30 sec
m = 1;

ttolW=datenum(0,0,0,0,10,0); %Within 30 sec
n = 1;
for ii = 1:length(datenums_LC)
    
    for jj = 1:length(datenums_cur)
        
        if abs(datenums_cur(jj)-datenums_LC(ii))<ttol        
            current_LCR(m) = velmag_f(jj);
            wL_LCR(m) = w_lev(jj);
            m = m+1;
        end
    end
    
    for kk = 1:length(datenums_Waves)
   
        if abs(datenums_Waves(kk)-datenums_LC(ii))<ttolW 
    
           Wave_LCR(n) = H13(kk);
           n = n+1;
           
        end
        
    end
        
end    


%Normalize data
stan_LC = (LC_filt-min(LC_filt))/(max(LC_filt)-min(LC_filt));
stan_WL = (wL_LCR-min(wL_LCR))/(max(wL_LCR)-min(wL_LCR));
stan_curmag = (current_LCR-min(current_LCR))/(max(current_LCR)-min(current_LCR));
stan_H13 = (Wave_LCR-min(Wave_LCR))/(max(Wave_LCR)-min(Wave_LCR));

RAO_tens_cur = stan_LC./stan_curmag;
RAO_tens_WL = stan_LC./stan_WL;
RAO_tens_Wave = stan_LC([1:28 30:end])./stan_H13;

%% Look at Tension Changes for different Current Values
  
% figure;
% dkk = 1;
% for curkk=2:length(CM_edges);
% 
%     LC_cukk = LC_filt(:,logical(cdn_CM(:,dkk).*cdn_all_RAO));
%     
%     for mkk=3%mkks
%         hold on
%         plot(datenums_LC(logical(cdn_CM(:,dkk).*cdn_all_RAO)),LC_cukk,'DisplayName',[num2str(CM_edges(curkk-1)) ' < Cur Mag < ' num2str(CM_edges(curkk)) ' m/s'])
%         % clrkk=clrkk+1; %counter to assign colors to plot
%     end
% %     pause
% dkk = dkk+1;
% end
% 
% datetickzoom
% %ylim([0 1.1])
% title('Mean Tensions for Different Current Speeds')
% ylabel('Tension (lbs)')
% legend(gca,'show')
% set(gca, 'ColorOrderIndex', 1)





 

