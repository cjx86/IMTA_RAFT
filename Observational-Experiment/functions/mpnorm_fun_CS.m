%Converted to function by Toby Dewhurst, 2014.

% EDITED BY COREY SULLIVAN 2016

% -------------------------------------------------------
%       UNH OOA - Motion Pack 
%
%       112 bytes of header
%				48 bytes of title
%				9 bytes of date
%				9 bytes of time
%				7 bytes of battery
%				38 bytes of free space filled with zeros
%				1 byte of fixed end "1A" hex
%       12 A/D channels sampled at 10 hz
%				1. Motion Pack Acceleration - z axis
%				2. Motion Pack Acceleration - y axis
%               3. Motion Pack Acceleration - x axis
%               4. Motion PAck Rotation - z axis
%				5. Motion Pack Rotation - y axis
%				6. Motion Pack Rotation - x axis
%				7. Optical Backscattering Sensor
%				8. Load Cell #8
%				9. Load Cell #9
%			  10. Battery #1
%			  11. Battery #2
%			  12. Scanning Sonar Signal
%
% fname = input('What is the filename  ','s');

%Ang in degrees, but
%Sj1Side(:,4:6) in rad(^2)/Hz


function [date, time, conf_cut, Disp20, Ang20, Disp_nf, Ang_nf, t, batt, a_mean20,a_mean_nf, rot_mean20, bandfreqs, Sj1side20_cut, Sj1side30_cut, Sj1side_nf_cut, bands, bandwidth, Sjarea20, tsvar20, Sjarea30, tsvar30, Sjarea_nf, tsvar_nf]=mpnorm_fun_CS(fname,ensembles,padyes,hannme,stds,passes,mean_bands,alph)

% fname='..\Data\Surfaced\14110216.00'
% fname='..\Data\Surfaced\14102214.00'

freq=10; %Sampling frequency

    %Open the file that was input
fid = fopen(fname,'r');
    %Read the header of the file to get information
hdr = fread(fid,112)';
hdr = setstr(hdr);
filetitle=hdr(1:48);
    %Date and time info
date=hdr(49:57);
time=hdr(58:66);
    %Battery log
batt=hdr(67:73);
filll=hdr(74:111);
finn=hdr(112);
%
% np=2739.*input('Enter number of blocks (65k pieces) in the file ');
% 4 blocks of data at 2730 points per block on each of 12 channels
np=4*2730;

% Read the actual data from the IMU (unit8 means unassigned 8 bit integers)
[dat,count] = fread(fid,2*12*np,'uint8');
fclose(fid);
%
rdat=[256*dat(1:2:length(dat)) + dat(2:2:length(dat))];
for i=1:length(rdat)
   if rdat(i) > 32768, rdat(i)=rdat(i)-65536 ; end
end
%
% digitizer is set to plus and minus 10 volts, or 20 volts over 12 bits
% therefore, volts output is 0.00488 volts/bit * counts

ndat=0.00488.*reshape(rdat,12,np);
%Truncate (at 17 minutes to elimate odd dip/spike/drift at end)
ndat=ndat(:,1:10200);

% calculate x-axis time base
[n,tim]=size(ndat);
if n~=12, disp('error in size'), end
t_2=0:(tim-1); %correcting the total time to start form 0
t=t_2./freq; %The time in seconds of each sample

% CONVERTING IMU Output into MEANINGFUL DATA FROM LAB TESTING OF OUTPUT

    % convert accelerations to g's
% factory calibration is 3.742 volts/g in z axis with 0.00083 g bias
zacc=ndat(1,:)./3.742 - 0.00083;
% factory calibration is 7.458 volts/g in y zxis with -0.00176 g bias
yacc=ndat(2,:)./7.485 + 0.00176;
% factory calibration is 7.504 volts/g in x axis with 0.00033 g bias
xacc=ndat(3,:)./7.504 - 0.00033;
%
% convert rates to deg/sec
% factory calibration is 0.050304 volts per deg/sec rate
% with 0.06 deg/s bias
% lab tests show offset of 0.1523 deg/s
zrate= ndat(4,:)./0.050304-0.06-0.1523;
% factory calibration is 0.049933 volts per deg/sec rate
% with 0.14 deg/s bias
% lab tests show offset of 0.1508 deg/s
yrate=ndat(5,:)./0.049933-0.14-0.1508;
% factory calibration is 0.050101 volts per deg/sec rate with 0.06 deg/s bias
% lab tests show offset of 0.1348 deg/s
xrate=ndat(6,:)./0.050101-0.06-0.1348 ;
%
% load cells
   lc1=ndat(7,:);
   lc2=ndat(8,:);
%
% optical bakscattering sensor
    obs=ndat(9,:);
%
% battery or bus voltage after diodes
    bat=1.565.*ndat(10,:);

% plot data IF DESIRED
% figure
% subplot(611),plot(t,zacc),grid
% title(['Mussel Raft Motion Package ',fname])
% ylabel('Z Acc')
% subplot(612),plot(t,yacc),grid
% ylabel('Y Acc')
% subplot(613),plot(t,xacc),grid
% ylabel('X Acc')
% subplot(614),plot(t,zrate),grid
% ylabel('Z Rot')
% subplot(615),plot(t,yrate),grid
% ylabel('Y Rot')
% subplot(616),plot(t,xrate),grid
% ylabel('X Rot')
% xlabel('Time in seconds')
%

% integrate acceleration to displacement
 [xdisp20,xvel20,a_xmean20]=accel2motion20(xacc,10,0,0);
 [ydisp20,yvel20,a_ymean20]=accel2motion20(yacc,10,0,0);
 [zdisp20,zvel20,a_zmean20]=accel2motion20(zacc,10,0,0);
 
 [xdisp30,xvel30,a_xmean30]=accel2motion30(xacc,10,0,0);
 [ydisp30,yvel30,a_ymean30]=accel2motion30(yacc,10,0,0);
 [zdisp30,zvel30,a_zmean30]=accel2motion30(zacc,10,0,0);
 
 [xdisp_nf,xvel_nf,a_xmean_nf]=accel2motion(xacc,10,0,0);
 [ydisp_nf,yvel_nf,a_ymean_nf]=accel2motion(yacc,10,0,0);
 [zdisp_nf,zvel_nf,a_zmean_nf]=accel2motion(zacc,10,0,0);

% subplot(311),plot(t,xdisp20),grid
% title(['Mussel Raft Rotation ',fname])
% ylabel('X Displacement - m')
 
% subplot(312),plot(t,ydisp20),grid
% ylabel('Y Displacement - m')

% subplot(313),plot(t,zdisp20),grid
% ylabel('Z Displacement - m')
% xlabel('Time in seconds')

% Integrate rate of rotation to an angle

[xang20, dgps_meanx20]=rate2angle20(xrate,10,0);
[yang20, dgps_meany20]=rate2angle20(yrate,10,0);
[zang20, dgps_meanz20]=rate2angle20(zrate,10,0);
 
[xang30, dgps_meanx30]=rate2angle30(xrate,10,0);
[yang30, dgps_meany30]=rate2angle30(yrate,10,0);
[zang30, dgps_meanz30]=rate2angle30(zrate,10,0);
 
[xang_nf, dgps_meanx_nf]=rate2angle(xrate,10,0);
[yang_nf, dgps_meany_nf]=rate2angle(yrate,10,0);
[zang_nf, dgps_meanz_nf]=rate2angle(zrate,10,0);


Disp20=[xdisp20 ydisp20 zdisp20];
Ang20=[xang20 yang20 zang20];

a_mean20=[a_xmean20 a_ymean20 a_zmean20];
rot_mean20=[dgps_meanx20 dgps_meany20 dgps_meanz20];
 
Disp30=[xdisp30 ydisp30 zdisp30];
Ang30=[xang30 yang30 zang30];
a_mean30=[a_xmean30 a_ymean30 a_zmean30];
rot_mean30=[dgps_meanx30 dgps_meany30 dgps_meanz30];

%zdisp_nf = zdisp_nf - L_x.*sind(yang20) - L_y.*sind(xang20);

Disp_nf=[xdisp_nf ydisp_nf zdisp_nf];
Ang_nf=[xang_nf yang_nf zang_nf];
a_mean_nf=[a_xmean_nf a_ymean_nf a_zmean_nf];
rot_mean_nf=[dgps_meanx_nf dgps_meany_nf dgps_meanz_nf];

% subplot(311),plot(t,xang),grid
% title(['Mussel Raft Rotation ',fname])
% ylabel('X Rotation in deg')
% subplot(312),plot(t,yang),grid
% ylabel('Y Rotation in deg')
% subplot(313),plot(t,zang),grid
% ylabel('Z Rotation in deg')
% xlabel('Time in seconds')

%% Calculate spectra

%Displacement
Ffreq=freq/(floor(length(ndat)/ensembles)) %Fourier frequency spacing
bands=round(0.5/Ffreq/64) %Band average to match wave data (64 freqs between 0 and 0.5 Hz)
    %This is how many bands you need to average to take all the info from
    %the IMU file and put it into 64 bands

for mkk=1:3
% [IMU_bandfreqs(:,mkk) IMU_Sj1side(:,mkk) IMU_con95(:,mkk) IMU_bandwidth(:,mkk) IMU_Sjarea(:,mkk) IMU_var(:,mkk)]=spectra(q(:,mkk),IMU_freq,ensembles,bands,padyes,hannme,stds,passes,mean_bands(mkk));
[~        , Sj1side20(1,:,mkk), con95(:,mkk), bandwidth, Sjarea20(1,mkk), tsvar20(1,mkk)]=spectra2(Disp20(:,mkk),freq,ensembles,bands,padyes,hannme,stds,passes,mean_bands,alph);
[bandfreqs, Sj1side20(1,:,3+mkk), con95(:,3+mkk), bandwidth, Sjarea20(1,3+mkk), tsvar20(1,3+mkk)]=spectra2(Ang20(:,mkk)*pi/180,freq,ensembles,bands,padyes,hannme,stds,passes,mean_bands,alph);

[~        , Sj1side30(1,:,mkk), con95(:,mkk), bandwidth, Sjarea30(1,mkk), tsvar30(1,mkk)]=spectra2(Disp30(:,mkk),freq,ensembles,bands,padyes,hannme,stds,passes,mean_bands,alph);
[~        , Sj1side30(1,:,3+mkk), con95(:,3+mkk), bandwidth, Sjarea30(1,3+mkk), tsvar30(1,3+mkk)]=spectra2(Ang30(:,mkk)*pi/180,freq,ensembles,bands,padyes,hannme,stds,passes,mean_bands,alph);

[~        , Sj1side_nf(1,:,mkk), ~, ~, Sjarea_nf(1,mkk), tsvar_nf(1,mkk)]=spectra2(Disp_nf(:,mkk),freq,ensembles,bands,padyes,hannme,stds,passes,mean_bands,alph);
[~        , Sj1side_nf(1,:,3+mkk), ~, ~, Sjarea_nf(1,3+mkk), tsvar_nf(1,3+mkk)]=spectra2(Ang_nf(:,mkk)*pi/180,freq,ensembles,bands,padyes,hannme,stds,passes,mean_bands,alph);

%[~        , Sjacc(1,:,mkk), ~, bandacc(:,mkk), ~, ~]=spectra2(amean20(:,mkk),freq,ensembles,bands,padyes,hannme,stds,passes,mean_bands,alph);

end

cutoff=0.5; %Hz. Truncate spectrum at 0.5 Hz.
cutkk=find(bandfreqs>cutoff,1);
bandfreqs=bandfreqs(1:cutkk);
Sj1side20_cut(1,:,:)=Sj1side20(1,1:cutkk,:);
Sj1side30_cut(1,:,:)=Sj1side30(1,1:cutkk,:);
Sj1side_nf_cut(1,:,:)=Sj1side_nf(1,1:cutkk,:);
conf_cut=con95(1:cutkk,:);

end

%% Plot
% semilogy(bandfreqs(2:length(bandfreqs)),Sj1side(2:length(Sj1side)),'cy')
% figure
% loglog(bandfreqs(2:end), squeeze(Sj1side_cut(1,2:end,1)))
% hold all
% loglog(bandfreqs(2:end), squeeze(Sj1side_nf_cut(1,2:end,1)))
% % 
% xrange_s=1;%0.5*freq; 
% yrange_s=1.2*max(Sj1side);
% title('One-sided Sample Spectral Density of Heave Displacement')
% xlabel('Frequency')
% ylabel('Spectral Density')
% clear param
% param(1)={['Ensembles: ' num2str(ensembles) ' ']};
% param(2)={['Bands: ' num2str(bands) ' ']};
% text(xrange_s-1,yrange_s-0.1,param,'HorizontalAlignment','right','VerticalAlignment','top')
% xlim([0 xrange_s])
% % ylim([0 yrange_s])
% %plot confidence
% hold all
% line([0.8*xrange_s 0.8*xrange_s],[0.3*yrange_s 0.3*yrange_s+con95],'Color',[0 1 0])
% line([0.8*xrange_s 0.8*xrange_s+bandwidth],[0.2*yrange_s 0.2*yrange_s],'Color',[0 1 0])
% legend('Spectrum', '95% Confidence','Bandwidth','Location','NorthEast')



%  %% Test by comparing with fft routine
%  Fs=freq;
% Nsamps = length(y);
% t = (1/Fs)*(1:Nsamps);          %Prepare time data for plot
% y_fft = abs(fft(y));            %Retain Magnitude
% y_fft = y_fft(1:Nsamps/2);      %Discard Half of Points
% f = Fs*(0:Nsamps/2-1)/Nsamps;   %Prepare freq data for plot

