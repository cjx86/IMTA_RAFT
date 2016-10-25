
clc, clear, close all
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))
addpath(genpath('../'));
addpath(genpath('../functions'));

%% Current information

load('IMTA Raft Study 20160113T152340_averaged.mat');

%Convert output time to mm/dd
unix_epoch = datenum(1970,1,1,0,0,0);
matlab_time = sens.time./86400 + unix_epoch; 
%test_t = sens.dnum;
date_time = matlab_time(57:21667);  %Truncate data to eliminate deployment/recovery

%Truncate depth data to eliminate deployment/recovery
depth = sens.pd(57:21667);

% Pull out the velocity information

vel_E = wt.vel(57:21667,1:19,1)';
vel_N = wt.vel(57:21667,1:19,2)';
vel_up = wt.vel(57:21667,1:19,3)';
vel_Err = wt.vel(57:21667,1:19,4)';

%depth ranges in m where to vertically average
    %Start avg .5 m below the surface
ind_T=floor(((depth-.5)-info.cell1)/info.cell); %Index for top to avg
    %Mussel Ropes are 13' Net is 12' set near one of those
ind_B=ceil(((depth-15*0.3048)-info.cell1)/info.cell); %Index for bott

% %For AVGING over constant depths, ignoring distance from surface
%  d_L = 2.67; 
%  d_H = 6.67;
% %     %Find the indices matching the limits
%  ind_B = find(wt.r(1,:) == d_L);  
%  ind_T = find(wt.r(1,:) == d_H);  

    
    %Shorten the t.s. 
v_E = vel_E(ind_B:ind_T,:);
v_N = vel_N(ind_B:ind_T,:);

    %Find magnitude and direction from full N,E t.s.
v_mag = sqrt(v_E.^2 + v_N.^2);
v_dir = atan2(v_E,v_N);


bot_weight = repmat(max(v_mag,[],1),8,1);

vmag_wt = v_mag./bot_weight;
v_dir_weighted = v_dir.*vmag_wt;

decl = -15.31667; %declination in degrees

v_dir = v_dir + deg2rad(decl);
v_dir_weighted = v_dir_weighted + deg2rad(decl);

    %Calculate avg mag and dir
velmag_avg = nanmean(v_mag(1:7,:),1);
vel_dir_avg = nanmean(v_dir(1:7,:),1);
vel_dir_Wavg = nanmean(v_dir_weighted,1);
vel_dir_deg = 180/pi*vel_dir_avg;
vel_dir_Wdeg = 180/pi*vel_dir_Wavg;

freq = 1;
T_lp = 120;    

    %Vertically average N,E over the desired space first
velN_avg = mean(v_N(1:7,:),1);
%[ates,velN_comp,~,~,~,~] = stddev_filter(date_time,velN_avg,3,3);
velE_avg = mean(v_E(1:7,:),1);
%[btes,velE_comp,~,~,~,~] = stddev_filter(date_time,velE_avg,3,3);

    %Calculate mag and dir from avged N and E components
velmag2 = velN_avg.^2 + velE_avg.^2;
vel_dir2 = atan2(velE_avg,velN_avg);
vel_dir_deg2 = 180/pi*vel_dir2;

% figure,plot(date_time,depth),grid,datetick
% title('Depth (m)')

    %These are smoothed sets of the averaged data

velmag_sg = sgolayfilt(velmag_avg,5,T_lp*freq+1); 
vmag2_sg = sgolayfilt(velmag2,5,T_lp*freq+1); 

veldir_sg=sgolayfilt(vel_dir_avg,5,T_lp*freq+1); 
veldir_sg_deg = sgolayfilt(vel_dir_deg,5,T_lp*freq+1);
veldir_Wsg_deg = sgolayfilt(vel_dir_Wdeg,5,T_lp*freq+1);

veldir2_sg = sgolayfilt(vel_dir_deg2,5,T_lp*freq+1); 

figure,plot(veldir_sg_deg),hold on,plot(veldir_Wsg_deg)

%% Plotting 

figure,plot(date_time(3725:13805),velmag_sg(3725:13805),'LineWidth',1.5),grid
datetickzoom('x','mm/dd','keepticks')
ylabel('Velocity Magnitude, [m/s]')
xlabel('Date, [mm/dd]')
set(gca,'FontSize',22)
title('Depth Averaged Current Magnitude During Deployment','FontSize',28)
xlim([date_time(3725) date_time(13805)])
ylim([0 0.675])

figure,plot(date_time(3725:13805),veldir_sg_deg(3725:13805),'Color',[.85 .33 .1],'Linewidth',1.5),grid
datetickzoom('x','mm/dd','keepticks')
set(gca,'FontSize',22)
title('Depth Averaged Current Direction During Deployment','FontSize',28)
xlabel('Date, [mm/dd]')
ylabel('Direction, [deg]')
xlim([date_time(3724) date_time(13806)])
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)

    % Velocity magnitude
figure,plot(date_time(3725:13805),velmag_sg(3725:13805),'LineWidth',1.5),grid
datetickzoom('x','mm/dd','keepticks')
set(gca,'FontSize',16)
ylim([0 .625])
title('Depth Averaged Velocity Magnitude During Deployment','FontSize',18)
xlabel('Date [mm/dd]')
ylabel('Velocity Magnitude [m/s]')
xlim([date_time(3724) date_time(13806)])
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)

    % Velocity Direction
figure,plot(date_time,veldir_sg_deg,date_time,veldir2_sg),grid
datetickzoom('x','mm/dd HH','keepticks')
title('Current Direction During Deployment as a Function of Depth (m)')
xlabel('Date (mm/dd)')
ylabel('Direction in deg')

figure
subplot(2,1,1),plot(date_time,velmag_sg,'LineWidth',1.5),grid
datetickzoom('x','mm/dd','keepticks')
set(gca,'FontSize',16)
xlabel('Date [mm/dd]')
ylabel('Velocity Magnitude [m/s]')
subplot(2,1,2),plot(date_time,veldir_sg_deg,'r'),grid
datetickzoom('x','mm/dd','keepticks')
ylabel('Current Direciton in Deg') % right y-axis
suptitle('Showing Current Direction and Magnitude for alignmment')
xlabel('Date (mm/dd)')

%% Confidence intervals
% alph=0.05;  % Setting desired confidence
% DOFs_adcp=2*bands_IMU*4; % Setting DOF's
% [chi2l, chi2h]=chi2_vals(alph,DOFs_IMU);    %Automated way to find CHI^2
%     %Error bar creation for entire plot
% Sl_IMU=DOFs_IMU*Sj1side20/chi2l; %Needs to be a function of Sj magnitude if not on log plot
% Sh_IMU=DOFs_IMU*Sj1side20/chi2h; %Ends of error bars
%     
%     %Actual sizing of error bars to plot
% errL=Sj1side20-Sl_IMU; %Length of error bars
% errH=Sh_IMU-Sj1side20; %Length of error bars


%% Wave Analysis portion

spectra = dlmread('spectrum_ts.txt');
data = dlmread('wave_output.txt',',');

f_step = .5/length(spectra(1,:));
spec_freqs = f_step:f_step:0.5;

%Pull the date out of the data file
yr = data(:,1);
mo = data(:,2);
day = data(:,3);
hr = data(:,4);
min = data(:,5);
sec = data(:,6);

datenums=datenum(2000+yr,mo,day,hr,min,00);
ttol=datenum(0,0,0,0,1,0); %Within a minute

stkk=find(abs(datenums-datenum(2016,1,15,16,00,00))<ttol,1); %Start time dictated by IMU
endkk=find(abs(datenums-datenum(2016,1,19,10,00,00))<ttol,1); %end time
endkk = endkk-1; %Placing the last piece at 1/19 14 hr 20 min

s_kk=find(abs(date_time-datenum(2016,1,15,16,00,00))<ttol,1)+1; %Start time dictated by IMU
e_kk=find(abs(date_time-datenum(2016,1,19,10,00,00))<ttol,1)+1; %end time

%Wave Height information
Hs = data(:,7);
H10 = data(:,8);
H13 = data(:,9);
Hmax = data(:,10);
%Wave Period Information
Tp = data(:,11);
T10 = data(:,12);
T13 = data(:,13);
Tmax = data(:,14);

Dp = zeros(length(data(:,15)),1);

for ii = 1:length(data(:,15))
    
   if data(ii,15) == 0
        
        Dp(ii) = 0;
    
    else
        Dp(ii) = data(ii,15)+decl;

    end
end

Depth = data(:,16);
Cur_mag = data(:,17);
Cur_dir = data(:,18)+decl;

figure,
plot(datenums,H13)
ylabel('Significant Wave Height, m')
datetick('x','keepticks')

% sname = '..\Combined\Data\First\Full_ADCP_data.mat';
% 
% save(sname,'data','spectra','datenums','date_time','velmag_sg','veldir_sg_deg','depth')

%Trim Current data to align with other instruments
datenums_cur = date_time(s_kk:e_kk);
velmag_f=velmag_sg(s_kk:e_kk);
veldir_f=veldir_sg_deg(s_kk:e_kk);

H13_nmask = ones(length(H13),1);
Dp_nmask = ones(length(Dp),1);
Tp_nmask = ones(length(Tp),1);

for i = 1:length(H13)
    
    if H13(i) == 0
        H13_nmask(i) = 0;
    end
    
    if Dp(i) == decl-1
        Dp_nmask(i) = 0;
    end
    
    if Tp(i) == -0.1
        Tp_nmask(i) = 0;
    end
    
end

figure,
plot(datenums(logical(H13_nmask)),H13(logical(H13_nmask)),'LineWidth',1.5)
ylabel('Significant Wave Height [m]')
datetickzoom('x','mm/dd','keepticks')
xlabel('Date [mm/dd]')
set(gca,'XLim',[datenums(1)-.05 datenums(end)])
set(gca,'FontSize',20)
title('Significant Wave Height','FontSize',24)
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)

Dp_tot = Dp;
Dp_mtot = Dp(logical(Dp_nmask));
date_dp_tot = datenums(logical(Dp_nmask));

figure,
plot(date_dp_tot(1:287),Dp_mtot(1:287),'LineWidth',1.5)
hold on,plot(date_dp_tot(288:end),Dp_mtot(288:end),'Color',[0 .45 .74],'LineWidth',1.5)
ylabel('Dominant Wave Direction, [deg]')
xlim([datenums(stkk) datenums(end)+0.01])
datetickzoom('x','mm/dd','keepticks')
xlabel('Date [mm/dd]')
set(gca,'FontSize',20)
title('Dominant Wave Direction During ADCP Deployment','FontSize',24)
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)

Tp_tot = Tp(logical(Tp_nmask));
date_tp_tot = datenums(logical(Tp_nmask));

figure,
plot(date_tp_tot(stkk:289),Tp_tot(stkk:289),'LineWidth',1.5)
hold on,plot(date_tp_tot(290:end),Tp_tot(290:end),'Color',[0 .45 .74],'LineWidth',1.5)
ylabel('Peak Wave Period [sec]')
xlim([datenums(stkk)+.5 datenums(end)+0.01])
datetickzoom('x','mm/dd','keepticks')
xlabel('Date [mm/dd]')
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
set(gca,'FontSize',20)
title('Peak Wave Period During ADCP Deployment','FontSize',24)
hold off

%Trim Waves data to align with other instruments
datenums_Waves=datenums(stkk:endkk);
Cur_dir=Cur_dir(stkk:endkk);
Cur_mag=Cur_mag(stkk:endkk);
H13=H13(stkk:endkk);
T13=T13(stkk:endkk);
Hs=Hs(stkk:endkk);
Tp=Tp(stkk:endkk);
Dp=Dp(stkk:endkk);
last_day_hr = data(stkk:endkk,3:6);
Waves_Spectra=spectra(stkk:endkk,:); %mm/sqrt(Hz)
w_lev = depth(s_kk:e_kk);

dropouts=isnan(datenums);

%Save all the desired data
% savename='..\Combined\Data\First\Waves_meta.mat';
% save(savename,'datenums_Waves','datenums_cur','dropouts','Cur_dir','veldir_f','Cur_mag','velmag_f','H13','T13','Hs','Tp','Dp','Waves_Spectra','w_lev')

figure,
plot(datenums_Waves(1:end-30),H13(1:end-30),'LineWidth',1.5)
title('Significant Wave Height During Deployment of ADCP','FontSize',22)
ylabel('Significant Wave Height [m]')
datetickzoom('x','mm/dd HH','keepticks')
xlabel('Date [mm/dd HH]')
set(gca,'XLim',[datenums_Waves(1)-.05 datenums_Waves(end-30)+.05])
set(gca,'FontSize',18)
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)

Dp_mask = ones(length(Dp)-50,1);
for i = 1:length(Dp)-50    

    if Dp(i) == decl-1
        Dp_mask(i) = 0;
    end
    
end

Dp2 = Dp(logical(Dp_mask));

figure,
plot(datenums_Waves(logical(Dp_mask)),Dp2,'LineWidth',1.5)
ylabel('Dominant Wave Direction [deg]')
datetickzoom('x','mm/dd HH','keepticks')
xlabel('Date [mm/dd HH]')
set(gca,'FontSize',20)
set(gca,'XLim',[datenums_Waves(1)-.05 datenums_Waves(end-50)+.05])
title('Dominant Wave Direction During Portion of ADCP Deployment','FontSize',24)
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)


figure,
plot(datenums_Waves(10:end-54),(Tp(10:end-54)+0.1),'LineWidth',1.5);
datetickzoom('x','mm/dd HH','keepticks')
xlabel('Date [mm/dd HH]')
ylabel('Peak Wave Period [sec]')
set(gca,'FontSize',20)
title('Peak Wave Period Height During Portion of ADCP Deployment','FontSize',24)
set(gca,'XLim',[datenums_Waves(6) datenums_Waves(end-50)])
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
q.LineWidth = 1.5;


figure,
subplot(2,1,1),plot(datenums_Waves,H13,'LineWidth',1.5)
ylabel('Significant Wave Height, [m]')
datetickzoom('x','mm/dd')
set(gca,'FontSize',20)
title('Significant Wave Height','FontSize',28)
xlim([datenums_Waves(1)-0.01 datenums_Waves(end)+0.01])
grid on
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
subplot(2,1,2),plot(datenums_Waves,Dp,'r','LineWidth',1.5)
ylabel('Dominant Wave Direction, [deg]')
datetickzoom('x','mm/dd')
xlim([datenums_Waves(1)-0.01 datenums_Waves(end)+0.01])
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
set(gca,'FontSize',20)
ylim([-25 360])
title('Dominant Period','FontSize',28)

% figure
% plot(datenums_Waves,Cur_mag)
% legend('Average Mag only','Avg N,E then Mag','Waves Output')
% 
% figure
% plot(datenums_Waves,wrapTo180(Cur_dir))
% legend('Average Dir only','Avg N,E first','Waves Output')


figure,
semilogy(spec_freqs,1/10^6*Waves_Spectra(39,:).^2,'Linewidth',1.5)
ylabel('Spectral Energy Density,  [m^2/Hz]')
xlabel('Frequency, [Hz]')
set(gca,'FontSize',20)
title(['Wave Energy Spectral Density for ' datestr(datenums_Waves(39),'mm/dd/yy HH:MM')],'FontSize',24)

figure,
semilogy(spec_freqs(2:end-1),(1/10^6.*mean(Waves_Spectra(:,2:end-1),1).^2),'Linewidth',1.5);
xlabel('Frequency, [Hz]')
ylabel('Spectral Energy Density, [m^2/Hz]')
set(gca,'FontSize',20)
title('Average of All Wave Energy Spectral Density','FontSize',24)
%Minorticks on
    %Compass Rose to show flow direction
    
L_ens = 2883; %Length of ts for 2 days   

j = 1;
for i = 1:length(veldir_sg)
    
    if velmag_sg(i) > .2
    
        if (0 <= veldir_sg(i) && veldir_sg(i) <= pi/12)
    
            flood_dir1(j) = veldir_sg(i);
            j = j+1;
            
        end
        
        if (355*pi/180) <= veldir_sg(i) &&  veldir_sg(i) < 2*pi
            
            flood_dir1(j) = veldir_sg(i);
            j = j+1;
        end
        
    end
end

flood_dir = flood_dir1;

for ii = 1:length(flood_dir1)
    
    if flood_dir1(ii) >= (355*pi/180)
        
        flood_dir(ii) = flood_dir1(ii)-2*pi;
        
    end
    
end

flood_dir_mn = mean(flood_dir)*180/pi;
    
figure,
a = polar(veldir_sg,velmag_sg,'.');
set(a,'MarkerSize',3)
th = findall(gcf,'Type','text');

for i = 1:length(th),
    set(th(i),'FontSize',22)
end

moorax = 108.660278;

graft = hgtransform;
raft = rectangle('position',0.05*[-1,-1.9,2,3.72],'Parent',graft);
graft.Matrix = makehgtform('zrotate',moorax*pi/180);
view([90 -90])
suptitle('Compass Plot of Current Direction')

dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@polarangs)

