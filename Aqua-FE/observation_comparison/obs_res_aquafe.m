clc,clear,close all

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))
addpath(genpath('..\..\..\Field_Testing'));
addpath(genpath('..\..\..\Field_Testing\TJD_MAT_funs'));
addpath(genpath('..\..\'));

pathdir = pwd;
dirData = dir(pathdir);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
fileList = {dirData(~dirIndex).name}';  %# Get a list of the files

    ids = 1:length(fileList);   
    
    %Area of  element
A_LC = (2.027*10^-3); %m^2 from .opt file

%% Currents to pretentison to 250 lbs

bear = 108.660278; %degrees

v_mag_pre = 0.60587*8/25;
v_dir_p = 5.9; %Visually found from ADCP data
v_dir_modp = 90 - (bear-v_dir_p);

%% General runs of single freq waves

D_sing = 8.9579;

dom_wav_dir = 135;

longside_raft_ang = 180-(dom_wav_dir-(bear - 90));
    %This may be slightly off by a couple degrees because I am using the 

ang_of_cur_modelg = longside_raft_ang - (bear - 90) + v_dir_p;

vconfz_p = v_mag_pre.*cosd(ang_of_cur_modelg);
vconfx_p = v_mag_pre.*sind(ang_of_cur_modelg);

T_wvs = [2 2.25 2.5 3 10/3 11/3 4:1:6 20/3 7 8 9 10 11 12 15 18 20 24]; %seconds

%T_wvs = [2 2.25 2.5 3 10/3 11/3 4:1:6 20/3 7 8 9 9 10 10 11 11 12 15 18 20 24]; %seconds

NT_wvs = [9 10 11];

L_sing = zeros(1,length(T_wvs));

g = 9.81;

dp = zeros(length(T_wvs),1);

%Calculating Wavelength
for i = 1:length(T_wvs)

    [L_sing(i),~] = find_L_dispersion(T_wvs(i),D_sing);

    if L_sing(i) < 2*D_sing %Check for deep water
        k(i) = 2*pi./L_sing(i);
                    
    else
        k(i)=(2*pi./T_wvs(i)).^2/g; %Deep water dispersion relation
        dp(i) = 1;
    end
    
end

%From 3 2/3 sec up is DEEP WATER

H_wvs = 0.2;
a = H_wvs/2;   

%FIRST RUN 1.8112 SEC WAVES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mid Water Waves 1/3 Max Cur

RAO_cog = zeros(1,length(k));
RAO_IMU = zeros(1,length(k));

for id=1:length(T_wvs)

   fil(id,:) = [fileList{id}];

mid_data(id,:,:) = importdata(fil(id,:));
pd1(id,:) = fil(id,13:14);
pd2(id,:) = fil(id,16);

t_wc(id,:) = mid_data(id,2:end-1,1);
mid_stress_shore(id,:) = mid_data(id,2:end-1,2);
mid_stress_chan(id,:) = mid_data(id,2:end-1,3);

    %Convert Pa to N and multiply by C.S. area (in inches) to get lbs  

mid_load_Nsh(id,:) = mid_stress_shore(id,:).*A_LC;
mid_load_Nch(id,:) = mid_stress_chan(id,:).*A_LC;
mid_load_lbsh(id,:) = mid_load_Nsh(id,:)*.22481;
mid_load_lbch(id,:) = mid_load_Nch(id,:)*.22481;

cog(id,:) = mean(mid_data(id,2:end-1,4:7),3);
IMU(id,:) = mid_data(id,2:end-1,8);

rol_dif(id,:) = mid_data(id,2:end-1,5)-mid_data(id,2:end-1,7);
rol_dif2(id,:) = mid_data(id,2:end-1,4)-mid_data(id,2:end-1,6);


roll1(id,:) = asin(rol_dif(id,:)./1.2192);
roll2(id,:) = asin(rol_dif(id,:)./1.2192);

roll(id,:) = (roll1(id,30:end) + roll2(id,30:end))/2;


pit_dif1(id,:) = mid_data(id,2:end-1,6)-mid_data(id,2:end-1,7);
pit_dif2(id,:) = mid_data(id,2:end-1,4)-mid_data(id,2:end-1,5);

pitch1(id,:) = asin(pit_dif1(id,:)./.2032);
pitch2(id,:) = asin(pit_dif2(id,:)./.2032);

pitch(id,:) = (pitch1(id,30:end) + pitch2(id,30:end))/2;


%% Shortern series for spectra analysis
cg1(id,:) = cog(id,30:end)-mean(cog(id,30:end));
IU1(id,:) = IMU(id,30:end)-mean(IMU(id,30:end));
sig(id) = 2*pi./T_wvs(id);


wave_load1(id,:) = a*cos(k(id)-sig(id).*t_wc(id,30:end-1));
wave_steep1(id,:) = a.*k(id);

%Avg Heave amplitude
[avgcog(id)] = amps_CS(t_wc(id,30:end),cg1(id,:),1,length(cg1(1,:)),T_wvs(id));
[avgIMU(id)] = amps_CS(t_wc(id,30:end),IU1(id,:),1,length(IU1(1,:)),T_wvs(id));    

%Heave RAOs    
RAO_cog(id) = avgcog(id)/a;
RAO_IMU(id) = avgIMU(id)/a;

%Avg Pitch and roll amplitudes
[avgp(id)] = amps_CS(t_wc(id,30:end),pitch(id,:),1,length(pitch(1,:)),T_wvs(id));
[avgr(id)] = amps_CS(t_wc(id,30:end),roll(id,:),1,length(roll(1,:)),T_wvs(id));

%pitch Roll RAOs
RAO_pitch(id) = avgp(id)./wave_steep1(id);
RAO_roll(id) = avgr(id)./wave_steep1(id);

%Dimensional Pitch Roll RAOs
RAO_pitch_dim(id) = avgp(id)./a;
RAO_roll_dim(id) = avgr(id)./a;


end

figure,plot(cg1(20,:)),hold on,plot(200,avgcog(20),'o'),plot(200,-avgcog(20),'o'),hold off
% Motion Check


for ii = 1:length(T_wvs)
    
    if ii <= 11

[avg_loadsh(ii)] = amps_CS(t_wc(ii,:),mid_load_Nsh(ii,200:end),1,length(mid_load_Nsh(1,200:end)),T_wvs(ii));
[avg_loadch(ii)] = amps_CS(t_wc(ii,:),mid_load_Nch(ii,200:end),1,length(mid_load_Nch(1,200:end)),T_wvs(ii));      

avg_loadsh(ii) = avg_loadsh(ii) + mean(mid_load_Nsh(ii,200:end));
avg_loadch(ii) = avg_loadch(ii) + mean(mid_load_Nch(ii,200:end));

mn_loadsh(ii) = mean(mid_load_Nsh(ii,200:end));
mn_loadch(ii) = mean(mid_load_Nch(ii,200:end));

std_loadsh(ii) = std(mid_load_Nsh(ii,200:end));
std_loadch(ii) = std(mid_load_Nch(ii,200:end));

    else
[avg_loadsh(ii)] = amps_CS(t_wc(ii,:),mid_load_Nsh(ii,150:end),1,length(mid_load_Nsh(1,150:end)),T_wvs(ii));
[avg_loadch(ii)] = amps_CS(t_wc(ii,:),mid_load_Nch(ii,150:end),1,length(mid_load_Nch(1,150:end)),T_wvs(ii));

avg_loadsh(ii) = avg_loadsh(ii) + mean(mid_load_Nsh(ii,300:end));
avg_loadch(ii) = avg_loadch(ii) + mean(mid_load_Nch(ii,300:end));

mn_loadsh(ii) = mean(mid_load_Nsh(ii,300:end));
mn_loadch(ii) = mean(mid_load_Nch(ii,300:end));

std_loadsh(ii) = std(mid_load_Nsh(ii,300:end));
std_loadch(ii) = std(mid_load_Nch(ii,300:end));

  end
     
    %Load RAOsclose all
RAO_shore(ii) = avg_loadsh(ii)./a;
RAO_chan(ii) = avg_loadch(ii)./a;

%/2.20462 for lbs

end

%% Longer trials

for ii=1:3

    id = length(T_wvs)+ii;
    Nfil = [fileList{id}];

Nmid_data = importdata(Nfil);
Npd1 = Nfil(13:14);
Npd2 = Nfil(16);

Nt_wc = Nmid_data(2:end-1,1);
Nmid_stress_shore = Nmid_data(2:end-1,2);
Nmid_stress_chan = Nmid_data(2:end-1,3);

    %Convert Pa to N and multiply by C.S. area (in inches) to get lbs  

Nmid_load_Nsh = Nmid_stress_shore.*A_LC;
Nmid_load_Nch = Nmid_stress_chan.*A_LC;
Nmid_load_lbsh = Nmid_load_Nsh*.22481;
Nmid_load_lbch = Nmid_load_Nch*.22481;


[Navg_loadsh] = amps_CS(Nt_wc,Nmid_load_Nsh(200:end),1,length(Nmid_load_Nsh(200:end)),NT_wvs(ii));
[Navg_loadch] = amps_CS(Nt_wc,Nmid_load_Nch(200:end),1,length(Nmid_load_Nch(200:end)),NT_wvs(ii));      

Navg_loadsh = Navg_loadsh + mean(Nmid_load_Nsh(200:end));
Navg_loadch = Navg_loadch + mean(Nmid_load_Nch(200:end));

Nmn_loadsh(ii) = mean(Nmid_load_Nsh(200:end));
Nmn_loadch(ii) = mean(Nmid_load_Nch(200:end));

Nstd_loadsh(ii) = std(Nmid_load_Nsh(200:end));
Nstd_loadch(ii) = std(Nmid_load_Nch(200:end));

figure,plot(Nt_wc,Nmid_load_Nch,100,Navg_loadch,'o')

clear('Nmid_data','Nt_wc')
clear('Nmid_stress_shore','Nmid_stress_chan','Nmid_load_Nsh','Nmid_load_Nch')
%clear ('rol_dif','rol_dif2','roll1','roll2','roll')
%clear('pit_dif1','pit_dif2','pitch1','pitch2','pitch')
%clear ('cog','IMU','cg1','IU1','load_shN_P','load_chN_P')

    %Load RAOsclose all
NRAO_shore(ii) = Navg_loadsh./a;
NRAO_chan(ii) = Navg_loadch./a;

%/2.20462 for lbs

end


%% Plotting Data

figure,plot(t_wc(3,:),(mid_load_Nch(3,:)),'LineWidth',1.5),hold on
%plot(mn_loadch(9),'o')
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)

f_cg1 = 1./T_wvs;

f_cga = f_cg1;
f_IUa = f_cga;

%lc_sg=sgolayfilt((mid_load_Nch(8,:)),1,T_wvs(8)*5); %Assuming 5 Hz. Number of frames, f must be odd.

figure,plot(t_wc(8,:),(mid_load_Nch(8,:)),'LineWidth',1.5)
%hold on,plot(t_wc(8,:),lc_sg)
xlabel('Time, [sec]')
ylabel('Force, [N]')
set(gca,'FontSize',20)
title('Aqua-FE 5 Second Wave: Channel Side Tension Time Series','FontSize',24)
xlim([0 40])
ylim([650 1400])

%% Load Cell Stuff

deploy_LC = load('zLC2_meta.mat');

% mean
figure,subplot(2,1,1),plot(f_cga,mn_loadch./2.20462,'o','Color',[0 .45 .74],'MarkerSize',12,'LineWidth',1.5)
set(gca,'FontSize',20)
title('Aqua-FE Waves and Current Test: Mean Mooring Line Loads','FontSize',24)
xlabel('Frequency, [Hz]')
set(gca,'YTick',[0:100:700])
ylabel('Force [lbs]')
%hold on,plot(1./NT_wvs,Nmn_loadch./2.20462,'s')
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca,'FontSize',20)
subplot(2,1,2),plot(deploy_LC.yd,deploy_LC.mn_lb,'LineWidth',1.5,'Color',[.85 .33 .1])
xlim([(deploy_LC.yd(1)-0.05) (deploy_LC.yd(end)+0.05)])
set(gca,'FontSize',20)
title('Deployment Mean Mooring Line Loads','FontSize',24)
datetickzoom('x','keepticks')
ylabel('Force [lbs]')
ylim([210 280])
xlabel('Date [mm/dd]')


%std dev
figure,subplot(2,1,1),plot(f_cga,std_loadch./2.20462,'o','MarkerSize',12,'LineWidth',1.5)
xlim([0.025 0.5])
set(gca,'YTick',[0:50:300])
set(gca,'FontSize',20)
title('Aqua-FE Waves and Current Test: Std Dev. Mooring Line Loads','FontSize',24)
ylabel('Force [lbs]')
xlabel('Frequency, [Hz]')
hold on,plot(1./NT_wvs,Nstd_loadch,'s')
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca,'FontSize',20)
subplot(2,1,2),plot(deploy_LC.yd,deploy_LC.std_lb,'LineWidth',1.5,'Color',[.85 .33 .1])
xlim([(deploy_LC.yd(1)-0.05) (deploy_LC.yd(end)+0.05)])
datetickzoom('x','keepticks')
ylabel('Force [lbs]')
xlabel('Date [mm/dd]')
set(gca,'FontSize',20)
title('Deployment Std Dev. Mooring Line Loads','FontSize',24)

%RAOs
dep_data = importdata('..\..\..\Field_Testing\Combined\Data\First\Full_ADCP_data.mat');
other_RAO = importdata('..\..\..\Field_Testing\Combined\Tools\RAO_dep.mat');

%Dates that the ADCP and LC's overlap IN EST
st_d = 14;
st_h= 18; %This is in EST
st_date = datenum(2016,01,st_d,st_h,00,00);

end_d = 29;
end_h = 09; %This is in EST
end_date = datenum(2016,01,end_d,end_h,00,00);

[~, indw] = min(abs(dep_data.datenums-st_date));
[~, indw2] = min(abs(dep_data.datenums-end_date));

[~, indL] = min(abs(deploy_LC.yd-st_date));
[~, indL2] = min(abs(deploy_LC.yd-end_date));

dates_LC = deploy_LC.yd(indL:indL2);
LC_mean = deploy_LC.mn_lb(1,indL:indL2);
LC_stddev = deploy_LC.std_lb(1,indL:indL2);
Sp_LC = deploy_LC.Sp_lc(indL:indL2,:);
freqs_LC = deploy_LC.bandf;

%Truncating data to line up
data1 = dep_data.data(indw:indw2,:);
dates_waves = dep_data.datenums(indw:indw2);
Wave_Spectra = dep_data.spectra(indw:indw2,:);

Waves_Spectra=1/10^6*Wave_Spectra(:,:).^2; 

Hs = data1(:,7);
H10 = data1(:,8);
H13 = data1(:,9);
Hmax = data1(:,10);
%Wave Period Information
Tp = data1(:,11);
T10 = data1(:,12);
T13 = data1(:,13);
Tmax = data1(:,14);
Depth = data1(:,16);
Cur_mag = data1(:,17);
Cur_dir = data1(:,18);

ttolw=datenum(0,0,0,0,0,10); %Within 30 sec
mw = 1;
for ii = 1:length(dates_LC)
    
    for jj = 1:length(dates_waves)
        
        if abs(dates_waves(jj)-dates_LC(ii))<ttolw
            H13_t(mw) = H13(jj);
            date_wv(mw) = dates_waves(jj);
            wvspec(mw) = .5*1027*9.8*trapz(Waves_Spectra(jj,:));
            mw = mw+1;
        end
    end    
    
end
% 
% for i = [1:(length(wvspec)/3)]
%        
%     if i ~= 29 && i<=28
%         m = i;
%     
%     wvspec(m,:) = wvspec(3*i-2); %This is 0 min output
%     
%     date_wv(m,:) = date_wv(3*i-2);
%     
%     Hs(m,:) = Hs(3*i-2);
%     
%     H13_t(m,:) = H13(3*i-2);
%         
%     Tp(m,:) = Tp(3*i-2); 
%     
%     T13(m,:) = T13(3*i-2); 
%     
%     elseif i ~= 29 && i>=30
%     m = i-1;
% 
%     wvspec(m,:) = wvspec(3*i-2); %This is 0 min output
%     
%     date_wv(m,:) = date_wv(3*i-2);
%     
%     Hs(m,:) = Hs(3*i-2);
%     
%     H13_t(m,:) = H13(3*i-2);
%     
%     Tp(m,:) = Tp(3*i-2); 
%     
%     T13(m,:) = T13(3*i-2); 
%     end
% 
% end

figure,plot(f_cga,RAO_shore./2.20462 ,'-','MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor',[0 .45 .74],'MarkerSize',12),hold on
plot(f_cga,RAO_chan./2.20462 ,'o','MarkerFaceColor',[.85 .33 .1],'MarkerEdgeColor',[.85 .33 .1],'MarkerSize',12,'LineWidth',1.5)
%plot(f_cga,RAO_shore_R,'--','Color',[.85 .33 .1],'LineWidth',1.5)
%plot(f_cga,RAO_chan_R,'x','MarkerFaceColor',[.85 .33 .1],'MarkerEdgeColor',[.85 .33 .1],'MarkerSize',10,'LineWidth',1.5)
hold on,plot(1./NT_wvs,NRAO_chan./2.20462,'s')
xlim([0.025 .5])
set(gca,'XTick',[0:.05:.5])
ylim([2000 14500])
set(gca,'YTick',[2000:1250:14500])
title('Aqua-FE Waves and Current Test: Mooring Line Load Dimensional RAOs','FontSize',22)
xlabel('Frequency, [Hz]')
ylabel('Response Amplitude, [lb/m]')
legend('Shore Side Load','Channel Side Load','Location','NorthWest')
set(gca,'OuterPosition',[-.002 -.018 1 .944])
set(gca,'FontSize',18)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca,'FontSize',18)

dim_load_dep = LC_mean./H13_t;

dim_load_dep(real(dim_load_dep)==Inf)=NaN;
dim_load_dep(real(dim_load_dep)==-Inf)=NaN;

% figure,plot(date_wv,dim_load_dep),hold on
% plot(date_wv(10),nanmean(dim_load_dep),'o',date_wv(20),mean(RAO_chan./2.20462),'o')
% datetickzoom('x','keepticks')

%% Plotting

% figure
% for id = 1:length(T_wvs)
% hold on
%     %Plot Channel Side Loading
% plot(t_wc(id,:),mid_load_lbch(id,:),'DisplayName',[num2str(T_wvs(id)),' Second Wave'])
% 
% title('Mid Water .2 m/s 4.75 Angled Cur & Wave Loading (Lbs)')
% 
% end
% legend(gca,'show')
% hold off

%% Plotting Each Heave Motion
% figure
% for id = 1:length(T_wvs)
%     hold on
% plot(t_wc(id,30:end-1),(cog(id,30:end)-nanmean(cog(id,30:end))),'DisplayName',[num2str(T_wvs(id)),' Second Wave'])
% title('Heave motion [m]')
% 
% end
% legend(gca,'show')
% hold off

% figure
% for id = 1:length(T_wvs)
%     hold on
% plot(t_wc(id,30:end-1),(IMU(id,30:end)-nanmean(IMU(id,30:end))),'DisplayName',[num2str(T_wvs(id)),' Second Wave'])
% 
% title('Vertical motion @ IMU location [m]')
% 
% end
% legend(gca,'show')
% hold off


%% AquaFE Spectral Information


% for ii = 1:length(k)
% 
%     [f_cg(ii,:),S_cg(ii,:)] = oneside_spec(cg1(ii,:),.2);
%     [f_IU(ii,:),S_IU(ii,:)] = oneside_spec(IU1(ii,:),.2);
%     [f_wave(ii,:),S_wave(ii,:)] = oneside_spec(wave_load1(ii,:),.2);
%     
% end
% 
%     S_c = nansum(S_cg,2);
% 
%     S_I = nansum(S_IU,2);
% 
%     S_w = nansum(S_wave,2);
% 
% 
% figure
% for ii = 1:length(k)
%     hold on
% semilogy(f_cg(ii,:),S_cg(ii,:),'DisplayName',[num2str(T_wvs(ii)),' Second Wave'])
% title('AquaFE Center of Gravity Spectra')
% xlabel('Wave Period')
% ylabel('Energy')
% 
% end
% hold off
% legend(gca,'show')
% 
% figure
% for ii = 1:length(k)
%     hold on
% semilogy(f_IU(ii,:),S_IU(ii,:),'DisplayName',[num2str(T_wvs(ii)),' Second Wave'])
% title('AquaFE IMU Location Spectra')
% xlabel('Wave Period')
% ylabel('Energy')
% 
% end
% hold off
% legend(gca,'show')
% 
% figure
% for ii = 1:length(k)
%     hold on
% semilogy(f_IU(ii,:),S_wave(ii,:),'DisplayName',[num2str(T_wvs(ii)),' Second Wave'])
% title('AquaFE Wave Spectra')
% xlabel('Wave Period')
% ylabel('Energy')
% 
% end
% hold off
% legend(gca,'show')
% 
%% AquaFE Spectral RAOs
% 
%     RAO_cogsp = sqrt(abs(S_c)./abs(S_w));
%     
%     RAO_IUsp = sqrt(abs(S_I)./abs(S_w));


%% AquaFE Amplitude RAOs

% figure,plot(f_cga,RAO_cog)
% xlim([0 .5])
% ylim([0 1.2])
% title('Aqua-FE Center of Gravity RAOs')
% xlabel('Frequency [Hz]')
% ylabel('Response')
% 
% figure,plot(f_IUa,RAO_IMU)
% xlim([0 .5])
% ylim([0 1.2])
% title('Aqua-FE IMU Location RAOs')
% xlabel('Frequency [Hz]')
% ylabel('Response')
% 
% figure,plot(f_cg40,RAO_cog40,'*',f_cg40,RAO_IMU40,'*')
% xlim([0 .5])
% ylim([0 1.2])
% title('Aqua-FE 40 cm Waves RAOs')
% xlabel('Frequency [Hz]')
% ylabel('Response')
% legend('Center of Gravity','IMU Location')

figure,
plot(T_wvs,L_sing)
title('Wavelength vs. Period')
xlabel('Wave Period')
ylabel('Wavelength')
grid

figure,plot(1./T_wvs,wave_steep1,'ok','MarkerSize',8,'LineWidth',1.5),hold on
plot(1./T_wvs,avgp,'Color',[0 .45 .74],'LineWidth',1.5)
plot(1./T_wvs,avgr,'Color',[.85 .33 .1],'LineWidth',1.5)
title('Aqua-FE Waves and Current Test: Wave Slope and Raft Angles','FontSize',22)
xlabel('Frequency, [Hz]')
ylabel('Wave Slope or Raft Angle, [rad]')
legend(' Wave Slope',' Pitch Angle',' Roll Angle','Location','NorthWest')
set(gca,'FontSize',18)
xlim([0.025 0.5])

%% Spectra and RAOs from Deployment


Hs_edges=0:max(other_RAO.H13)/6:max(other_RAO.H13);

heave_rao = squeeze(other_RAO.out_heave);

UCheave_rao = squeeze(other_RAO.out_heaveUC);

% Spectras
figure; 
plot(other_RAO.freqs(2:end-1),mean(other_RAO.Waves_Spectra(:,2:end-1,3),1));
hold on
plot(other_RAO.freqs(2:end-1),mean(other_RAO.Sj_z_uncorr(:,2:end-1),1));
plot(other_RAO.freqs(2:end-1),nanmean(other_RAO.Sj1side20(:,2:end-1,3),1));
title('Average of all Z-Displ Spectra and Wave Vertical Spectra')
xlabel('Frequency, Hz')
ylabel('Spectral Energy Density, m^2/Hz')
legend('Wave Spectra','Z Displacement Spectra','Corrected Z Displ Spectra')
hold off


figure,plot(other_RAO.freqs(2:end-1),nanmean(other_RAO.RAOs(:,:,3),1),'LineWidth',1.5)
hold on
plot(f_cga,RAO_cog,'o','MarkerSize',10,'MarkerFaceColor',[.85 .33 .1])
%plot(bfreqtest(2:end-1),nanmean(weird_RAO(:,2:end-1),1),'--')
ylabel('Response Amplitude')
xlabel('Frequency, [Hz]')
xlim([0 .5])
legend('Deployment Heave','Aqua-FE Heave')
set(gca, 'FontSize', 22)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca, 'FontSize', 22)
ylim([0 1.4])
set(gca, 'YTick', [0:.2:1.4])
title('Heave RAOs','FontSize', 28)
hold off


figure,plot(other_RAO.freqs(2:end-1),nanmean(other_RAO.RAOs_uncorr(:,:,3),1),'LineWidth',1.5)
hold on
plot(f_IUa,RAO_IMU,'x','MarkerSize',10,'MarkerFaceColor',[.85 .33 .1],'LineWidth',1.5)
ylabel('Response Amplitude')
xlabel('Frequency, [Hz]')
xlim([0.05 .5])
ylim([0 1.1])
legend('Deployment IMU Heave','Aqua-FE IMU Location Heave')
set(gca, 'FontSize', 22)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca, 'FontSize', 22)
title('IMU RAOs','FontSize', 28)
hold off

%% Pitch RAOs

dim_pitchdep_RAO = sqrt(other_RAO.Sj1side20(:,2:end-1,4)./other_RAO.Waves_Spectra(:,2:end-1,3));

dim_rolldep_RAO = sqrt(other_RAO.Sj1side20(:,2:end-1,5)./other_RAO.Waves_Spectra(:,2:end-1,3));

figure, 
plot(other_RAO.freqs(2:end-1),nanmean(other_RAO.Waves_Spectra(:,2:end-1,4),1),'ok','MarkerSize',10,'LineWidth',1.5);
hold on
plot(other_RAO.freqs(2:end-1),nanmean(other_RAO.Sj1side20(:,2:end-1,4),1),'LineWidth',1.5);
plot(other_RAO.freqs(2:end-1),nanmean(other_RAO.Sj1side20(:,2:end-1,5),1),'LineWidth',1.5);
xlabel('Frequency, Hz')
ylabel('Spectral Energy Density, [rad^2/Hz]')
set(gca,'FontSize',22)
% LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
% set(gca, 'FontSize', 22)
xlim([0.05 .5])
legend('Measured Deployment Wave Pitch','IMU Pitch Spectra','IMU Roll Spectra')
title('Date-Averaged Deployment Pitch and Roll Spectra and Wave Pitch Spectra','FontSize',28)
hold off


figure,plot(other_RAO.freqs(2:end-1),nanmean(other_RAO.RAOs(:,:,4),1),'Color',[0 .45 .74],'LineWidth',1.5)
hold on
plot(other_RAO.freqs(2:end-1),nanmean(other_RAO.RAOs(:,:,5),1),'--','Color',[.85 .33 .1],'LineWidth',1.5)
plot(f_cga,RAO_pitch,'o','MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor',[0 .45 .74],'MarkerSize',10)
plot(f_cga,RAO_roll,'x','LineWidth',1.5,'MarkerFaceColor',[.85 .33 .1],'MarkerEdgeColor',[.85 .33 .1],'MarkerSize',10)
ylabel('Response Amplitude')
xlabel('Frequency, [Hz]')
xlim([0.05 .5])
legend('Deployment Pitch','Deployment Roll','Aqua-FE Pitch','Aqua-FE Roll')
set(gca, 'FontSize', 22)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca, 'FontSize', 22)
ylim([0 3])
set(gca,'YTick',[0:0.5:4])
title('Pitch and Roll RAOs','FontSize',28)
hold off

plot_mask = ones(length(other_RAO.RAOs(:,1,3)),1);

for ii = 1:length(other_RAO.RAOs(:,1,3))
    
    for jj = 1:63
        
        if dim_pitchdep_RAO(ii,jj) == inf
            
            dim_pitchdep_RAO(ii,jj) = nan;
            plot_mask(ii) = 0;
        end
        
        if dim_rolldep_RAO(ii,jj) == inf
            
            dim_rolldep_RAO(ii,jj) = nan;
            
        end
        
    end
    
end
            
pitch_depdimplot = nanmean(dim_pitchdep_RAO,1);
roll_depdimplot = nanmean(dim_rolldep_RAO,1);
freq_plot = other_RAO.freqs(2:end-1);


%Dimensional Spectra Comparisons
figure; 
h_1 = semilogy(freq_plot,mean(other_RAO.Waves_Spectra(:,2:end-1,3),1));
hold on
h2 = semilogy(freq_plot,nanmean(other_RAO.Sj1side20(:,2:end-1,4),1));
hh = semilogy(freq_plot,nanmean(other_RAO.Sj1side20(:,2:end-1,5),1));
xlabel('Frequency, Hz')
ylabel('Spectral Energy Density, m^2/Hz')
legend('Wave Spectra','IMU Pitch Spectra','IMU Roll Spectra')
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca,'YTick',[0:0.5:4])
set(gca, 'FontSize', 20)
title('Average of all Z-Displ Spectra and Wave Vertical Spectra','FontSize',24)
hold off

%% Looking at WAVE spectra over each date

% figure,
% for ii = 1:90
%     
%     hold on
%     semilogy(other_RAO.freqs(2:end-1),(other_RAO.Waves_Spectra(ii,2:end-1,3)))
%     
% end
% hold off
% 
% title('Wave Spectra over each of the 90 days')
% xlabel('Frequency [Hz]')
% ylabel('Spectral Energy Density [m^2/Hz]')

figure,plot(freq_plot,pitch_depdimplot,'Color',[0 .45 .74],'LineWidth',1.5)
hold on
plot(freq_plot,roll_depdimplot,'--','Color',[.85 .33 .1],'LineWidth',1.5)
plot(1./T_wvs,RAO_pitch_dim','o','MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor',[0 .45 .74],'MarkerSize',8)
plot(1./T_wvs,RAO_roll_dim,'x','LineWidth',1.5,'MarkerFaceColor',[.85 .33 .1],'MarkerEdgeColor',[.85 .33 .1],'MarkerSize',8)
legend('Deployment Dimensional Pitch','Deployment Dimensional Roll','Aqua-FE Dimensional Pitch','Aqua-FE Dimensional Roll','Location','NorthWest')
xlabel('Frequency, [Hz]')
ylabel('Response Amplitude, [rad/m]')
set(gca, 'FontSize', 20)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca, 'FontSize', 20)
title('Dimensional Pitch and Roll RAOs','FontSize',24)
xlim([0 0.5])
ylim([0 0.6])
set(gca,'YTick',[0:.1:.6])
hold off

