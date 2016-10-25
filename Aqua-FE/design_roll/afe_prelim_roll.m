clc,clear,close all

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))
addpath(genpath('..\..\Field_Testing'));
addpath(genpath('..\..\Field_Testing\TJD_MAT_funs'));
addpath(genpath('..\..'));

pathdir = pwd;
dirData = dir(pathdir);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
fileList = {dirData(~dirIndex).name}';  %# Get a list of the files

    ids = 1:length(fileList);   

disp(' ')
disp('    THIS SCRIPT IS PREDICTED VALUES FROM AQUAFE IF THE WAVES HIT THE RAFT ON ITS LONG SIDE')    
pause(2)

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

T_wvs = [2 2.5 3 10/3 4 5 6 7 8 10 15 20];

t_trial = T_wvs.*25;

LOOP_L = t_trial./.005;
LOOP_S = t_trial./.004;

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
    %Choose any data you want to isolate
iso = 2;

RAO_cogroll = zeros(1,length(k));
RAO_IMUroll = zeros(1,length(k));

avg_loadshR = zeros(length(T_wvs));
avg_loadchR = zeros(length(T_wvs));

figure
for id=1:length(T_wvs)

   fil(id,:) = [fileList{id}];

mid_data(id,:,:) = importdata(fil(id,:));
%pd1(id,:) = fil(id,13:14);
%pd2(id,:) = fil(id,16);

t_wc(id,:) = mid_data(id,2:end-1,1);
mid_stress_shore(id,:) = mid_data(id,2:end-1,2);
mid_stress_chan(id,:) = mid_data(id,2:end-1,3);

    %Convert Pa to N and multiply by C.S. area (in inches) to get lbs  

mid_load_Nsh(id,:) = mid_stress_shore(id,:).*A_LC;
mid_load_Nch(id,:) = mid_stress_chan(id,:).*A_LC;
%mid_load_lbsh(id,:) = mid_load_Nsh(id,:)*.22481;
%mid_load_lbch(id,:) = mid_load_Nch(id,:)*.22481;


cog(id,:) = mean(mid_data(id,2:end-1,4:7),3);
IMU(id,:) = mid_data(id,2:end-1,8);

rol_dif(id,:) = mid_data(id,2:end-1,5)-mid_data(id,2:end-1,7);
rol_dif2(id,:) = mid_data(id,2:end-1,4)-mid_data(id,2:end-1,6);

roll1(id,:) = asin(rol_dif(id,:)./1.2192);
roll2(id,:) = asin(rol_dif2(id,:)./1.2192);

roll(id,:) = (roll1(id,30:end) + roll2(id,30:end))/2;

pit_dif1 = mid_data(id,2:end-1,6)-mid_data(id,2:end-1,7);
pit_dif2 = mid_data(id,2:end-1,4)-mid_data(id,2:end-1,5);

pitch1= asin(pit_dif1./.2032);
pitch2 = asin(pit_dif2./.2032);

pitch(id,:) = (pitch1(30:end) + pitch2(30:end))/2;

%% Shortern series for spectra analysis
cg1(id,:) = cog(id,30:end)-mean(cog(id,30:end));
IU1(id,:) = IMU(id,30:end)-mean(IMU(id,30:end));

load_shN_R(id,:) = mid_load_Nsh(id,30:end);% - mean(mid_load_Nsh(id,30:end));
load_chN_R(id,:) = mid_load_Nch(id,30:end);% - mean(mid_load_Nch(id,30:end));

sig(id) = 2*pi./T_wvs(id);

wave_steep1(id,:) = a.*k(id);

[avgcog(id)] = amps_CS(t_wc(id,30:end),cg1(id,:),1,length(cg1(1,:)),T_wvs(id));
[avgIMU(id)] = amps_CS(t_wc(id,30:end),IU1(id,:),1,length(IU1(1,:)),T_wvs(id));    
    
RAO_cogroll(id) = avgcog(id)/a;
RAO_IMUroll(id) = avgIMU(id)/a;

[avgp(id)] = amps_CS(t_wc(id,30:end),pitch(id,:),1,length(pitch(1,:)),T_wvs(id));
[avgr(id)] = amps_CS(t_wc(id,30:end),roll(id,:),1,length(roll(1,:)),T_wvs(id));

RAO_pitchroll(id) = avgp(id)./wave_steep1(id);
RAO_rollroll(id) = avgr(id)./wave_steep1(id);

RAO_pitch_dimroll(id) = avgp(id)./a;
RAO_roll_dimroll(id) = avgr(id)./a;

if id <= 2
    
[avg_loadshR(id)] = amps_CS(t_wc(id,79:end),load_shN_R(id,50:end),1,length(load_shN_R(1,50:end)),T_wvs(id));
[avg_loadchR(id)] = amps_CS(t_wc(id,79:end),load_chN_R(id,50:end),1,length(load_chN_R(1,50:end)),T_wvs(id));

avg_loadshR(id) = avg_loadshR(id) + mean(mid_load_Nsh(id,50:end));
avg_loadchR(id) = avg_loadchR(id) + mean(mid_load_Nch(id,50:end));

else
    
[avg_loadshR(id)] = amps_CS(t_wc(id,129:end),load_shN_R(id,100:end),1,length(load_shN_R(1,100:end)),T_wvs(id));
[avg_loadchR(id)] = amps_CS(t_wc(id,129:end),load_chN_R(id,100:end),1,length(load_chN_R(1,100:end)),T_wvs(id));

avg_loadshR(id) = avg_loadshR(id) + mean(mid_load_Nsh(id,100:end));
avg_loadchR(id) = avg_loadchR(id) + mean(mid_load_Nch(id,100:end));  
    
end


RAO_shore_R(id) = avg_loadshR(id)./a;
RAO_chan_R(id) = avg_loadchR(id)./a;

clear('mid_data','t_wc')
clear('mid_stress_shore','mid_stress_chan','mid_load_Nsh','mid_load_Nch')
clear ('rol_dif','rol_dif2','roll1','roll2','roll')
clear('pit_dif1','pit_dif2','pitch1','pitch2','pitch')
clear ('cog','IMU','cg1','IU1','load_shN_R','load_chN_R')

end

%% AquaFE Amplitude RAOs

%figure,plot(load_chN_R(2,:)),hold on,plot(avg_loadchR(2),'o'),hold off

%RAO_cog(1) = .2750;

f_cg1 = 1./T_wvs;

f_cga = f_cg1;
f_IUa = f_cga;

figure,plot(f_cga,avgcog,'o',f_cga,a*ones(1,length(avgcog)),':')

figure,plot(f_cga,RAO_shore_R,'x','MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor',[0 .45 .74],'MarkerSize',10),hold on,hold on
plot(f_cga,RAO_chan_R,'x','MarkerFaceColor',[.85 .33 .1],'MarkerEdgeColor',[.85 .33 .1],'MarkerSize',10,'LineWidth',1.5)
xlim([0 .5])
set(gca,'XTick',[0:.05:.5])
% ylim([0 1.2])
title('Aqua-FE Mooring Line Load Dimensional RAOs','FontSize',18)
xlabel('Frequency [Hz]')
ylabel('Response Amplitude [N/m]')
legend('Shore Side Load','Channel Side Load')
set(gca,'OuterPosition',[-.002 -.018 1 .944])
set(gca,'FontSize',16)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca,'FontSize',16)


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

figure,h = plot(1./T_wvs,wave_steep1,1./T_wvs,avgp);
hold on
plot(1./T_wvs,avgr)
title('Long Side Wave Incidence AquaFE Wave Steepness and Predicted Raft Angles')
xlabel('Frequency [Hz]')
ylabel('Wave Steepness or Raft Angle')
legend('Wave Steepness','Pitch Angle','Roll Angle')

%set(h(1),'LineColor','blue')

%% Spectra and RAOs from Deployment

dep_data = importdata('RAO_dep.mat');

Hs_edges=0:max(dep_data.H13)/6:max(dep_data.H13);

heave_rao = squeeze(dep_data.out_heave);

UCheave_rao = squeeze(dep_data.out_heaveUC);

% Spectras
figure; 
plot(dep_data.freqs(2:end-1),mean(dep_data.Waves_Spectra(:,2:end-1,3),1));
hold on
plot(dep_data.freqs(2:end-1),mean(dep_data.Sj_z_uncorr(:,2:end-1),1));
plot(dep_data.freqs(2:end-1),nanmean(dep_data.Sj1side20(:,2:end-1,3),1));
title('Average of all Deployment Z-Displ Spectra and Wave Vertical Spectra')
xlabel('Frequency, Hz')
ylabel('Spectral Energy Density, m^2/Hz')
legend('Wave Spectra','Z Displacement Spectra','Corrected Z Displ Spectra')
hold off


figure,plot(dep_data.freqs(2:end-1),nanmean(dep_data.RAOs(:,:,3),1))
hold on
plot(dep_data.freqs(2:end-1),nanmean(dep_data.RAOs_uncorr(:,:,3),1))
plot(f_cga,RAO_cogroll,'*')
plot(f_IUa,RAO_IMUroll,'*')
title('Mean Heave RAOs from deployment and Amplitude RAOs from Predicted LS Wave Incidence')
ylabel('RAOs')
xlabel('Frequency [Hz]')
xlim([0 .5])
legend('Observed COG Corrected','Observed Uncorrected','AquaFE Center of Gravity','AquaFE IMU location')
set(gca, 'ColorOrderIndex', 1)
hold off


%% Pitch RAOs

dim_pitchdep_RAO = sqrt(dep_data.Sj1side20(:,2:end-1,4)./dep_data.Waves_Spectra(:,2:end-1,3));

dim_rolldep_RAO = sqrt(dep_data.Sj1side20(:,2:end-1,5)./dep_data.Waves_Spectra(:,2:end-1,3));

figure, 
plot(dep_data.freqs(2:end-1),mean(dep_data.Waves_Spectra(:,2:end-1,5),1));
hold on
plot(dep_data.freqs(2:end-1),mean(dep_data.Sj1side20(:,2:end-1,4),1));
plot(dep_data.freqs(2:end-1),mean(dep_data.Sj1side20(:,2:end-1,5),1));
title('Average of all Deployment IMU Pitch and Roll Spectra and Wave Pitch Spectra')
xlabel('Frequency, Hz')
ylabel('Spectral Energy Density, rad^2/Hz')
legend('Measured Deployment Wave Pitch','IMU Pitch Spectra','IMU Roll Spectra')
hold off

figure,plot(dep_data.freqs(2:end-1),nanmean(dep_data.RAOs(:,:,4),1))
hold on
plot(dep_data.freqs(2:end-1),nanmean(dep_data.RAOs(:,:,5),1))
plot(f_cga,RAO_pitchroll,'*b')
plot(f_cga,RAO_rollroll,'*r')
title('Mean Pitch RAOs from deployment and Predicted LS Wave Incidence AquaFE Amplitude RAOs')
ylabel('RAOs')
xlabel('Frequency [Hz]')
xlim([0 .5])
legend('Observed Pitch','Observed Roll','AquaFE Pitch','AquaFE Roll')
set(gca, 'ColorOrderIndex', 1)
hold off

plot_mask = ones(length(dep_data.RAOs(:,1,3)),1);

for ii = 1:length(dep_data.RAOs(:,1,3))
    
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
freq_plot = dep_data.freqs(2:end-1);


%Dimensional Spectra Comparisons
figure; 
h_1 = semilogy(freq_plot,mean(dep_data.Waves_Spectra(:,2:end-1,3),1));
hold on
h2 = semilogy(freq_plot,nanmean(dep_data.Sj1side20(:,2:end-1,4),1));
hh = semilogy(freq_plot,nanmean(dep_data.Sj1side20(:,2:end-1,5),1));
title('Average of all Deployment Z-Displ Spectra and Wave Vertical Spectra')
xlabel('Frequency, Hz')
ylabel('Spectral Energy Density, m^2/Hz')
legend('Wave Spectra','IMU Pitch Spectra','IMU Roll Spectra')
hold off

%% Looking at WAVE spectra over each date

% figure,
% for ii = 1:90
%     
%     hold on
%     semilogy(dep_data.freqs(2:end-1),(dep_data.Waves_Spectra(ii,2:end-1,3)))
%     
% end
% hold off
% 
% title('Wave Spectra over each of the 90 days')
% xlabel('Frequency [Hz]')
% ylabel('Spectral Energy Density [m^2/Hz]')

figure,plot(freq_plot,pitch_depdimplot,'--b')
hold on
plot(freq_plot,roll_depdimplot,'--r')
plot(1./T_wvs,RAO_pitch_dimroll','ob')
plot(1./T_wvs,RAO_roll_dimroll,'or')
legend('Deployment Dimensional Pitch RAO','Deployment Dimensional Roll RAO','Predicted AquaFE Dimensional Pitch RAO','Predicted AquaFE Dimensional Roll RAO','Location','NorthWest')
title('Dimensional Pitch and Roll RAOs','FontSize',20)
xlabel('Frequency [Hz]','FontSize',17.6)
ylabel('Response Amplitude [rad/m]','FontSize',17.6)
hold off

savename=['..\pred_S\roll_prediction'];

save(savename,'avgr','RAO_cogroll','RAO_IMUroll','RAO_rollroll','RAO_roll_dimroll','RAO_shore_R','RAO_chan_R')

clc
disp(' ')
disp('    THIS SCRIPT IS PREDICTED VALUES FROM AQUAFE IF THE WAVES HIT THE RAFT ON ITS LONG SIDE')
