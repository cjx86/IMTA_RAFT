clc,clear,close all

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))
addpath(genpath('..\Field_Testing'));
addpath(genpath('..\Field_Testing\TJD_MAT_funs'));
addpath(genpath('..\..'));

pathdir = pwd;
dirData = dir(pathdir);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
fileList = {dirData(~dirIndex).name}';  %# Get a list of the files

    ids = 1:length(fileList);   

disp(' ')
disp('  THIS SCRIPT IS PREDICTED VALUES FROM AQUAFE IF THE WAVES HIT THE RAFT ON ITS SHORT SIDE')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mid Water Waves 1/3 Max Cur

RAO_cogpitch = zeros(1,length(k));
RAO_IMUpitch = zeros(1,length(k));

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

load_shN_P(id,:) = mid_load_Nsh(id,30:end);
load_chN_P(id,:) = mid_load_Nch(id,30:end);

sig(id) = 2*pi./T_wvs(id);

wave_steep1(id,:) = a.*k(id);

[avgcog(id)] = amps_CS(t_wc(id,30:end),cg1(id,:),1,length(cg1(1,:)),T_wvs(id));
[avgIMU(id)] = amps_CS(t_wc(id,30:end),IU1(id,:),1,length(IU1(1,:)),T_wvs(id));    
    
RAO_cogpitch(id) = avgcog(id)/a;
RAO_IMUpitch(id) = avgIMU(id)/a;

[avgp(id)] = amps_CS(t_wc(id,30:end),pitch(id,:),1,length(pitch(1,:)),T_wvs(id));
[avgr(id)] = amps_CS(t_wc(id,30:end),roll(id,:),1,length(roll(1,:)),T_wvs(id));

RAO_pitchpitch(id) = avgp(id)./wave_steep1(id);
RAO_rollpitch(id) = avgr(id)./wave_steep1(id);

RAO_pitch_dim(id) = avgp(id)./a;
RAO_roll_dimpitch(id) = avgr(id)./a;

if id <= 2
    
[avg_loadshP(id)] = amps_CS(t_wc(id,79:end),load_shN_P(id,50:end),1,length(load_shN_P(1,50:end)),T_wvs(id));
[avg_loadchP(id)] = amps_CS(t_wc(id,79:end),load_chN_P(id,50:end),1,length(load_chN_P(1,50:end)),T_wvs(id));

avg_loadshP(id) = avg_loadshP(id) + mean(mid_load_Nsh(id,50:end));
avg_loadchP(id) = avg_loadchP(id) + mean(mid_load_Nch(id,50:end));

else
    
[avg_loadshP(id)] = amps_CS(t_wc(id,129:end),load_shN_P(id,100:end),1,length(load_shN_P(1,100:end)),T_wvs(id));
[avg_loadchP(id)] = amps_CS(t_wc(id,129:end),load_chN_P(id,100:end),1,length(load_chN_P(1,100:end)),T_wvs(id));

avg_loadshP(id) = avg_loadshP(id) + mean(mid_load_Nsh(id,100:end));
avg_loadchP(id) = avg_loadchP(id) + mean(mid_load_Nch(id,100:end));  
    
end

RAO_shore_P(id) = avg_loadshP(id)./a;
RAO_chan_P(id) = avg_loadchP(id)./a;

clear('mid_data','t_wc')
clear('mid_stress_shore','mid_stress_chan','mid_load_Nsh','mid_load_Nch')
clear ('rol_dif','rol_dif2','roll1','roll2','roll')
clear('pit_dif1','pit_dif2','pitch1','pitch2','pitch')
clear ('cog','IMU','cg1','IU1','load_shN_P','load_chN_P')

end

f_cg1 = 1./T_wvs;

f_cga = f_cg1;
f_IUa = f_cga;

figure,plot(f_cga,avgcog,'o',f_cga,a*ones(1,length(avgcog)),':')

load('roll_prediction.mat')

figure,plot(f_cga,RAO_shore_P./2.20462,'-','MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor',[0 .45 .74],'MarkerSize',12,'LineWidth',2),hold on
plot(f_cga,RAO_chan_P./2.20462,'o','MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor',[0 .45 .74],'MarkerSize',12,'LineWidth',2)
plot(f_cga,RAO_shore_R./2.20462,'--','Color',[.85 .33 .1],'LineWidth',2)
plot(f_cga,RAO_chan_R./2.20462,'x','MarkerFaceColor',[.85 .33 .1],'MarkerEdgeColor',[.85 .33 .1],'MarkerSize',12,'LineWidth',2)
xlim([0.05 .5])
%set(gca,'XTick',[0:.05:.5])
ylim([1500 3250])
set(gca,'YTick',[1500:250:3250])
xlabel('Frequency, [Hz]')
ylabel('Response Amplitude, [lb/m]')
legend('Shore Side Load Pitch','Channel Side Load Pitch','Shore Side Load Roll','Channel Side Load Roll','Location','NorthWest')
set(gca,'OuterPosition',[-.002 -.018 1 .944])
set(gca,'FontSize',22)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca,'FontSize',22)
title('Aqua-FE Mooring Line Load Dimensional RAOs','FontSize',28)

%Cneter of Gravity Heave
figure,plot(f_cga,RAO_cogpitch,'o','MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor',[0 .45 .74],'MarkerSize',12),hold on
plot(f_cga,RAO_cogroll,'x','MarkerFaceColor',[.85 .33 .1],'MarkerEdgeColor',[.85 .33 .1],'MarkerSize',12,'LineWidth',1.5)
xlabel('Frequency [Hz]')
ylabel('Response Amplitude')
xlim([0.05 .5])
set(gca,'XTick',[0:.05:.5])
ylim([0 1.2])
set(gca,'YTick',[0:.2:1.2])
legend('Short Side Wave Incidence Heave','Long Side Wave Incidence Heave')
%set(gca,'OuterPosition',[-.002 -.018 1 .944])
set(gca,'FontSize',22)
title('Aqua-FE Heave RAOs','FontSize',28)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca,'FontSize',22)

figure,
plot(f_cga,RAO_pitchpitch,'o','MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor',[0 .45 .74],'MarkerSize',12),hold on;
plot(f_cga,RAO_rollroll,'x','MarkerFaceColor',[.85 .33 .1],'MarkerEdgeColor',[.85 .33 .1],'MarkerSize',12,'LineWidth',1.5)
ylabel('RAOs')
xlabel('Frequency [Hz]')
xlim([0.05 .5])
ylim([0 3.75])
set(gca,'XTick',[0:.05:.5])
set(gca,'YTick',[0:.5:4])
legend('Aqua-FE Pitch','Aqua-FE Roll','Location','NorthWest')
set(gca, 'ColorOrderIndex', 1)
set(gca,'FontSize',22)
title('Aqua-FE Predicted Pitch and Roll RAOs','FontSize',28)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca,'FontSize',22)
hold off

figure,
plot(T_wvs,L_sing)
title('Wavelength vs. Period')
xlabel('Wave Period')
ylabel('Wavelength')
grid

figure,plot(1./T_wvs,wave_steep1,'o','LineWidth',1.5,'MarkerSize',12),hold on
plot(1./T_wvs,avgp,'--','LineWidth',1.5)
plot(1./T_wvs,avgr,'LineWidth',1.5)
xlabel('Frequency, [Hz]')
ylabel('Wave Slope or Raft Angle, [rad]')
legend('Wave Slope','Pitch Angle','Roll Angle')
xlim([.05 0.5])
ylim([0 0.12])
set(gca,'XTick',[0:.05:.5])
set(gca,'YTick',[0:.02:.12])
set(gca,'FontSize',22)
set(gca,'OuterPosition',[-.002 -.018 1 .944])
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]')
set(gca,'FontSize',22)
%set(h(1),'LineColor','blue')
title('Aqua-FE Wave Slope and Aqua-FE Raft Angles','FontSize',28)
hold off

%% Pitch RAOs

fig4 = figure;
plot(f_cga,RAO_pitch_dim,'o','MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor',[0 .45 .74],'MarkerSize',12),hold on;
plot(f_cga,RAO_roll_dimroll,'x','MarkerFaceColor',[.85 .33 .1],'MarkerSize',12,'LineWidth',1.5)
ylabel('Response Amplitude [rad/m]')
xlabel('Frequency [Hz]')
xlim([0.05 .5])
ylim([0 .45])
set(gca,'YTick',[0:.05:.45])
legend('AquaFE Dimensional Pitch','AquaFE Dimensional Roll','Location','NorthWest')
set(gca, 'ColorOrderIndex', 1)
set(gca,'OuterPosition',[-.002 -.018 1 .944])
set(gca,'FontSize',22)
title('Aqua-FE Predicted Dimensional Pitch and Roll RAOs','FontSize',28)
LinkTopAxisData(1./[20 15 10 8 6 4 3 2 1],[20 15 10 8 6 4 3 2 1],'Period [sec]');
set(gca,'FontSize',22)
hold off
