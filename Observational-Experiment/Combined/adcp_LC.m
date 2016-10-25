% Comparison of full ADCP and Loadcell Data Sets

clear;clc;close all;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))
addpath(genpath('../'));
addpath(genpath('../functions'));

load('..\Full_ADCP_data.mat');
LC_meta=load('..\LC2_meta.mat');

LC_meta.yd = LC_meta.yd-datenum(0,0,0,5,0,0); %Corrects LC dates to EST

%Dates that the ADCP and LC's overlap IN EST
st_d = 14;
st_h= 18; %This is in EST
st_date = datenum(2016,01,st_d,st_h,00,00);

end_d = 29;
end_h = 09; %This is in EST
end_date = datenum(2016,01,end_d,end_h,00,00);

    %Find closest this date in the datenums array (min,sec could be off)
[~, indc] = min(abs(date_time-st_date));
[~, indc2] = min(abs(date_time-end_date));

[~, indw] = min(abs(datenums-st_date));
[~, indw2] = min(abs(datenums-end_date));

[~, indL] = min(abs(LC_meta.yd-st_date));
[~, indL2] = min(abs(LC_meta.yd-end_date));

%Truncating data to line up
data1 = data(indw:indw2,:);
dates_waves = datenums(indw:indw2);
Wave_Spectra = spectra(indw:indw2,:);

decl=-15.31667;

velmag = velmag_sg(indc:indc2);
veldir = veldir_sg_deg(indc:indc2);
dates_cur = date_time(indc:indc2);
w_lev = depth(indc:indc2);
MSL = mean(w_lev);

dates_LC = LC_meta.yd(indL:indL2);
LC_mean = LC_meta.mn_lb(1,indL:indL2);
LC_stddev = LC_meta.std_lb(1,indL:indL2);
Sp_LC = LC_meta.Sp_lc(indL:indL2,:);
freqs_LC = LC_meta.bandf;

LC_mn_conf_up = LC_mean + LC_stddev*2.92/sqrt(6000);
LC_mn_conf_low = LC_mean - LC_stddev*2.92/sqrt(6000);

figure,plot(dates_LC,LC_mean,dates_LC,LC_mn_conf_up,'--',dates_LC,LC_mn_conf_low,'--')

%Wave Height information
Hs = data1(:,7);
H10 = data1(:,8);
H13 = data1(:,9);
Hmax = data1(:,10);
%Wave Period Information
Tp = data1(:,11);
T10 = data1(:,12);
T13 = data1(:,13);
Tmax = data1(:,14);

Dp = data1(:,15)+decl;
Depth = data1(:,16);
Cur_mag = data1(:,17);
Cur_dir = data1(:,18);

Waves_Spectra=1/10^6*Wave_Spectra(:,:).^2; 

%% Plotting of all Current Data

figure,plot(dates_cur,w_lev,'LineWidth',1.75)
%plot(dates_waves,Depth,'LineWidth',1.75)
title('Tidal Elevation','FontSize',24)
datetickzoom('x','keepticks')
xlim([dates_cur(1)-.1 dates_cur(end)+.1])
ylabel('Water Depth [m]')
ylim([7 11])
set(gca,'FontSize',20)
xlabel('Date [mm/dd]')
linkaxes(findall(gcf,'type','axes'), 'x');
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)

%Plot of water level and current direction
figure,subplot(3,1,1),plot(dates_cur,velmag,'LineWidth',1.5,'Color',[0 .45 .74])
hold on,plot(dates_cur(9420).*ones(length(0:.01:.8)),[0:.01:.8],'k--')
plot(dates_cur(10198).*ones(length(0:.01:.8)),[0:.01:.8],'k--')
plot(dates_cur(10900).*ones(length(0:.01:.8)),[0:.01:.8],'k--')
plot(dates_cur(11692).*ones(length(0:.01:.8)),[0:.01:.8],'k--')
ylim([0 0.65])
title('Current Magnitude','FontSize',24)
ylabel('Velocity Magnitude [m/s]')
datetickzoom('x','keepticks')
set(gca,'FontSize',20)
subplot(3,1,2),plot(dates_cur,veldir,'LineWidth',1.5,'Color',[.85 .33 .1]),datetickzoom('x','keepticks')
hold on,plot(dates_cur(9420).*ones(length(-175:1:75)),[-175:1:75],'k--')
plot(dates_cur(10198).*ones(length(-175:1:75)),[-175:1:75],'k--')
plot(dates_cur(10900).*ones(length(-175:1:75)),[-175:1:75],'k--')
plot(dates_cur(11692).*ones(length(-175:1:75)),[-175:1:75],'k--')
ylabel('Velocity Direction [deg]')
title('Current Direction','FontSize',24)
ylim([-175 75])
set(gca,'FontSize',20)
subplot(3,1,3),plot(dates_cur,w_lev,'LineWidth',1.5,'Color',[0.9290 0.6940 0.1250])
hold on,plot(dates_cur(9420).*ones(length(6:1:12)),[6:1:12],'k--')
plot(dates_cur(10198).*ones(length(6:1:12)),[6:1:12],'k--')
plot(dates_cur(10900).*ones(length(6:1:12)),[6:1:12],'k--')
plot(dates_cur(11692).*ones(length(6:1:12)),[6:1:12],'k--')
title('Tidal Elevation','FontSize',24)
datetickzoom('x','keepticks')
ylabel('Water Depth [m]')
ylim([7 11])
set(gca,'YTick',[7 8 9 10 11])
set(gca,'FontSize',20)
xlabel('Date [mm/dd]')
linkaxes(findall(gcf,'type','axes'), 'x');
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
%suptitle('Water Depth, Current Direction and Magnitude','FontSize',22)

%EBB FLOW HIGHLIGHTS
figure,subplot(3,1,1),plot(dates_cur,velmag,'Color',[0 .45 .74],'Linewidth',1.5)
hold on,plot(dates_cur(9661).*ones(length(0:.01:.8)),[0:.01:.8],'k','Linewidth',1.5)
plot(dates_cur(9801).*ones(length(0:.01:.8)),[0:.01:.8],'k--')
plot(dates_cur(10431).*ones(length(0:.01:.8)),[0:.01:.8],'k','Linewidth',1.5)
plot(dates_cur(10622).*ones(length(0:.01:.8)),[0:.01:.8],'k--')
plot(dates_cur(11155).*ones(length(0:.01:.8)),[0:.01:.8],'k','Linewidth',1.5)
plot(dates_cur(11312).*ones(length(0:.01:.8)),[0:.01:.8],'k--')
plot(dates_cur(11920).*ones(length(0:.01:.8)),[0:.01:.8],'k','Linewidth',1.5)
plot(dates_cur(12091).*ones(length(0:.01:.8)),[0:.01:.8],'k--')
plot(dates_cur(8931).*ones(length(0:.01:.8)),[0:.01:.8],'k','Linewidth',1.5)
plot(dates_cur(9067).*ones(length(0:.01:.8)),[0:.01:.8],'k--')
ylim([0 0.65])
title('Current Magnitude','FontSize',24)
ylabel('Velocity Magnitude [m/s]')
datetickzoom('x','keepticks')
set(gca,'FontSize',20)
subplot(3,1,2),plot(dates_cur,veldir,'Color',[.85 .33 .1],'Linewidth',1.5),datetickzoom('x','keepticks')
hold on,plot(dates_cur(9661).*ones(length(-175:1:75)),[-175:1:75],'k','Linewidth',1.5)
plot(dates_cur(9801).*ones(length(-175:1:75)),[-175:1:75],'k--')
plot(dates_cur(10431).*ones(length(-175:1:75)),[-175:1:75],'k','Linewidth',1.5)
plot(dates_cur(10622).*ones(length(-175:1:75)),[-175:1:75],'k--')
plot(dates_cur(11155).*ones(length(-175:1:75)),[-175:1:75],'k','Linewidth',1.5)
plot(dates_cur(11312).*ones(length(-175:1:75)),[-175:1:75],'k--')
plot(dates_cur(8931).*ones(length(-175:1:75)),[-175:1:75],'k','Linewidth',1.5)
plot(dates_cur(9067).*ones(length(-175:1:75)),[-175:1:75],'k--')
plot(dates_cur(11920).*ones(length(-175:1:75)),[-175:1:75],'k','Linewidth',1.5)
plot(dates_cur(12091).*ones(length(-175:1:75)),[-175:1:75],'k--')
ylabel('Velocity Direction [deg]')
title('Current Direction','FontSize',24)
ylim([-175 75])
set(gca,'FontSize',20)
subplot(3,1,3),plot(dates_cur,w_lev,'Color',[0.9290 0.6940 0.1250],'Linewidth',1.5)
hold on,plot(dates_cur(9661).*ones(length(6:1:12)),[6:1:12],'k','Linewidth',1.5)
plot(dates_cur(9801).*ones(length(6:1:12)),[6:1:12],'k--')
plot(dates_cur(10431).*ones(length(6:1:12)),[6:1:12],'k','Linewidth',1.5)
plot(dates_cur(10622).*ones(length(6:1:12)),[6:1:12],'k--')
plot(dates_cur(11155).*ones(length(6:1:12)),[6:1:12],'k','Linewidth',1.5)
plot(dates_cur(11305).*ones(length(6:1:12)),[6:1:12],'k--')
plot(dates_cur(8931).*ones(length(6:1:12)),[6:1:12],'k','Linewidth',1.5)
plot(dates_cur(9067).*ones(length(6:1:12)),[6:1:12],'k--')
plot(dates_cur(11920).*ones(length(6:1:12)),[6:1:12],'k','Linewidth',1.5)
plot(dates_cur(12091).*ones(length(6:1:12)),[6:1:12],'k--')
title('Tidal Elevation','FontSize',24)
datetickzoom('x','keepticks')
ylabel('Water Depth [m]')
ylim([7 11])
set(gca,'YTick',[7 8 9 10 11])
set(gca,'FontSize',20)
xlabel('Date [mm/dd]')
linkaxes(findall(gcf,'type','axes'), 'x');
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)


%% Check Alignment of data: Check LC, ADCP alignment and plot
 
%      %Filtering LC data
[LC_again,LC_filt,LC_filt_mean,~,numgg,~]=stddev_filterCS(dates_LC,LC_mean,6,3);

% 
figure; %Plot current effects  %FIGURE (4)
subplot(311),[a1x,~,~] = plotyy(dates_cur,velmag,dates_cur,w_lev);
datetickzoom%('keepticks')
ylabel('Velocity Mag (m/s)')
set(a1x(2),'YLim',[7 11])
set(a1x(1),'YLim',[0 0.7])
legend('Current Magnitude','Water Level','Location','NorthWest')
subplot(312),plot(dates_cur,veldir,'DisplayName','Current Direction'),datetickzoom('x','keepticks'),legend('Show')
ylabel('Velocity Dir (degrees)')
subplot(313),plot(dates_LC,LC_filt,'DisplayName','LC2 Mean')
ylim([190 285])
hold on
datetickzoom('x','keepticks')
legend('Show')
xlabel('Date')
ylabel('Mean Tension, lb')
linkaxes(findall(gcf,'type','axes'), 'x');
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
suptitle('Current Effects on Tension')
hold off

figure,[a11,~,~]=plotyy(dates_cur,w_lev,dates_LC,LC_filt);
set(a11(1),'YLim',[7 11])
set(a11(2),'YLim',[210 280])
set(a11(2),'YTick',210:5:280)
datetickzoom('x','keepticks')
xlabel('Date')
ylabel(a11(1),'Water Level, m')
ylabel(a11(2),'Mean Tension, lb')
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)


%% Check How Different Factors Affect Tensions

figure; %Plot wave effects FIGURE (5)
subplot(311),plot(dates_waves,Hs,dates_waves,H13)
legend('Hs','H_{1/3}')
ylim([0 1.25])
datetickzoom('x','keepticks')
ylabel('Sig Wave Ht (m)')
subplot(312),plot(dates_waves,Dp,'DisplayName','Dp')
legend('Show')
datetickzoom('x','keepticks')
ylabel('Dom Wave Dir (degrees)')
subplot(313),plot(dates_LC,LC_stddev,'DisplayName','LC2 Std Dev');
datetickzoom('x','keepticks')
ylim([0 25])
hold on
legend('Show')
ylabel('Std Dev of Tension (lbs)')
xlabel('Date')
linkaxes(findall(gcf,'type','axes'), 'x');
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
suptitle('Wave Effects on Tension')
hold off

%Plots of waves and current overlayed
figure,subplot(311)
[ax1,~,~] = plotyy(dates_waves,H13,dates_cur,velmag);
datetickzoom('x','keepticks')
legend('H_{1/3}','Velocity Magnitude')
ylabel(ax1(1),'Sig Wave Ht')
ylabel(ax1(2),'Current Magnitude')
subplot(312)
xlim([dates_cur(1) dates_cur(end)])
[ax2,~,~] = plotyy(dates_waves,Dp,dates_cur,veldir);
datetickzoom('x','keepticks')
legend('Dp','Cur Dir')
ylabel(ax2(1),'Wave Dir')
ylabel(ax2(2),'Cur Dir')
subplot(3,1,3),[hm,~,~] = plotyy(dates_LC,LC_filt,dates_LC,LC_stddev);
legend('LC2 Mean','LC2 Std Dev')
datetickzoom('x','keepticks')
set(hm(1),'YLim',[190 300])
set(hm(2),'YLim',[0 25])
ylabel(hm(1),'Mean Tension (lbs)')
ylabel(hm(2),'Std Dev of Tension (lbs)')
suptitle('Date Overlay of all Events')
linkaxes(findall(gcf,'type','axes'), 'x');
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)

%Plot of current directly overlayed with mean tension
figure,[axu,~,~] = plotyy(dates_cur,velmag,dates_LC,LC_filt);
datetickzoom('x','keepticks')
set(axu(1),'XLim',[dates_LC(50) dates_LC(end-100)])
set(axu(1),'YLim',[0 0.7])
set(axu(1),'YTick',0:.1:.7)
set(axu(2),'YLim',[210 350])
set(axu(2),'YTick',210:10:350)
set(axu(2),'XLim',[dates_LC(50) dates_LC(end-100)])
legend('Velocity Magnitude','Mean Tension')
ylabel(axu(1),'Current Magnitude,m/s')
ylabel(axu(2),'Mean Tension,lbs')
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
title('Mean Tension on Velocity magnitude')

%% Finding High Tensions with LOW WAVES

low_wavst = 19;
low_wavst_hr = 00;
low_wavend = 23;
low_wavend_hr = 00;
dat1 = datenum(2016,01,low_wavst,low_wavst_hr,00,00);
dat2 = datenum(2016,01,low_wavend,low_wavend_hr,00,00);

[~, in] = min(abs(dates_cur-dat1));
[~, in2] = min(abs(dates_cur-dat2));

[~, inwv] = min(abs(dates_waves-dat1));
[~, in2wv] = min(abs(dates_waves-dat2));

[~, inLC] = min(abs(dates_LC-dat1));
[~, in2LC] = min(abs(dates_LC-dat2));

figure; %Plot current effects during LOW WAVE EVENTS 1/18 09 to 1/24 09 
subplot(311),[axx,~,~] = plotyy(dates_cur(in:in2+700),velmag(in:in2+700),dates_cur(in:(in2+700)),w_lev(in:(in2+700)));
datetickzoom
ylabel('Velocity Mag (m/s)')
set(axx(1),'YTick',(0:.1:1))
set(axx(2),'YLim',[7 11])
set(axx(2),'YTick',(7:.5:11))
legend('Current Magnitude','Water Level','Location','NorthWest')
subplot(312),[a_x,~,~] = plotyy(dates_cur(in:in2+700),veldir(in:in2+700),dates_cur(in:(in2+700)),w_lev(in:(in2+700)));
datetickzoom
set(a_x(2),'YLim',[7 11])
set(a_x(2),'YTick',[7:.5:11])
legend('Current Direction','Water Level','Location','NorthWest')
ylabel('Velocity Dir (degrees)')
subplot(313),[a_x2,~,~] = plotyy(dates_LC(inLC:in2LC+20),LC_filt(inLC:in2LC+20),dates_LC(inLC:in2LC+20),LC_stddev(inLC:in2LC+20));
legend('LC2 Mean','LC2 Std Dev','Location','NorthWest')
datetickzoom
legend('Show')
linkaxes(findall(gcf,'type','axes'), 'x');
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
suptitle('Low Wave Event 1/18 to 1/24 Overlay')

m = 1;
for ii = in:(in2+700)
    
    if abs(w_lev(ii)-8.008) < .003
        date_mllw(m) = dates_cur(ii);
        low_wat_mag(m) = velmag(ii);
        low_wat_dir(m) = veldir(ii);
        m=m+1;
    end
    
end

%% Standardization and Cross-Correlation

ttol=datenum(0,0,0,0,0,10); %Within 30 sec
m = 1;
for ii = 1:length(dates_LC)
    
    for jj = 1:length(dates_cur)
        
        if abs(dates_cur(jj)-dates_LC(ii))<ttol && (jj <= (length(dates_cur)-9) && jj > 10)   
            current_cross(m) = mean(velmag((jj-9):(jj+9)));
            wL_cross(m) = mean(w_lev((jj-9):(jj+9)));
            m = m+1;
                
        elseif abs(dates_cur(jj)-dates_LC(ii))<ttol && (jj >= (length(dates_cur)-9) || jj < 10)
            current_cross(m) = velmag(jj);
            wL_cross(m) = w_lev(jj);
            m = m+1;        
     
        end
        
    end    
    
end

%NORMALIZE IT
stan_LC = (LC_filt(1:250)-min(LC_filt(1:250)))/(max(LC_filt(1:250))-min(LC_filt(1:250)));
stan_LCsd = (LC_stddev(1:250)-min(LC_stddev(1:250)))/(max(LC_stddev(1:250))-min(LC_stddev(1:250)));
stan_WL = (wL_cross(1:250)-min(wL_cross(1:250)))/(max(wL_cross(1:250))-min(wL_cross(1:250)));
stan_curmag = (current_cross(1:250)-min(current_cross(1:250)))/(max(current_cross(1:250))-min(current_cross(1:250))); 

%Try checking the cross correlation of WL and Tension
lag_t = 6;
    %TIDE AND MEAN TNESION
crossc_wl_T = cross_corr(0,lag_t,stan_WL,stan_LC);
crossc_T_wl = cross_corr(0,lag_t,stan_LC,stan_WL);
crosscTIDE_TENS = [fliplr(crossc_T_wl), crossc_wl_T(2:end)];

    %CURRENT AND STD OF TEnSION
crossc_cur_stdT = cross_corr(0,lag_t,stan_curmag,stan_LCsd);
crossc_stdT_cur = cross_corr(0,lag_t,stan_LCsd,stan_curmag);
crosscCUR_stdTENS = [fliplr(crossc_stdT_cur), crossc_cur_stdT(2:end)];

% CURRENT AND MEAN TENSION
crossc_cur_T = cross_corr(0,lag_t,stan_curmag,stan_LC);
crossc_T_cur = cross_corr(0,lag_t,stan_LC,stan_curmag);
crosscCUR_TENS = [fliplr(crossc_T_cur), crossc_cur_T(2:end)];

l1 = 37;
l2 = 53;

% Calculating the number of independent variables
[N_wlT, ~] = Num_ind_var(l1,l2,stan_WL,stan_LC);
[N_cur_T, ~] = Num_ind_var(l1,l2,stan_curmag,stan_LC);
[N_cur_stdT, ~] = Num_ind_var(l1,l2,stan_curmag,stan_LCsd);

% Chi squared value
    X_a = 3.84; % 1 DOF 95 % conf
    %X_a = 2.71; % 1 DOF 90%  conf

    %Cross Correlation
confwlT = (X_a/N_wlT)^.5*ones(1,length(-lag_t:1:lag_t));
confcur_T = (X_a/N_cur_T)^.5*ones(1,length(-lag_t:1:lag_t));
confcur_stdT = (X_a/N_cur_stdT)^.5*ones(1,length(-lag_t:1:lag_t));

if lag_t == 0
    disp(['The Correlation of Water Level to Tension is',crosscTIDE_TENS])
else
    figure,plot((-lag_t:1:lag_t),crosscTIDE_TENS,(-lag_t:1:lag_t),confwlT,(-lag_t:1:lag_t),-confwlT)
    set(gca,'XTick',-lag_t:1:lag_t)
    title('Cross Correlation of Tension and Water Level')
    xlabel('lag (hrs)')
    ylabel('Correlation')
    text(2,0.5,'Negative Lags means')
    text(2,0.45,'         Tension LEADS')
    text(2,0.4,'               Water Level')
end


if lag_t == 0
    disp(['The Correlation of Vel Mag to Tension is',crosscCUR_TENS])
else
figure,plot((-lag_t:1:lag_t),crosscCUR_TENS,(-lag_t:1:lag_t),confcur_T,(-lag_t:1:lag_t),-confcur_T)
set(gca,'XTick',-lag_t:1:lag_t)
title('Cross Correlation of Tension and Current Magnitude')
xlabel('lag (hrs)')
ylabel('Correlation')
text(2,0.6,'Negative Lags means')
text(2,0.55,'         Tension LEADS')
text(2,0.5,'               Vel mag')
end

%FOR THESIS

    figure,
    subplot(2,1,1)
    plot((-lag_t:1:lag_t),crosscCUR_TENS,'LineWidth',1.5)
set(gca,'XTick',-lag_t:1:lag_t)
title('Cross Correlation of Tension and Current Magnitude','FontSize',22)
ylim([-0.5 0.7])
ylabel('Correlation')
set(gca,'FontSize',18)
text(4,0.55,'Negative Lags means')
text(4,0.475,'         Tension LEADS')
text(4,0.4,'              Velocity mag')
    subplot(2,1,2),
    plot((-lag_t:1:lag_t),crosscTIDE_TENS,'LineWidth',1.5)
    set(gca,'XTick',-lag_t:1:lag_t)
    title('Cross Correlation of Tension and Water Level','FontSize',22)
    xlabel('Lag Time [Hrs]')
    ylim([-0.5 0.7])
    ylabel('Correlation')
    set(gca,'FontSize',18)
    text(4,0.5,'Negative Lags means')
    text(4,0.425,'         Tension LEADS')
    text(4,0.35,'               Water Level')

%Standard Deviation and Current
figure,
    plot((-lag_t:1:lag_t),crosscCUR_stdTENS,'LineWidth',1.5),hold on
    plot((-lag_t:1:lag_t),confcur_stdT,(-lag_t:1:lag_t),-confcur_stdT)
set(gca,'XTick',-lag_t:1:lag_t)
set(gca,'XLim',[-lag_t lag_t])
title('Cross Correlation of Tension Standard Deviation and and Current Magnitude','FontSize',22)
ylim([-.5 0.7])
xlabel('Date [mm/dd]')
ylabel('Correlation')
set(gca,'FontSize',18)    
    
%% Cross Correlation betwen H_13 and standard deviation of tension


ttolw=datenum(0,0,0,0,0,10); %Within 30 sec
mw = 1;
for ii = 1:length(dates_LC)
    
    for jj = 1:length(dates_waves)
        
        if abs(dates_waves(jj)-dates_LC(ii))<ttolw
            wave_cross(mw) = H13(jj);
            date_wv(mw) = dates_waves(jj);
            wvspec_cross(mw) = .125*1027*9.8*trapz(Waves_Spectra(jj,:),2);
            mw = mw+1;
        end
    end    
    
end

%NORMALIZE IT

stan_LCwv = (LC_filt-min(LC_filt))/(max(LC_filt)-min(LC_filt));
stan_wv = (wave_cross-min(wave_cross))/(max(wave_cross)-min(wave_cross)); 
stan_LCstd = (LC_stddev-min(LC_stddev))/(max(LC_stddev)-min(LC_stddev));
stan_wvspec = (wvspec_cross-min(wvspec_cross))/(max(wvspec_cross)-min(wvspec_cross)); 

%Try checking the cross correlation of WL and Tension
lag_t = 6;
crossc_T_wvstd = cross_corr(0,lag_t,stan_wv,stan_LCstd);
crossc_wv_Tstd = cross_corr(0,lag_t,stan_LCstd,stan_wv);
crosscWAVE_TENSstd = [fliplr(crossc_T_wvstd), crossc_wv_Tstd(2:end)];

crossc_T_wv = cross_corr(0,lag_t,stan_wv,stan_LCwv);
crossc_wv_T = cross_corr(0,lag_t,stan_LCwv,stan_wv);
crosscWAVE_TENS = [fliplr(crossc_T_wv), crossc_wv_T(2:end)];

crossc_T_wvsp = cross_corr(0,lag_t,stan_wvspec,stan_LCstd);
crossc_wvsp_T = cross_corr(0,lag_t,stan_LCstd,stan_wvspec);
crosscWAVEsp_TENS = [fliplr(crossc_T_wvsp), crossc_wvsp_T(2:end)];

%Confidence Intervals
l1 = 37;
l2 = 53;

% Calculating the number of independent variables
[N_wvLCstd, ~] = Num_ind_var(l1,l2,stan_wv,stan_LCstd);
[N_wvLCmn, ~] = Num_ind_var(l1,l2,stan_wv,stan_LCwv);
[N_wvLCwsp, ~] = Num_ind_var(l1,l2,stan_wvspec,stan_LCstd);

% Chi squared value
    X_a = 3.84; % 1 DOF 95 % conf
    
    %Cross Correlation
confwvLCstd = (X_a/N_wvLCstd)^.5*ones(1,length(-lag_t:1:lag_t));
confwvLCmn = (X_a/N_wvLCmn)^.5*ones(1,length(-lag_t:1:lag_t));
confwvLCwsp = (X_a/N_wvLCwsp)^.5*ones(1,length(-lag_t:1:lag_t));

%Standard Deviation and Sig Wave Height and Spectra
    figure,
subplot(2,1,1),plot((-lag_t:1:lag_t),crosscWAVE_TENSstd,'LineWidth',1.5),hold on
plot((-lag_t:1:lag_t),confwvLCwsp,(-lag_t:1:lag_t),-confwvLCwsp)
set(gca,'XTick',-lag_t:1:lag_t)
set(gca,'XLim',[-lag_t lag_t])
title('Cross Correlation of Tension Standard Deviation and Significant Wave Height','FontSize',22)
ylim([0 0.6])
set(gca,'YTick',0:.1:.6)
ylabel('Correlation')
set(gca,'FontSize',18)
grid
subplot(2,1,2)
    plot((-lag_t:1:lag_t),crosscWAVEsp_TENS,'LineWidth',1.5),hold on
    plot((-lag_t:1:lag_t),confwvLCmn,(-lag_t:1:lag_t),-confwvLCmn)
set(gca,'XTick',-lag_t:1:lag_t)
set(gca,'XLim',[-lag_t lag_t])
title('Cross Correlation of Tension Standard Deviation and Wave Energy','FontSize',22)
ylim([0 0.6])
set(gca,'YTick',0:.1:.6)
xlabel('Lag, [Hrs]')
ylabel('Correlation')
set(gca,'FontSize',18)
grid

%Mean Tension and Sig Wave Height
figure,
    plot((-lag_t:1:lag_t),crosscWAVE_TENS,'LineWidth',1.5),hold on
    %plot((-lag_t:1:lag_t),confwvLCmn,(-lag_t:1:lag_t),-confwvLCmn)
set(gca,'XTick',-lag_t:1:lag_t)
title('Cross Correlation of Mean Tension and Significant Wave Height','FontSize',28)
ylim([0 0.3])
ylabel('Correlation')
xlabel('Lag, [Hrs]')
set(gca,'FontSize',22)
grid
hold off

%% Plot of wave energy flux directly overlayed with mean tension

figure,[axu,p1,p2] = plotyy(date_wv,stan_wvspec,dates_LC,stan_LCstd);
datetickzoom('x','keepticks')
set(axu(1),'XLim',[dates_LC(1) dates_LC(end)])
set(axu(1),'FontSize',22)
set(p1,'LineWidth',2)
set(p2,'LineStyle','--')
%set(p2,'LineWidth',1.5)
set(axu(2),'YLim',[0 1])
set(axu(2),'YTick',0:.2:1)
set(p2,'LineWidth',1.2)
set(axu(2),'XLim',[dates_LC(1) dates_LC(end)])
set(axu(2),'FontSize',22)
legend('  Wave Energy Density','  \sigma_{Tension}','Location','NorthWest')
ylabel(axu(1),'Wave Energy per Area [J/m^2]')
ylabel(axu(2),'Tension Standard Deviation [lbs]')
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
title('Tension Standard Deviation and Wave Energy Density','FontSize',28)
xlabel('Date [mm/dd]')    


%Plot of signfinicant wave height overlayed with mean tension
figure,
[axu,p1,p2] = plotyy(dates_waves,H13,dates_LC,LC_stddev);
datetickzoom('x','keepticks')
set(axu(1),'XLim',[dates_LC(1) dates_LC(end)])
set(axu(1),'YLim',[0 1])
set(axu(1),'FontSize',22)
set(p1,'LineStyle','-')
set(p1,'LineWidth',1.5)
set(axu(1),'YTick',0:.2:1)
%set(axu(2),'YLim',[0 1])
%set(axu(2),'YTick',0:.2:1)
set(p2,'LineWidth',1.2)
set(p2,'LineStyle','--')
set(axu(2),'XLim',[dates_LC(1) dates_LC(end)])
set(axu(2),'FontSize',22)
legend('  H_{1/3}','  \sigma_{Tension}','Location','NorthWest')
ylabel(axu(1),'Signifcant Wave Height, [m]')
ylabel(axu(2),'Std. Dev Tension. [lb]')
xlabel('Date, [mm/dd]')   
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
title('Tension Standard Deviation and Significant Wave Height','FontSize',28)



figure,
[axu,p1,p2] = plotyy(dates_waves,H13,dates_LC,LC_filt);
datetickzoom('x','keepticks')
ylabel(axu(1),'Significant Wave Height, [m]')
ylabel(axu(2),'Mean Tension, [lb]')
xlabel('Date, [mm/dd]')
set(axu(1),'XLim',[dates_LC(1) dates_LC(end)])
set(axu(1),'YLim',[0 1])
set(axu(1),'FontSize',22)
set(p1,'LineStyle','-')
set(p1,'LineWidth',1.5)
set(axu(1),'YTick',0:.2:1)%
set(axu(2),'YLim',[215 280])
set(axu(2),'YTick',[210:10:280])
set(p2,'LineWidth',1.2)
set(p2,'LineStyle','--')
legend('  H_{1/3}',' Mean Tension','Location','NorthWest')
set(axu(2),'XLim',[dates_LC(1) dates_LC(end)])
set(axu(2),'FontSize',22)
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
title('Mean Tension and Significant Wave Height','FontSize',28)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Relative Times and Values for Channel

chan_cur = importdata('tide_curjan16.xlsx');
chan_date = datenum(datetime(chan_cur.textdata));
chan_vel = chan_cur.data;

chan_cur2 = importdata('tide2_curjan16.xlsx');
chan_date2 = datenum(datetime(chan_cur2.textdata));
chan_vel2 = chan_cur2.data;

i = 1;
g = 1;
for j = 1:length(chan_vel)

    if  chan_vel(j) > 0
    chn_date(i) = chan_date(j);
    chn_vel(i) = chan_vel(j);
    i=i+1;
    end
    
    if chan_vel2(j) > 0
    chn_date2(g) = chan_date2(j);
    chn_vel2(g) = chan_vel2(j);
    g=g+1;
    end
    
end

figure,plot(chn_date,chn_vel,'o'),datetickzoom

ttc = datenum(0,0,0,0,0,1);

k = 1;
for i = 1:length(velmag)
    
    for j = 1:length(chn_date)
        
        if abs(dates_cur(i)-chn_date(j)) <= ttc
           
           peak_cur(k) = max(velmag(i-50:i+50));
           [a(k) b(k)] = ind2sub(size(velmag),find(velmag==peak_cur(k)));
           k = k+1;                   
        end
        
    end
    
end

figure,plot(chn_date2,chn_vel2/100,'o','MarkerSize',10,'LineWidth',1.5),hold on
plot(dates_cur(b),peak_cur,'x','MarkerSize',10,'LineWidth',1.5)
dcm_obj = datacursormode(figure(gcf));
set(dcm_obj,'UpdateFcn',@datestrcurs)
datetickzoom('x','keepticks')
xlabel('Date [mm/dd]')
ylabel('Velocity Magnitude [m/s]')
title('Maximum Flood Currents in Channel and at Site','FontSize',24)
legend('Channel','Site')
set(gca,'FontSize',20)
%xlim([chn_date(1)-0.25 chn_date(end)+.25])


Cur_ratio = peak_cur./(chn_vel/100);
Cur_ratio2 = peak_cur./(chn_vel2/100);
mean(Cur_ratio2)

adcp_mac_cur_dates = datevec(dates_cur(b));
tide_cur_dates = datevec(chn_date);
tide_cur_dates2 = datevec(chn_date2);

time_dif = etime(adcp_mac_cur_dates,tide_cur_dates);
time_dif2 = etime(adcp_mac_cur_dates,tide_cur_dates2);
%This is the time elapsed from max cur in cahnnel to max cur seen by ADCP

time_dif_min = time_dif./60
time_dif_min2 = time_dif2./60

Mean_time_dif = mean(time_dif_min)
Mean_time_dif2 = mean(time_dif_min2)
