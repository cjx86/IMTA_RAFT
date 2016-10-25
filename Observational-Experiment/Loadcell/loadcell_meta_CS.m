% NOTE: TIME SERIES PLOTTING IS ON LINE 190 You must input the month,
% day and the year
 
clear;clc;close all;
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))
addpath(genpath('../'));
addpath(genpath('../functions'));

%% Input desired data to read

% %Load filenames
lcnum='2';

%% Find the data desired in the pwd 

dirName =['..\Loadcell\L' lcnum];
dirData = dir(dirName);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
%COPIED

   %Checking which deployment and loadcell  data is from
    if strcmp(lcnum,'2')
        
    st_id=1; %First file of interest
%     ids=st_id:(st_id+20) %For testing
    ids=[st_id:length(fileList)-1];
    endList=length(fileList)-1;
    
    else 
        disp(['L', lcnum, ' Data from ', dep, ' Deployment is bad'])
        disp(' ')
        disp('Individual record investigation is recommended')
           st_id=1; %First file of interest
           ids=[st_id:length(fileList)-25];
           endList=length(fileList)-25;
        
    end
%         st_id=120+61; %First file of interest
%         ids=[st_id:st_id+1];

%st_id=28; %First file of interest 01/15/16 HR21 UTC; (hr 16 est)     

lc_all = zeros(endList,6000);
lc1_all = zeros(endList,6000);

%%      Looping through the Loadcell Data 
for id=st_id:endList
    
kk=id+1-st_id;
disp(num2str(id))    
%    
% sname = 'files.dat';
% fiz = fopen(sname,'r');
% fname=fscanf(fiz,'%s');
% fclose(fiz)
% %

% 
% fix=fopen(fid,'r');
fname=fileList{id};
% fid=fname(((id-1)*12+1):((id-1)*12+12))
fix=fopen([dirName filesep fname]);
[dat,count] = fread(fix,24000,'uint8');
lc_pre=2.5.*[256*dat(1:4:24000) + dat(2:4:24000)]./32768; 
bat=6.*2.5.*[256*dat(3:4:24000) + dat(4:4:24000)]./32768;
%2.5 is for A/D board
fclose(fix); 
%

% normalize to standard values for 20,000 lb load cell
%  sen= 4258;  % 8,000 lbs/volt * 2.5 v = 20,000 lb full scale
%  For SN 010 and LC2 The zero reading was 0.192 in the lab.

freq=5;
del = 1/freq;
sen=20000./(2.069*18/2); %2.069 mV/V. Excitement voltage: 18V/2?
zer1=0.192;
zer = 0;
lc_a=sen.*(lc_pre-zer);
%lc1 = sen.*(lc_pre-zer1);

month(kk)=str2num(fname(3:4));
day=str2num(fname(5:6));
hour=str2num(fname(7:8));

yr = 2016;

yd(kk)=datenum(yr,month(kk),day,hour,00,00);
%
mn_lb(kk)=mean(lc_a);
std_lb(kk)=std(lc_a);
bt(kk)=mean(bat);
%

%% Calculate spectra
padyes=0;
hannme=1;
n_ens=1;
alph=0.05;
Ffreq=freq/(floor(length(lc_a)/n_ens)); %Fourier frequency spacing
n_bands=round(0.5/Ffreq/64); %Band average to match wave data (64 freqs between 0 and 0.5 Hz)
num_std=100; %Standard deviations to filter--high # means don't filter
pass=1;

%Set number of indicies over which to find a running mean. To just demean: mean_bands=0
mean_bands(1)=0; 

%[bandfT S_pT con95T bandwidthT S_areaT(kk) lc_varT(kk)]=spectra2(lc,freq,n_ens,n_bands,padyes,hannme,num_std,pass,mean_bands,alph);
            % Toby's function
                % Band avg
                % Ensemble avg
                % Confidence Intervals
                % Variance
%length(lc)
                
[lc, lc_mean, lc_var(kk), lc_std(kk), N_lc, lc_sq_mean] = stddev_filter(lc_a,num_std,pass);               
%[lc1, lc1_mean, lc1_var(kk), lc1_std, N_lc1, lc1_sq_mean] = stddev_filter(lc1,num_std,pass); 

lc_all(kk,:) = lc;
%lc1_all(kk,:) = lc1;

[bandf, S_p, conf, S_area(kk), DOF] = spectra_cs(lc,del,n_ens,n_bands,alph);
mean_bands(1)=0; 

[bandff, S_pp, conff, S_areaa(kk), DOFF] = spectra_cs_han(lc,del,n_ens,n_bands,alph);

[bandfreqs_lc, Sj1_, con95, ~, Sjarea_lc(kk), var_lc(kk)]=spectra2(lc_a,freq,n_ens,n_bands,padyes,hannme,num_std,pass,mean_bands,alph);

bandwidth=1/N_lc/del*n_bands*(n_ens+1/n_bands);
               
%Plot
% semilogy(bandfreqs(2:length(bandfreqs)),Sj1side(2:length(Sj1side)),'cy')
cutoff=0.5; %Hz. Truncate spectrum at 0.5 Hz.
cutkk=find(bandf>cutoff,1); 
bandf=bandf(1:cutkk);
bandf_lc=bandfreqs_lc(1:cutkk);
Sp_lc(kk,:)=S_p(1:cutkk);
Spp_lc(kk,:) = S_pp(1:cutkk);
Sj1_lc(kk,:) = Sj1_(1:cutkk);

% figure
% plot(bandfreqs_lc(2:end),Sj1side_cut_lc(kk,2:end))
% hold on
% 
% % xrange_s=0.5; %max(staff_bandfreqs);
% % yrange_s=0.3+0.1;
% title('One-sided Sample Spectral Density of load cell tension')
% xlabel('Frequency')
% ylabel('Spectral Density, lb^2/Hz')
% clear param
% param(1)={['Ensembles: ' num2str(ensembles) ' ']};
% param(2)={['Bands: ' num2str(bands) ' ']};
% text(max(bandfreqs_lc)-1,max(Sj1side_lc(kk,:))-0.1,param,'HorizontalAlignment','right','VerticalAlignment','top')
% xlim([0 xrange_s])
% % ylim([0 yrange_s])
% %plot confidence
% line([0.4 0.4],[0.1 0.1+con95])
% %Find variance by area
% hold on
% %line([0.1 0.1+staff_bandwidth],[0.1 0.1])
% legend('Spectrum', '95% Confidence','Bandwidth','Location','NorthEast')

% conf_cut=con95(1:cutkk,:);
% %% Plot
% % semilogy(bandfreqs(2:length(bandfreqs)),Sj1side(2:length(Sj1side)),'cy')
% figure
% plot(bandfreqs(2:length( bandfreqs)), Sj1side(2:length(Sj1side)))
% hold on
% 
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

%
end

%% Plot the time series of the desired day MUST INPUT DATA HERE

    %INPUT MONTH DAY AND HOUR YOU WANT TO PLOT
mon_plot = '01';
day_plot = '16';
hour_plot = '08'; % 8 reflects UTC time, it is 3 am EST
    %Compile inputs into one string
dh_plot = ['L', lcnum mon_plot day_plot hour_plot];

    %Find the index of the file you want to plot
index_p = strfind(fileList, dh_plot);
index_plot = find(not(cellfun('isempty', index_p)));
    %Use the found index to choose the data set from the entire 
lc_plot = lc_all(index_plot,:);
tim=(1:6000)./freq;

    %High-pass filter,check for low-freq oscillations (Savitzky-Golay Filter)
T_lp=16; %s. Period for low-pass filter
    %Applying the actual filter
lc1_sg=sgolayfilt(lc_plot,5,T_lp*freq+1); %Assuming 5 Hz. Number of frames, f must be odd.

figure,
plot(tim/60,lc_all(index_plot,:),'LineWidth',1.25'),grid,hold on
%plot(tim/60,lc1_sg),grid
xlabel('Time [min]','FontSize',20)
ylabel('Tension [lbs]','FontSize',20)
set(gca,'FontSize',20)
title(['Loadcell ',lcnum, ' Tension 01/',day_plot, '16 ', hour_plot,':00 am'],'FontSize',28)
hold off
%% Calculate and plot the mean and std dev of the loadcell data

for n = 1:length(lc_all(:,1))
    
   [mean_lc_all(n),~,std_lc(n)] = stats(lc_all(n,1100:4800));
          
end

% figure,
% for n = 1:length(lc_all(:,1))
%     subplot(2,1,1)
%     plot(lc_all(n,1100:4800))
%     subplot(2,1,2)
%     plot(lc_all(n,1100:end),'r'),hold off
%     pause
%     
% end


    %Plotting of the Mean Force
figure,plot(yd(28:end),mean_lc_all(28:end),'LineWidth',1.25),hold on
%plot(yd(28:end),mn_lb(28:end))
datetickzoom('x','mm/dd','keepticks')
title(['Mean Tension from Loadcell ',lcnum,' Samples'],'FontSize',28)
ylabel('Force, [lb]')
xlabel('Date [mm/dd]')
set(gca,'FontSize',18)
ylim([210 280])

    %Plotting of the Standard Deviation
figure,plot(yd(28:end),std_lc(28:end),'LineWidth',1.25),hold on
%plot(yd(28:end),lc_std(28:end))
datetickzoom('x','mm/dd','keepticks')
title(['Standard Deviation of Tension from Loadcell ',lcnum,' Samples'],'FontSize',28)
ylabel('Standard Deviation of Force, [lb]')
xlabel('Date [mm/dd]')
set(gca,'FontSize',20)
hold off

figure,plot(bandf,nanmean(Sp_lc(28:end,:),1),'LineWidth',1)
title(['Tension Spectra from Loadcell ',lcnum,' Samples'],'FontSize',22)
ylabel('Standard Deviation Tension Spectra, [lbs^2]')
xlim([0.025 0.5])
xlabel('Date [mm/dd]')
set(gca,'FontSize',18)

%% Save

% [what,when]=max(mn);
% datestr(when)
% if exist(['..\Results\' dep],'file') == 0
%     mkdir('..\Results\', dep)
%     savename=['..\Results\' dep '\LC' lcnum '_meta.mat']
%     save(savename,'yd','mn_lb','std_lb','sen','bt','fileList','bandf','Sp_lc','S_area','lc_var','bandf_lc', 'Sj1_lc')
%     
% else
%     
% savename=['..\Results\' dep '\LC' lcnum '_meta.mat']
% save(savename,'yd','mn_lb','std_lb','sen','bt','fileList','bandf','Sp_lc','S_area','lc_var','bandf_lc', 'Sj1_lc')
% 
% end
% save('..\Results\Submerged\LC3_meta.mat','yd','mn_lb','sd_lb','sen','bt','fileList')



