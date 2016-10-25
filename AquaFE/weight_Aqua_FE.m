
%% Aqua-FE INPUTS: ALL MEASUREMENTS should be in or converted to SI

clc,clear,close all
%expected draft is 
Draft_exp = 0.455*.3048;  %in m

% IN GENERAL I NEED TO DECREASE VOLUME INCREASE DENSITY TO MATCH MASS 
% BUT INCREASE/MATCH DRAFT

% 1) check volumes and masses of pipes and align them
%       ~3,469.4 LBS 10 IN PIPE @ 150.8858 ft^3 (4.2726 m^3)
    
% 2) do the same for CROSS plates, floats

%       ~523.7190 lbs @ 15404 in^3 (8.9141 ft^3) plates
%       ~651 lbs @ 8 ft^3 floates

% 3) mussels lines, using in water weight of out of water?


% 4) figure out decking and stanchions add weight from real raft
%    4.5) figure out how to add it to aquafe raft
%            ~ for stanchions make a stanchions element, make heavier
%              plug in, instead of plate element



%% Real Raft information
    
    %Pipes
               %lbs/ft / A
eff_denpipe = 13.14/((5.375/12)^2*pi); %lbs/ft^3

actual_den10in = eff_denpipe*(1/2.20462)*(1/.3048)^3; %conversion to kg/m^3

tot_L10 = .3048*(209.3889 + 4*2.5 + 8*30/12);
d_10_m = (10.75*.0254);

actual_WPA_pipe = tot_L10*d_10_m; %m^2

actual_V_pipe = 150.8858*.3048^3; %Comes from actual lengths from SW model

actual_mass_pipe = 3416/2.20462; %kg
    
    %mussels
actual_mass_muss = 3663.4*3/2.20462; %lbs to kg

%actual_ROUGH_mass_net = .005*

    
    %mussel ballast
actual_mass_mussB = 2/2.20462; %for a single ballast
act_mass_mussB_tot = actual_mass_mussB*110;

actual_CSA_mussB = 4.125*2.625*.0254^2;
actual_V_mussB = (2.625/2*.0254)^2*pi*(4.125*.0254);
                      %Area               Height  

%% HDPE Pipe Pieces Aqua_FE

denH = actual_den10in; %We want to match the pipe density
denV = denH;

A_H = .0173; %This area was calculated 4/1/16 it comes from using the actual density
A_V = A_H;   % of the pipe, the actual weight of all the 10 in pipe, and the total
             % length of all the AquaFE pipe elements.
             % I apply the density to each element that makes up the pipe
             % and find the corresponding area that, when multiplied by the
             % total length gives the same mass as the actual model roughly.

N_H = 318;
N_V = 192;

%Horizontal Pipe sections
LHlongout = 12.6111*3*2; % 3 layers deep and 2 long outer sides
LHshortout = 6.7056*3*2; % 3 layers deep and 2 short outer sides
LHpen1 = 5.39115*3*4; % three b/c 3 layers 4 b/c 4 pen1 sides
LHpen2 = 5.4864*3*4; % three b/c 3 layers 4 b/c 4 pen2 sides

L_Htot = LHlongout+LHshortout+LHpen1+LHpen2; %m

V_Htot = L_Htot*A_H;

mass_Htot = V_Htot*denH; %This mass is in kg

%Vertical Pipe sections
L_V1 = 0.135; %m

L_Vtot = L_V1*N_V;
V_Vtot = L_Vtot*A_V;

mass_Vtot = V_Vtot*denV; % in kg

d_horiz = sqrt(4*A_H/pi);
WPA_pipe = L_Htot/3*d_horiz; %Water plane area of ALL the Pipe elements

%Total of the 10 in pipes
V_pipe_tot = (V_Vtot+V_Htot); %m^3


%% New Calcs to match things 4/1/16

Good_A_pipe = actual_mass_pipe/((L_Vtot+L_Htot)*actual_den10in);

Good_mass_pipes = denH*(L_Vtot+L_Htot)*Good_A_pipe; %kgs

%NOTE, THIS "NEW" area has already been plugged in!


%% HDPE Plates

denPlate = 943.5; %kg/m^3 %This matches actual plate density

N_plate = 126;
L_plate = .2032;

%Mass below is taken from SW model calcs (Weight.m) in lbs
actual_plate_mass = 588.1531/2.20462; %lbs from estimate/2.20462 for kg

Good_A_plate = actual_plate_mass./(denPlate*L_plate*N_plate);

Good_plate_mass = denPlate*L_plate*N_plate*Good_A_plate;

%% plates that account for mass of stanchions

stanch_weight_tot = 1402.5/2.20462;
W_deck = 1961.7/2.20462;

fat_plate_mass = 588.1531/2.20462 + stanch_weight_tot + W_deck;

Good_A_fatplate = fat_plate_mass./(denPlate*L_plate*N_plate);

Good_fatplate_mass = denPlate*L_plate*N_plate*Good_A_fatplate;


%% Floats

actual_float_V = 8*.3048^3; %m^3
actual_float_den = 31*(1/2.20462)/actual_float_V;
                 %31 lbs/2.2->kg/m^3 for dens
                 
actual_WPA_float = (4*.3048*2*.3048);
A_float_div4 = actual_WPA_float/4;

denFloat = actual_float_den; %this matches the actual density of the floats
A_float = .185; %this is the old A_float

N_float = 84;
L_float = .3048; %1 foot in depth (which matches the true depth)

V_floattot_b4 = L_float*N_float*A_float;

WPA_float_b4 = (A_float/pi)^.5*2;

mass_float_tot = denFloat*V_floattot_b4;

Good_A_float = 21*actual_float_V/(L_float*N_float);

Good_V_float = L_float*N_float*Good_A_float;

Good_mass_float = denFloat*Good_A_float*L_float*N_float;

%% Mussel Lines & their ballast

N_full_ropes = 22;
N_actual_ropes = 110; %5 ropes per 1 consolidated

denMuss = 1794; %This is the actual density of the mussel ropes
A_muss = .02356; 

    %ROPES SHOULD CHANGE DIAMETER TO ~ 2in ropes

%One rope represents 5 consolidated ones    
N_muss = 264;
L_muss = .3048;

V_musstot = N_muss*L_muss*A_muss;

mass_muss_tot = V_musstot*denMuss;

%Good_V_muss = Good_A_muss*L_muss*N_muss;


%Empty Mussels

act_Lmuss = 3.6576;
act_mass_empty_muss1 = 5/2.20462; %Rough saturated mass of 1 line
act_A_empt_muss = pi*((1.05)*.0254)^2/4;
act_V_muss_empt1 = act_Lmuss*act_A_empt_muss;

act_den_muss_emp = act_mass_empty_muss1/act_V_muss_empt1;

mass_empt_cons = act_mass_empty_muss1*5;

mass_empt_cons_tot = mass_empt_cons*22;

Good_A_empt_cons = mass_empt_cons/(3.6576*act_den_muss_emp);


%Mussel Ballast

N_mussB = 22;
L_mussB = .104775; %4.125" tall

A_mussBfive = (5*(actual_CSA_mussB/pi)^.5)^2*pi; %5x radius

den_concrete = 2400; %actual density of concrete

L_test_mussB = 5*actual_mass_mussB/(den_concrete*.1750);

Good_V_mussB = N_mussB*A_mussBfive*L_test_mussB;

Good_mass_mussB_tot = N_mussB*den_concrete*L_test_mussB*.1750;



%% Nets & their ballast (Net coeff are 27.14 and 1.2)
denNet = 1200;
A_net = 4.910E-6;

%N_net = 432;
L_ratio = 22.07;

%there are(28v + 30h) = 58 elements in a panel short side panels are LnetV x LnetV
% long side panels are L_netH x L_netV

L_netV = .9144;
L_netH = .9144;
L_netH2 = .898523;

L_net_tot = 2*(L_netV*28+L_netH*30)+2*(L_netV*20+L_netH2*30);
                %Short Side Panels     %Long Side Panels
V_nettot = L_ratio*L_net_tot*A_net;

mass_net_tot = V_nettot*denNet;

denNetB = 2168000;
A_netB = 1.287E-5;

N_netB = 16;
L_netB = .4064;

V_netBtot = N_netB*L_netB*A_netB;

mass_netB_tot = V_netBtot*denNetB;

%% Totals 

Total_mass_emp_muss = (Good_mass_pipes+Good_fatplate_mass+Good_mass_float ...
                    +mass_empt_cons_tot+Good_mass_mussB_tot+mass_net_tot+mass_netB_tot)

Total_mass = (Good_mass_pipes+Good_plate_mass+Good_mass_float ...
                    +mass_muss_tot+Good_mass_mussB_tot+mass_net_tot+mass_netB_tot)                
                
Total_Vol = Good_A_pipe*(L_Htot+L_Vtot)+Good_A_plate*N_plate*L_plate+ ...
            Good_V_float+V_musstot+Good_V_mussB+V_nettot+V_netBtot
                
                
Total_weight_lbs = 2.20462*Total_mass


Total_weight_emp_lbs = 2.20462*Total_mass_emp_muss
     
                
if Total_weight_emp_lbs > 9804.7
    disp(' ')
    disp(['Model is ',num2str(Total_weight_emp_lbs-9804.7),' lbs over your actual weight estimate'])
    disp(' ')
end

if Total_weight_emp_lbs < 9804.7
    disp(' ')
    disp(['Model is ',num2str(-Total_weight_emp_lbs+9804.7),' lbs under your actual weight estimate'])
    disp(' ')
end

    
%% USING WPA IS CORRECT

N_V_new = 190-56;

Location_toDrop_Vert = (N_V - N_V_new)/2;
    %can drop 12 from corners
    % or 28 loactions from inner float sections (insert stiffeners instead)

Av_wp = 2*actual_WPA_pipe/N_V;

V_Hw = actual_V_pipe - Av_wp*L_V1*N_V_new;

A_Hwp = V_Hw/L_Htot;

%% Polysteel Calcs

% 2 in polysteel weight per foot = 790 lbs per 1200 ft

A_poly2 = 2.027e-3; %2 inch diam in m^2

poly2_den = (790*.453592)/(1200*.3048*A_poly2);

%% Chain 

ch_lin_den = 1.45*.453592/.3048; %kg/m

A_1_link = (.0254*(.375/2))^2*pi; %m^2
L_1_link = (4*.0254);

link_per_m = ceil(3.2808/(1.22/12)); %number of links in a meter

V_per_m = A_1_link*link_per_m*L_1_link;

actual_den_chain = ch_lin_den/V_per_m;

L_ch_el = .333375;

Good_A_chain = (ch_lin_den*L_ch_el)/(actual_den_chain*L_ch_el);



                
