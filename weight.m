%% Mussel Raft Weight Calculations
clc, clear, close all
%% Number of each part
N_10t = 4;
N_10el = 8;
N_billets = 21;
N_plates = 63; %number of cross plates
N_stanchions = 33;

%% Lengths of pipe ALL LENGTHS ARE INITIALLY IN FEET

L_Lin = 458.34/12;
L_Lout = 458.34/12;
L_Sout = 180/12;
L_Sin = 249.5/12;
L_cross = L_Sin;

    %Total length of the 10 in HDPE pipe FEET
L_10 = L_Lout*2 + 13.7275*4 + (10.125/12)*2 + (137.2/12)*2 + (77.5/12)*2 + ...
       (102.45/12)*2 + (42.75/12)*2 + (105.12/12)*2 + (44.51/12)*2;
       
   
    % Total length of nubs for stanchions
L_nub = 22; 


%% Vols of diff plate pieces (in^3) b/c of given density units

    %Volume of cross plates
V_cross = 40.75*8*.75; %in^3

    % Volume of extra HDPE plates
V_ext_plates = 2*(.75*42*18) + 4*(.75*54*5); %in^3
                %Large middle      %Corner plates
   
V_gusset = (N_plates*2+16)*5.41;
    % the 16 is extra gussets for net hanging and corner plates
    % the 2 is because there are 2 gussets for each cross plate
    % 5.41 is the volume of the gusset

%% Weight of single or weight per unit length or density

w_billet = 31; %lbs
lin_den_10 = 13.14; %lbs/ft
w_elbow = 27; %lbs
w_tee = 40*2/3; %40 lbs is at full size One that was ordered is smaller
w_flange = 12+10.5;   

w_45_elbow_6 = 5.5; %lbs
w_wye = 27; %lbs
w_lil_T = 10; %lbs

den_plate = 0.034; %lbs/in^3

%Hardware

w_bolts = .3125 + .173 + .355 + .2348;
w_nut_h = .038;
w_nut_sm = .0162;

w_pipe_sup = 0.97; %lbs
L_den_muss_pipe = 1.6; %lb/ft

v_pipe_sup = .97/.291;

%% Calculate the total weight of each type of material

W_pipe = lin_den_10*L_10 + 2.3*L_nub;

W_tees = N_10t*w_tee;

W_elbows = w_elbow*N_10el;

W_flanges = w_flange*12; 

           %cross plates                add'l support plates     %gusset
            %estimate
W_plates = N_plates*V_cross*den_plate + V_ext_plates*den_plate + V_gusset*den_plate;

W_floats = N_billets*w_billet;

W_wyes = N_stanchions*w_wye;

W_45s = N_stanchions*w_45_elbow_6;

W_little_Ts = N_stanchions*w_lil_T;

W_stanch = W_wyes + W_45s + W_little_Ts;

   A_stanch = pi*((6.625/12)^2-(5.348/12)^2)/4; %lb/ft
   L_stanch1 = 2*(7.5/12) + 18/12   + 36.75/12 + (6.5/12 + 16/12);  
                %45 deg      %45 leg    %up leg    %small T

V_stanch = A_stanch*L_stanch1*N_stanchions;

W_mussels = 3663.44; %this is in lbs and is in-water weight

%W_mussels = 3333.33;

W_muss_dry = W_mussels*3;

Fb_muss = W_muss_dry - W_mussels;

W_empty_muss = 5*110;

      n_boards = 239; %number of deck boards
      l_boards = 26.5; %length of boards (in.)
      lin_den_brd = 13/8; %linear density of deck board 13 lb per 8 ft section
      
W_deck =     23*48     +   (n_boards*l_boards/12)*lin_den_brd;
        %23, 10' 4 x 4 weight      total board length

A_4by = (3.5/12)^2; %ft^3     
A_board = (1/12*5.5/12); %5/4 x 6 in ft^3

V_deck = 23*A_4by*10   +    (n_boards*l_boards/12)*A_board;
        %4 x 4 Volume         %Deck Board Volume

W_hardware = (6*w_bolts)*21 + (12*(w_nut_h+w_nut_sm))*21;

W_muss_pipe = 42*w_pipe_sup + 120*L_den_muss_pipe;

%% Various Total Weights and Res. Buoyancy

    %Total weight and volume of the 10in HDPE
W_10in_HDPE = lin_den_10*L_10 + W_tees + W_flanges + W_elbows;

V_10in_HDPE = pi*(5.375/12)^2*(L_10 + 4*3.5 + 8*30/12); %ft^3
            %est lengths large fittings ^--------^

W_lift = W_pipe + W_elbows + W_tees + W_flanges + W_plates + W_floats + 23*48;            
            
W_mat = W_pipe + W_elbows + W_tees + W_flanges + W_plates + W_floats + W_wyes + W_45s + W_little_Ts;

W_all_mat = W_mat + W_deck + W_hardware + W_muss_pipe

W_net = 200; %Assuming Net is 200 lbs

W_net_bal = 400;

W_muss_bal = 2*110;

W_total = W_all_mat + W_mussels + W_net + W_net_bal + W_muss_bal;

W_total2 = W_all_mat + W_muss_dry + W_net + W_net_bal + W_muss_bal;

W_total_empty = W_all_mat + W_empty_muss + W_net + W_net_bal + W_muss_bal

Fb_net = .43*100;
Fb_muss_bal = .43*440;
Fb_muss_emp = 550*.43;
Fb_floats = 64*7.95*21;

Fb_vol = Fb_floats + Fb_net + Fb_muss_bal + 64*(V_10in_HDPE + (V_cross+V_ext_plates+V_gusset+42*v_pipe_sup)/12^3 + V_deck + V_stanch);

%Fb_vol_emp = Fb_floats + Fb_net + Fb_muss_emp + Fb_muss_bal + 64*(V_10in_HDPE + (V_cross+V_ext_plates+V_gusset+42*v_pipe_sup)/12^3 + V_deck + V_stanch);    
    
Res_Buoy = 100*((Fb_vol - W_total)/Fb_vol)

%Fb_vol2 = Fb_floats + Fb_net + Fb_muss + Fb_muss_bal + 64*(V_10in_HDPE + (V_cross+V_ext_plates+V_gusset+42*v_pipe_sup)/12^3 + V_deck + V_stanch);    
%Res_Buoy2 = 100*((Fb_vol2 - W_total2)/Fb_vol2)

% cost_10 = 16*L_10 + (325*8+250*4) + 71*21;
