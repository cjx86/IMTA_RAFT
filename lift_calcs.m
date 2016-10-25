clc, clear, close all

W_actual = 5.8775e+03;

L_full_str = 30;

L_to_lift = L_full_str - 9.1;
    %30' strap used to the farthest ends of the raft
    
d_short_side = 15 + 2*10.75/12 + 2;
    %distance to the short end of the raft

    d_long_side = 132.77/12;

% Full lift calculation

To = W_actual/(4*cosd(44.09));
To2 = W_actual/(4*cosd(38.91));
To3 = W_actual/(4*cosd(36));
    
COG_wide = 65.56/12;

COG_nar = (271-226)/12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WIDE

W_wide = 4543.61-42.5*20;

L1 = 132.77/12;
L2 = d_short_side;

d1 = 42.81/12;
d2 = 33.62/12;
d3 = d2;

h_wide = sqrt(L_to_lift^2 - L2^2 - d2^2);
    %Height of the lift poin in the wide lift.
        %Dictated by use of 30' strap from  chosen,short-end lift point
    
L_Lside_wide = sqrt(h_wide^2 + L1^2 + d1^2);
    %Length from lift point to chosen long side lift location    
    
%Calculation  to reuse a 30' strap for the long-side lift point   
num_wraps_Lside_wide = (L_full_str - L_Lside_wide)/9.1;
    circ_pipe = pi*10.75/12;
    dum = num2str(num_wraps_Lside_wide);
    single_wr_Lsidewide = (9.1*str2double(dum(2:6)))/circ_pipe;

disp(['Number of Wraps for Long edge Wide lift = ',num2str(single_wr_Lsidewide)])    
disp(' ')    

    %Length of a center strap
L_mid_wide = sqrt(h_wide^2 + d3^2) + 9.1;    

disp(['Length of Extra strap for middle = ',num2str(L_mid_wide),' ft'])

%% Wide Tensions
theta1_W = atan(sqrt(L1^2 + d1^2)/h_wide);
theta2_W = atan(sqrt(L2^2 + d2^2)/h_wide);
theta3_W = atan(d3/h_wide);

gam1_W = atan(d1/L1);
gam2_W = atan(d2/L2);

k_W = [2*cos(theta1_W) 2*cos(theta2_W) cos(theta3_W); ...
       2*sin(theta1_W)*sin(gam1_W) -2*sin(theta2_W)*sin(gam2_W) -sin(theta3_W);...
       2*cos(theta1_W)*(d1) -2*cos(theta2_W)*(d2) -cos(theta3_W)*(d3)];
  
A_W = [W_wide;0;0];

%rref([k_W A_W])

chk_pinvW = A_W./(k_W*pinv(k_W)*A_W);

x_W = pinv(k_W)*A_W;

check_W = W_wide./sum([2*cos(theta1_W);2*cos(theta2_W);cos(theta3_W)].*x_W)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NARROW

W_nar = 3526.45-42.5*12;

L4 = L1;
L5 = d_short_side;

d4 = 22.25/12;
d5 = 56/12;
d6 = d5;

h_nar = sqrt(L_to_lift^2 - L5^2 - d5^2);


L_Lside_nar = sqrt(h_nar^2 + L4^2 + d4^2);

    %Calculation  to reuse a 30' strap 
num_wraps_Lside_nar = (L_full_str - L_Lside_nar)/9.1;
    dum2 = num2str(num_wraps_Lside_nar);
    single_wr_Lside_nar = (9.1*str2double(dum2(2:6)))/circ_pipe;

L_mid_nar = sqrt(h_nar^2 + d6^2) + 9.1; 

%% Narrow Tensions
theta1_N = atan(sqrt(L4^2 + d4^2)/h_nar);
theta2_N = atan(sqrt(L5^2 + d5^2)/h_nar);
theta3_N = atan(d6/h_nar);

gam1_N = atan(d4/L4);
gam2_N = atan(d5/L5);

k_N = [2*cos(theta1_N) 2*cos(theta2_N) cos(theta3_N); ...
    2*sin(theta1_N)*sin(gam1_N) -2*sin(theta2_N)*sin(gam2_N) -sin(theta3_N); ...
    2*cos(theta1_N)*(d4) -2*cos(theta2_N)*(d5) -cos(theta3_N)*(d6)];

A_N = [W_nar;0;0];

chk_pinvN = A_N./(k_N*pinv(k_N)*A_N);

x_N = pinv(k_N)*A_N;

check_N = W_nar./sum([2*cos(theta1_N);2*cos(theta2_N);cos(theta3_N)].*x_N)

