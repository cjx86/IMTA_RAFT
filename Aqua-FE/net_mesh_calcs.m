%Net Element Calculcation
clc, clear,close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Actual Values from the real world nets
 
d = 0.08/12*.3048; %twine diameter in meters (0.08")
half_mesh = 1.125/12*.3048; %half mesh size of nets (1.125")
Area = 13*.3048*16*.3048; %Net panel outside area (12' x 15')
Area_bottom = 16*16*.3048^2;

L_tactual = 2*Area/half_mesh;

L_tactual_bot = 2*Area_bottom/half_mesh;

S = 2*d/half_mesh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Values from AquaFE model. The depth and width of the net panels
  %are different because of how the gemoetry of the raft is technically
  %different becaue of the use of trusses as the structure.
 
a = 3.6576; %Net Depth in m from AquaFE
b = 5.39115; %Net Panel width in m from AquaFE
el_spc = 0.9144; %Element spacing in AquaFE nets

%Calculation of the total length of the elements in the model
L_els = (6*5+4*7)*el_spc;

%Calculation of the length ratio of the real world and the model
ratio = L_tactual/L_els;


%Calculation of the normal drag coefficient using d and half mesh
% These are with empirical solutions


C_d = 1 + 2.73*(d/half_mesh) + 3.12*(d/half_mesh)^2;

C_d2 = 1 - 1.24*S + 13.7*S^2;

%There ARE alternatives that are based on reynolds number or current
%velocity (SEE BELOW)


%% Normal Drag Coefficient

max_v = .65; %m/s
min_v = .03;

d_el = .0026; %m

minRe_n = 1027*d_el*min_v/.00162;

maxRe_n = 1027*d_el*max_v/.00162;

Cn = 1.1 + 4*maxRe_n^(-.5);

