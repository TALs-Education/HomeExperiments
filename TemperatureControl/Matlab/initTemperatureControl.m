%% initTemperatureControl
clc
clear
%% Numerical parameters (all SI units)
T_amb  = 25;           % °C ambient (for plotting offset)
C_L    = 28.8;         % J/K   local thermal mass  (heater+heatsink)
C_R    = 0.85;         % J/K   remote pad thermal mass
R_h    = 32;           % K/W   heatsink → ambient
R_cond = 71;           % K/W   copper trace
R_conv = 173;          % K/W   remote pad convection
P_step = 1;            % W     step magnitude (change later if needed)

%% State-space matrices  (states = ΔT_L, ΔT_R)
A = [-1/(C_L*R_h)-1/(C_L*R_cond),   1/(C_L*R_cond);
      1/(C_R*R_cond),             -1/(C_R*R_cond)-1/(C_R*R_conv)];
B = [1/C_L; 0];
C = eye(2);
D = [0; 0];