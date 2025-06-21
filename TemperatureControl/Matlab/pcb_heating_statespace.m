%==========================================================================
% pcb_heating_statespace.m
%
% A fully commented MATLAB script that *derives* and then *implements*
% the 2-node state-space model for the PCB-heating experiment.
%
% ─────────────────────────  Thermal picture  ────────────────────────────
%
%          ┌──────────────┐            ┌──────────────┐
%   P(t) → │  node L      │── R_cond ─▶│  node R      │── R_conv ─▶ Tamb
%          │ (heater)     │            │ (remote pad) │
%          └───┬────┬─────┘            └──────────────┘
%              │R_h │
%              ▼    ▼
%             Tamb  Tamb      (both conduct to ambient)
%
%   States  : x = [ΔT_L ; ΔT_R]   temperatures above ambient (°C or K)
%   Input   : u = P(t)            heater power (W)
%   Outputs : y = x               we read both temperatures
%
%==========================================================================

clc
clear

%% 1.  USER PARAMETERS  (all values in SI units)
Tamb  = 25;      % Ambient temperature  [°C]  (reference for ΔT)
C_L   = 28.8;    % Thermal capacitance, local node (heater+heatsink) [J/K]
C_R   = 0.85;    % Thermal capacitance, remote copper pad            [J/K]
R_h   = 32;      % Thermal resistance, heatsink → ambient            [K/W]
R_cond= 71;      % Thermal resistance, copper trace (L ↔ R)          [K/W]
R_conv= 173;     % Thermal resistance, remote pad → ambient          [K/W]

% Explanation of the symbols:
%   C_L    = m_L * c_p_L   (mass * specific-heat) of heater + heatsink
%   C_R    = m_R * c_p_R   of the small copper/FR-4 pad
%   R_h    = 1 / hA        natural-convection path through the finned Al HS
%   R_cond = L / (k_cu * A) one-dimensional Fourier conduction in Cu trace
%   R_conv = 1 / hA        natural convection from bare remote pad
%
% NOTE:  Replace any of the numbers with refined measurements if available
%        (e.g. from a datasheet or CFD).  Units *must* stay in SI.

%% 2.  DERIVATION OF THE STATE EQUATIONS  (energy balance)
%
%   Node L (heater under heatsink):
%     C_L * d(ΔT_L)/dt =
%         + P(t)                                      (electrical power in)
%         – (ΔT_L) / R_h                              (to ambient)
%         – (ΔT_L – ΔT_R) / R_cond                    (to node R)
%
%   Node R (remote pad):
%     C_R * d(ΔT_R)/dt =
%         + (ΔT_L – ΔT_R) / R_cond                    (from node L)
%         – (ΔT_R) / R_conv                           (to ambient)
%
%  Rearrange into canonical ẋ = A·x + B·u :
%
A = [ -(1/R_h + 1/R_cond)/C_L,    ( 1/R_cond )/C_L ;
       ( 1/R_cond )/C_R,        -(1/R_cond + 1/R_conv)/C_R ];

B = [ 1/C_L ; 0 ];

%  Output matrix – we want both temperatures:
C = eye(2);
D = zeros(2,1);

%% 3.  BUILD THE STATE-SPACE OBJECT
plant = ss(A,B,C,D);

%% 4.  Step simulation  — with 10-s delay and 1-s sampling
Pstep    = 1;      % [W]
t_end    = 2500;   % [s]
t        = (0:1:t_end)';          % 1-second grid  (column vector)

% Build delayed input:  0 W for t<10 s, then 1 W
u        = zeros(size(t));
u(t >= 10) = Pstep;

% Simulate continuous-time plant on the same grid
y = lsim(plant, u, t);            % y(:,1)=ΔT_L, y(:,2)=ΔT_R

% Convert to absolute temperature (°C)
T_L = Tamb + y(:,1);
T_R = Tamb + y(:,2);

%% 5.  Save data 
T = table(t, u, T_L, T_R,'VariableNames', {'t_s', 'Power_W', 'T_local_C', 'T_remote_C'});
writetable(T,'pcb_step_response.csv');
save('pcb_step_response.mat', 'T', 'plant', 'Tamb', 'u');

%% 6.  Plot quick look
figure, plot(t,T_L,'LineWidth',1.3), hold on
plot(t,T_R,'LineWidth',1.3), grid on
xlabel('Time  [s]'), ylabel('Temperature  [°C]')
title('1-W step response'), legend('Local','Remote','Location','southeast')

%% 7.  Analytical transfer functions
G_L = tf(plant(1))     % ΔT_L / P
G_R = tf(plant(2))     % ΔT_R / P

