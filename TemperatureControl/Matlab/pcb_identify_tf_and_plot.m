%======================================================================
% pcb_identify_tf_and_plot.m  —  use logged input vector
% pcb_step_response.csv
%======================================================================

clc
clear

initTemperatureControl
%% 1.  Load data -------------------------------------------------------
csvFile = 'pcb_step_response.csv';         % adjust path if needed
T       = readtable(csvFile);
t       = T.t_s;                           % [s] time vector
u       = T.Power_W;                       % [W] logged power input
y_L     = T.T_local_C;                     % [°C] local sensor
y_R     = T.T_remote_C;                    % [°C] remote sensor

Ts      = median(diff(t));                % sample time (should be 1 s)

%% 2.  Create iddata objects ------------------------------------------
Tamb    = y_L(1);                          % use first point as ambient
y_Ldz   = y_L - Tamb;                      % zero bias
y_Rdz   = y_R - Tamb;

data_L  = iddata(y_Ldz, u, Ts, 'Name','Local');
data_R  = iddata(y_Rdz, u, Ts, 'Name','Remote');

%% 3.  Estimate transfer functions -------------------------------------
np_L = 2;  nz_L = 1;                       % 2 poles, 1 zero for local
np_R = 2;  nz_R = 0;                       % 2 poles, 0 zeros for remote

sys_L = tfest(data_L, np_L, nz_L);        % ΔT_L / P
sys_R = tfest(data_R, np_R, nz_R);        % ΔT_R / P

disp('Identified transfer function – Local sensor (ΔT_L / P):')
sys_L
disp('Identified transfer function – Remote sensor (ΔT_R / P):')
sys_R



%% 3.1  Get state-space matrices of the identified models --------------
% Convert TF → minimal state-space (continuous-time)
ss_L = minreal( ss(sys_L) );      % Local sensor model
ss_R = minreal( ss(sys_R) );      % Remote sensor model


%% 4.  Simulate identified models over the same time vector -----------
yhat_L = lsim(sys_L, u, t);               % simulated local temperature
yhat_R = lsim(sys_R, u, t);               % simulated remote temperature

% Add ambient back in:
yhat_L = yhat_L + Tamb;
yhat_R = yhat_R + Tamb;

%% 5.  Plot comparison -------------------------------------------------
figure('Name','Step response comparison','Color','w')

subplot(2,1,1)
plot(t, y_L, 'k', 'LineWidth',1.2), hold on
plot(t, yhat_L, 'r--', 'LineWidth',1.2)
grid on, xlabel('Time [s]'), ylabel('Temperature [°C]')
title('Local sensor – data vs. identified TF')
legend('Measured','Identified model','Location','southeast')

subplot(2,1,2)
plot(t, y_R, 'k', 'LineWidth',1.2), hold on
plot(t, yhat_R, 'r--', 'LineWidth',1.2)
grid on, xlabel('Time [s]'), ylabel('Temperature [°C]')
title('Remote sensor – data vs. identified TF')
legend('Measured','Identified model','Location','southeast')

%% 6.  Delay estimation  (local → remote) ------------------------------
% Uses the identified transfer functions (sys_L, sys_R) and same time grid t

% --- 6.1  Simulate 1-W ideal step (ΔT only, no ambient) ---------------
u_id  = ones(size(t));              % 1-W step at t = 0
yL_id = lsim(sys_L, u_id, t);       % local model output
yR_id = lsim(sys_R, u_id, t);       % remote model output

% --- 6.2  Steady-state gains ------------------------------------------
K_L = dcgain(sys_L);                % °C/W
K_R = dcgain(sys_R);

% --- 6.3  50 % and 63 % rise delay ------------------------------------
idx50_L = find(yL_id >= 0.50*K_L, 1, 'first');
idx50_R = find(yR_id >= 0.50*K_R, 1, 'first');
idx63_L = find(yL_id >= 0.632*K_L, 1, 'first');
idx63_R = find(yR_id >= 0.632*K_R, 1, 'first');

dT50 = (idx50_R - idx50_L) * Ts;    % seconds
dT63 = (idx63_R - idx63_L) * Ts;

fprintf('\n--- Estimated Local→Remote Delay ---\n');
fprintf('Δt50  (50 %% rise):      %.1f  s\n', dT50);
fprintf('Δt63  (63 %% rise):      %.1f  s\n', dT63);

% --- 6.5 overlay plot ---------------------------------------
figure('Name','Model-based delay illustration','Color','w')
plot(t, yL_id, 'b', 'LineWidth',1.2), hold on
plot(t, yR_id, 'r', 'LineWidth',1.2)
xline(t(idx50_L), ':b', '50% L'); xline(t(idx50_R), ':r', '50% R');
xline(t(idx63_L), '--b', '63% L'); xline(t(idx63_R), '--r', '63% R');
xlabel('Time  [s]'), ylabel('\DeltaT  [°C]')
title(sprintf('Model step responses – Δt50 = %.1f s,  Δt63 = %.1f s',dT50,dT63))
legend('Local (model)','Remote (model)','Location','southeast')
grid on