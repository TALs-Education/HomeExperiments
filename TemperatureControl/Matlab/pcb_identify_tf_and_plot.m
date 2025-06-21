%======================================================================
% pcb_identify_tf_and_plot.m
%======================================================================

clc
clear

%% 1.  Load data -------------------------------------------------------
csvFile = 'pcb_step_response.csv';          % adjust path if needed
T   = readtable(csvFile);
t   = T.t_s;                                % [s] time vector
y_L = T.T_local_C;                          % [°C] local sensor
y_R = T.T_remote_C;                         % [°C] remote sensor

%% 2.  Build the input vector (1-W step) ------------------------------
Pstep = 1;                                  % [W] magnitude used in sim
u      = ones(size(t))*Pstep;
u(1)   = 0;                                 % step occurs at t = 0+

Ts = median(diff(t));                       % sample time (1 s here)

%% 3.  Create iddata objects ------------------------------------------
Tamb  = y_L(1);                 % first sample ≈ ambient
y_Ldz = y_L - Tamb;             % de-biased outputs
y_Rdz = y_R - Tamb;

data_L = iddata(y_Ldz, u, Ts, 'Name','Local');
data_R = iddata(y_Rdz, u, Ts, 'Name','Remote');



%% 4.  Estimate transfer functions ------------------------------------
% Orders chosen from physics: 2 poles, 1 zero for local; 2 poles, 0 zeros for remote
np_L = 2;  nz_L = 1;                        % poles, zeros (Local)
np_R = 2;  nz_R = 0;                        % poles, zeros (Remote)

sys_L = tfest(data_L, np_L, nz_L);          % Local TF  ΔT_L(s)/P(s)
sys_R = tfest(data_R, np_R, nz_R);          % Remote TF ΔT_R(s)/P(s)

disp('Identified transfer function – Local sensor (ΔT_L / P):')
sys_L
disp('Identified transfer function – Remote sensor (ΔT_R / P):')
sys_R

%% 5.  Simulate identified models over the same time vector -----------
yhat_L = lsim(sys_L, u, t);            % ΔT_L_hat  (detrended)
yhat_R = lsim(sys_R, u, t);            % ΔT_R_hat

% Add back the ambient offset that was removed by detrend:
Tamb = mean(y_L(1));                        % ambient ≈ first sample
yhat_L = yhat_L + Tamb;
yhat_R = yhat_R + Tamb;

%% 6.  Plot comparison -------------------------------------------------
figure('Name','Step response comparison','Color','w')

subplot(2,1,1)
plot(t, y_L,'k', 'LineWidth',1.2), hold on
plot(t, yhat_L,'r--', 'LineWidth',1.2)
grid on, xlabel('Time  [s]'), ylabel('Temperature  [°C]')
title('Local sensor – data vs. identified TF')
legend('Measured','Identified model','Location','southeast')

subplot(2,1,2)
plot(t, y_R,'k', 'LineWidth',1.2), hold on
plot(t, yhat_R,'r--', 'LineWidth',1.2)
grid on, xlabel('Time  [s]'), ylabel('Temperature  [°C]')
title('Remote sensor – data vs. identified TF')
legend('Measured','Identified model','Location','southeast')

