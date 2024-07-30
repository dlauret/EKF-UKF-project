%% Smoother with RTS algorithm

load('EKF_smoother.mat')

% Smoother Initialisation
% z_n|n & P_n|n equal to the last correction
log_EKF.z_hat_smoothed(size(log_EKF.z_hat_prediction,1),:) = log_EKF.z_hat_correction(size(log_EKF.z_hat_correction,1),:);
log_EKF.P_smoothed(size(log_EKF.z_hat_prediction,1),:,:) = [log_EKF.P_correction(size(log_EKF.z_hat_correction,1),:,1);
log_EKF.P_correction(size(log_EKF.z_hat_correction,1),:,2)];


% Smoother execution
% recursive in reverse order


for k = size(log_EKF.z_hat_prediction,1)-1:-1:1 % k index of PREDICTION
    
    % Ck = P_k|k * F_k+1' * P_k+1|k^(-1)
    Ck = [log_EKF.P_correction(k,:,1);log_EKF.P_correction(k,:,2)] ...
        * [log_EKF.F_matrix(k+1,:,1); log_EKF.F_matrix(k+1,:,2)]' / ...
        [log_EKF.P_prediction(k+1,:,1);log_EKF.P_prediction(k+1,:,2)];
    
    % z_k|n = z_k|k + Ck*(z_k+1|n - z_k+1|k)
	log_EKF.z_hat_smoothed(k,:) = log_EKF.z_hat_correction(k,:)' + Ck*(log_EKF.z_hat_smoothed(k+1,:)' - log_EKF.z_hat_prediction(k+1,:)');
    
    
    % P_k|n = P_k|k + Ck*(P_k+1|n - P_k+1|k)*Ck'
    log_EKF.P_smoothed(k,:,:) = [log_EKF.P_correction(k,:,1);log_EKF.P_correction(k,:,2)] ...
        + Ck*([log_EKF.P_smoothed(k+1,:,1);log_EKF.P_smoothed(k+1,:,2)] - [log_EKF.P_prediction(k+1,:,1);log_EKF.P_prediction(k+1,:,2)])*Ck';
end



%% Plots

x_ekf = log_EKF.z_hat_tot(:, 1);
phi_ekf = log_EKF.z_hat_tot(:, 2);
time_ekf = log_EKF.time_ekf;


% Plot each column in a separate figure
figure();
% x
subplot(2, 1, 1);
% EKF
plot(time_ekf, x_ekf, 'g');
hold on;
% Smoothed
plot(time_ekf, log_EKF.z_hat_smoothed(:,1), 'b');
hold on;
% Real system
x = getElement(out.yout,'x').Values.Data;
x_time = getElement(out.yout,'x').Values.Time; % time vector (timeserie)
plot(x_time, x, 'r');
title('x: Smoother vs EKF vs real');

% phi
subplot(2, 1, 2);
% EKF
plot(time_ekf, phi_ekf, 'g');
hold on;
% Smoothed
plot(time_ekf, log_EKF.z_hat_smoothed(:,2), 'b');
hold on;
% Real system
phi = getElement(out.yout,'phi').Values.Data;
phi_time = getElement(out.yout,'phi').Values.Time; % time vector (timeserie)
plot(phi_time, phi, 'r');
title('phi: Smoother vs EKF vs real');

