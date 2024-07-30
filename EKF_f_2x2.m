%% EKF

% clear stored variable log_EKF
clear log_EKF;

%% Take data from simulation (simulink)
L  = center_of_mass;
% z_hat dim 2: x, phi
% h dim 2

% output of sensor omega for the f
omega = out.omega.signals.values;

% output of sensors values and times for the h
D = out.D.signals.values; % values of the sensor D
D_times = out.D.time; % times vector of sensor D
D_flags = zeros(size(D_times, 1), 1); % vector of zeros: flags to see if the measure has been used

a1 = out.a1.signals.values; % values of the sensor a1
a1_times = out.a1.time; % times vector of sensor a1
a1_flags = zeros(size(a1_times, 1), 1); % vector of zeros: flags to see if a measure has been used

% Initial values
std_dev_init = 1e-20;
z_hat = [x_0; phi_0] + std_dev_init*randn(2,1);
P = std_dev_init^2 * eye(2);

% Def. Q and R
% Q: noises matrix in f
Q = [std_dev_omega^2 0;
     0 std_dev_phi^2];
% R: sensors' variances in h
R = [std_dev_D^2 0;
     0 std_dev_a1^2];

% Time of simulation
t_max = (out.SimulationMetadata.ModelInfo.StopTime);

% dt = fastest sensor -> omega of the wheels, that is present in the f
dt = T_sensor_omega; % step for
n = int32(t_max / dt); % number of iterations
k = 0;

z_hat_tot = zeros(n,2); % x_hat total, all steps inside, matrix nx2
P_tot = zeros(n,2,2); % P total, 3 dimensions (n x 3 x 3)


%% Loop for prediction and correction

i=0; % counter to check how many corrections

for dt_sum = dt:dt:t_max
    
    k = k + 1;

    % here k is good index for omega beacause we are using dt = its
    % sampling time -> so when we increment k is a dt more in the loop,
    % meaning that there are the right number of values of omega

    % prediction
    [z_hat, P, F] = predict_EKF(z_hat, P, Q, dt, omega(k), wheel_radius);
    
    % save data for smoother
    log_EKF.z_hat_prediction(k,1:2) = z_hat;
    log_EKF.P_prediction(k,:,:) = P;
    log_EKF.F_matrix(k,:,:) = F;

        
    % find target index
    canCorrect = 0; % variable to check if correction can be done
    i_D = find(D_times <= dt_sum, 1, "last"); % 1 last index that verify the condition
    
    if D_flags(i_D) == 0
        canCorrect = 1; % measure of D has not been used
    end

    i_a1 = find(a1_times <= dt_sum, 1, "last");
    
    if ((a1_flags(i_a1) == 0) && (canCorrect == 1)) % check is value of a1 has not been used already
        a1_flags(i_a1) = 1;
        D_flags(i_D) = 1;
        canCorrect = 2; % we can correct
    end

    if canCorrect == 2
        % correction
        i = i + 1; % how many corrections
        % possiamo fare correzione di solo una delle uscite
        [z_hat, P] = correct_EKF(z_hat, P, R, D(i_D), a1(i_a1), L, L_a1, wheel_radius, S_a1);
        % different index for correction, otherwise array much bigger
    end
    % save for smoothing
    log_EKF.z_hat_correction(k,1:2) = z_hat;
    log_EKF.P_correction(k,:,:) = P;


    z_hat_tot(k,1:2) = z_hat;
    log_EKF.z_hat_tot(k,1:2) = z_hat; % for smoothing plot
    P_tot(k,:,:) = P;
    
end



%% PLOTS

% Extract columns
x_ekf = z_hat_tot(:, 1);
phi_ekf = z_hat_tot(:, 2);

time_ekf = linspace(0, t_max, n);
log_EKF.time_ekf = time_ekf;

save('EKF_smoother', 'log_EKF');

% Plot each column in a separate figure
figure();
subplot(2, 1, 1);
plot(time_ekf, x_ekf);
hold on;

x = getElement(out.yout,'x').Values.Data;
x_time = getElement(out.yout,'x').Values.Time; % time vector (timeserie)
plot(x_time, x);
title('x: EKF vs real');

subplot(2, 1, 2);
plot(time_ekf, phi_ekf);
hold on;

phi = getElement(out.yout,'phi').Values.Data;
phi_time = getElement(out.yout,'phi').Values.Time; % time vector (timeserie)
plot(phi_time, phi);
title('phi: EKF vs real');

%% Functions definitions
function [z_hat, P, F] = predict_EKF(z_hat, P, Q, dt, omega, r)

    F = [1 0;
         0 1];
    
    % noise vector w has dim = 2
    D = [r*dt 0;
         0 dt];
    
    z_hat = z_hat + dt * [omega*r; 0];
    P = F*P*F' + D*Q*D';
end

function [z_hat, P] = correct_EKF(z_hat, P, R, D, a1, L, L_a1, r, S_a1)
    
    x_curr = z_hat(1);
    phi_curr = z_hat(2);

    % function of sensor a1:
    % h(2) = atan2(x_curr - S_a1 + L*sin(phi_curr), L_a1 - L*cos(phi_curr) - r)
    
    % dh(2)/dx
    K1 = -(r - L_a1 + L*cos(phi_curr))/((r - L_a1 + L*cos(phi_curr))^2 + (x_curr - S_a1 + L*sin(phi_curr))^2);
    
    % dh(2)/dphi
    K2 = -(L*(L - L_a1*cos(phi_curr) - S_a1*sin(phi_curr) + r*cos(phi_curr) + x_curr*sin(phi_curr)))/...
        ((r - L_a1 + L*cos(phi_curr))^2 + (x_curr - S_a1 + L*sin(phi_curr))^2);
    
    H = [1 0;
         K1 K2];

    M = [1 0;
         0 1];
    
    e = [D; a1] - [x_curr; atan2(x_curr - S_a1 + L*sin(phi_curr), L_a1 - L*cos(phi_curr) - r)];
    e(2) = wrapToPi(e(2));
    
    S = H*P*H' + M*R*M';
    L_mat = P*H'*S^(-1);

    z_hat = z_hat + L_mat*e;
    P = (eye(2) - L_mat*H)*P*(eye(2) - L_mat*H)' + L_mat*M*R*M'*L_mat';
end