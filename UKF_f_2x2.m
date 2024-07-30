%% Take data from simulation (simulink)

% z_hat dim 2: x, phi
% h dim 2

% out sensori per la f
omega = out.omega.signals.values;

% output of sensor omega for the f
D = out.D.signals.values; % values of the sensor D
D_times = out.D.time; % times vector of sensor D
D_flags = zeros(size(D_times, 1), 1); % vector of zeros: flags to see if the measure has been used

a1 = out.a1.signals.values; % values of the sensor a1
a1_times = out.a1.time; % times vector of sensor a1
a1_flags = zeros(size(a1_times, 1), 1); % vector of zeros: flags to see if the measure has been used

% note: the value of sensors are in struct format

% Initial values
std_dev_init = 1e-5;       %---> std_dev_init < 0.004
z_hat = [x_0; phi_0] + std_dev_init*randn(2,1);
P = std_dev_init^2 * eye(2);

% Def. Q and R
% Q: matrice varianza rumore di processo nella f
Q = [std_dev_omega^2 0;    %<------------
     0 std_dev_phi^2];
% R: matrice varianze dei sensori nella h
R = [std_dev_D^2 0;
     0 std_dev_a1^2];

% Time of simulation
t_max = (out.SimulationMetadata.ModelInfo.StopTime);

% dt = fastest sensor -> omega of the wheels, that is present in the f
dt = T_sensor_omega; % step for
n = int32(t_max / dt); % number of iterations
k = 0;

z_hat_tot = zeros(n,2); % x_hat totale di tutti i passi, matrice nx2
P_tot = zeros(n,2,2); % P total, 3 dimensions (n x 2 x 2)


%% Inizio ciclo per predizione e correzione
i = 0;

for dt_sum = dt:dt:t_max
    
    k = k + 1;


    % prediction
    [z_hat, P] = predict_UKF(z_hat, P, Q, dt, omega(k), wheel_radius);
        
    % find target index
    canCorrect = 0; % variable to check if correction can be done
    i_D = find(D_times <= dt_sum, 1, "last"); % 1 last index that verify the condition

    if D_flags(i_D) == 0
        canCorrect = 1; % measure of D not been used
    end

    i_a1 = find(a1_times <= dt_sum, 1, "last");
    if ((a1_flags(i_a1) == 0) && (canCorrect == 1)) % check is value of D has not been used already
        a1_flags(i_a1) = 1;
        D_flags(i_D) = 1;
        canCorrect = 2; % we can correct
    end

    if canCorrect == 2
        % correction
        i = i+1; % how many corrections
        % possiamo fare correzione di solo una delle uscite
        [z_hat, P] = correct_UKF(z_hat, P, R, D(i_D), a1(i_a1), L, L_a1, wheel_radius, S_a1);
    end

    z_hat_tot(k,1:2) = z_hat; % riga k, colonne tutte e due
    P_tot(k,:,:) = P;
    
end


%% PLOTS

% Extract columns
x_ukf = z_hat_tot(:, 1);
phi_ukf = z_hat_tot(:, 2);

time_ukf = linspace(0, t_max, n);

% Plot each column in a separate figure
figure;
subplot(2, 1, 1);
plot(time_ukf, x_ukf);
hold on;

x = getElement(out.yout,'x').Values.Data;
x_time = getElement(out.yout,'x').Values.Time; % time vector (timeserie)
plot(x_time, x);
title('x: UKF vs real');

subplot(2, 1, 2);
plot(time_ukf, phi_ukf);
hold on;

phi = getElement(out.yout,'phi').Values.Data;
phi_time = getElement(out.yout,'phi').Values.Time; % time vector (timeserie)
plot(phi_time, phi);
title('phi: UKF vs real');

%% Functions definitions

% predict function
function [z_hat, P] = predict_UKF(z_hat, P, Q, dt, omega, r)
    [z_hat, P] = UnscentedTransform_F_Function([z_hat; 0; 0], blkdiag(P, Q), dt, omega, r);

end

% correct function
function [z_hat, P] = correct_UKF(z_hat, P, R, D, a1, L, L_a1, r, S_a1)
    [hat_y, S, Pxy] = UnscentedTransform_H_Function(z_hat, P, L, L_a1, r, S_a1);
    y = [D; a1];
    S = S + R;
    e = y - hat_y;
    e = atan2(sin(e), cos(e));
    L = Pxy*S^(-1);
    
    z_hat = z_hat + L*e;
    P = P - L*S*L';
end



% UT f function
function [propagated_mean, propagated_cov] = UnscentedTransform_F_Function(prior_mean, prior_cov, dt, omega, r)
    % define paramenters
    alpha = 1;            % not scaled
    beta = 2;             % optimum for the Gaussian case
    kappa = 0;            % lambda = 0 [weight of the central sigma_point]
    
    enne = size(prior_mean,1); % prior_mean=[x; phi; 0; 0], enne=4
    
    lambda = alpha^2*(enne + kappa) - enne;
    
    % compute weights
    w0 = lambda/(lambda + enne);
    wc(1) = w0 + 1 - alpha^2 + beta;
    wc(2:2*enne+1) = 1/2/(enne+lambda);
    
    wm = wc;
    wm(1) = w0;
    
    % factorise covariance matrix
    % SVD
    [U,S,~] = svd(prior_cov);
    GAMMA = U*S^(1/2); %GAMMA*GAMMA'=prior_cov
    
    % generate 2n+1 sigma points
    sigma_points = prior_mean;
    for i = 1:size(GAMMA,2)
        sigma_points(:,i+1)      = prior_mean + sqrt(enne+lambda)*GAMMA(:,i);
        sigma_points(:,i+1+enne) = prior_mean - sqrt(enne+lambda)*GAMMA(:,i);
        % sigma_point=matrice 4x(2n+1), ogni colonna è un sigma point
    end
    
    % apply function f to the 2n+1 sigma points
    for i = 1:size(sigma_points,2)
        propagated_sigma_points(1:2,i) = sigma_points(1:2,i) + ...
            dt*[(omega + sigma_points(3,i))*r; sigma_points(4,i)];
        % propagated_sigma_point=matrice 3x(2n+1)
    end
    
    % compute propagated mean
    propagated_mean = propagated_sigma_points*wm';
    
    % propagated ~
    propagated_tilde = propagated_sigma_points - propagated_mean;

    % compute virtual measurements second order
    propagated_cov = zeros();%zeros(4,4)
    for i = 1:size(sigma_points,2)
        propagated_cov = propagated_cov + wc(i)*propagated_tilde(:,i)*propagated_tilde(:,i)';
    end
end




% UT h function
function [virtual_measurements_mean, virtual_measurements_cov, cross_cov] = UnscentedTransform_H_Function(prior_mean, prior_cov, L, L_a1, r, S_a1)

    % define paramenters
    alpha = 1;            % not scaled
    beta = 2;             % optimum for the Gaussian case
    kappa = 0;            % lambda = 0 [weight of the central sigma_point]
    
    enne = size(prior_mean,1);
    
    lambda = alpha^2*(enne + kappa) - enne;
    
    % compute weights
    w0 = lambda/(lambda + enne);
    wc(1) = w0 + 1 - alpha^2 + beta;
    wc(2:2*enne+1) = 1/2/(enne+lambda);
    
    wm = wc;
    wm(1) = w0;
    
    % factorise covariance matrix
    % SVD
    [U,S,~] = svd(prior_cov);
    GAMMA = U*S^(1/2);
    
    % generate 2n+1 sigma points
    sigma_points = prior_mean;
    for i = 1:size(GAMMA,2)
        sigma_points(:,i+1)      = prior_mean + sqrt(enne+lambda)*GAMMA(:,i);
        sigma_points(:,i+1+enne) = prior_mean - sqrt(enne+lambda)*GAMMA(:,i);
    end
    
    % apply function h to the 2n+1 sigma points
    for i = 1:size(sigma_points,2)
         propagated_sigma_points(1,i) = sigma_points(1,i);
         propagated_sigma_points(2,i) = atan2(L_a1 - (L*cos(sigma_points(2,i)) + r), S_a1 - (sigma_points(1,i) + L*sin(sigma_points(2,i)))) + (3/2)*pi;
      
    end
  
    % compute virtual measurements mean (NOTE: mean of angles)
    virtual_measurements_mean(1,1) = propagated_sigma_points(1,:)*wm';
    virtual_measurements_mean(2,1) = atan2(sin(propagated_sigma_points(2,:))*wm', cos(propagated_sigma_points(2,:))*wm');

    
    % virtual measurements error ~y (NOTE: difference of angles)
    virtual_measurements_tilde = atan2(sin(propagated_sigma_points - virtual_measurements_mean),cos(propagated_sigma_points - virtual_measurements_mean));
    
    % compute virtual measurements second order
    virtual_measurements_cov = zeros(); %zeros(2,2)
    cross_cov = zeros(); %zeros(2,2)
    for i = 1:size(sigma_points,2)
        virtual_measurements_cov = virtual_measurements_cov + wc(i)*virtual_measurements_tilde(:,i)*virtual_measurements_tilde(:,i)';
        cross_cov = cross_cov + wc(i)*(sigma_points(:,i) - prior_mean)*virtual_measurements_tilde(:,i)';
    end
end