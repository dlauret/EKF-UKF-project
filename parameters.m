%% In this section the parameters of the equation on motion (EoM) are specified

% Parameters
wheel_radius = 0.2;             % [m]
wheel_mass = 2;                 % [kg]
body_mass = 0.5;               % [kg]
semiaxis_wheels = 0.2;          % [m]
gravity = 9.81;                 % [m/s^2]
center_of_mass = 0.4;          % [m]

% Pack parameters into a structure for better organization
eom_params = struct('wheel_radius', wheel_radius, ...
                'body_mass', body_mass, ...
                'wheel_mass', wheel_mass, ...
                'semiaxis_wheels', semiaxis_wheels, ...
                'gravity', gravity, ...
                'center_of_mass', center_of_mass);

% Initial conditions
x_0 = 0;  % Initial position
dx_0 = 0;  % Initial velocity
phi_0 = pi/2;  % Initial angular position
dphi_0 = 0;  % Initial angular velocity

% Initial conditions vector
V0 = [x_0; dx_0; phi_0; dphi_0];

%% Input parameters
C_input = 0.01; % torque [N*m]
freq_input = 1; % [Hz]

%% Parameters of the 3 sensors

% Parameters of the sensor that measures D
% Hokuyo UTM-30LX-EW LiDAR
% Frequency = 40Hz
std_dev_D = 1e-10;%0.015; % std deviation [m]
T_sensor_D = 0.025; % Sampling time [s]

% Parameters of the sensor that measures omega
% EPC MA3 Miniature Absolute Magnetic Shaft Encoder
% Frequency ≈ 683,733Hz
std_dev_omega = 1e-10;%0.001; % std deviation [rad/s]
T_sensor_omega = 0.001; % Sampling time [s]

% Parameters of the sensor that measures a1
% Camera Model: Intel RealSense D435
% with frame rate of 30 FPS, then T = 1/30
L_a1 = 1.6; % height of a1 [m]
S_a1 = 0.5; % distance from x = 0 [m]
std_dev_a1 = 1e-30;%0.07; % std deviation [rad], noises of different types
T_sensor_a1 = 0.033; % Sampling time [s]

% dev std for the noise of the 2nd equation in the f() of the EKF:
% Noise of phi
std_dev_phi = 1e-4;


