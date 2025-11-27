%2022-2023 GRADUATION PROJECT : SINGLE FRAME and KALMAN FILTERING BASE
%METHODS for ATTITUDE DETERMINATION

%This code contains caculations of Earth's magnetic field, sun and nadir
%direction vectors in body and reference frame.
%Second part contains TRIAD method, q method and Kalman Filtering.

%%%Written by Neslihan Gülsoy

% To find a Earth-pointing satellite
% LVLH frame is the main reference frame

%%%%%%%%%% NOTE : Gyro bias update kısmına bak 

% Attitude is body wrt. orbit frame.

%Update : The qtoEuler and CtoEuler give same results
% 
clear; clc; close all;
tic

%Period of satellite is almost 1.5 hours, so iterarion number is 2*T
%Parameters for simulation
p = params();

%Initial State
t = 0;
q = [0;0;0;1];
w = [0.0001; 0.0002; 0.0001];

bias = deg2rad(0.1)/(60*60);  % initial bias is 1 deg/hr
bias = 0;
bias = [bias; bias; bias];

% Initial Orbit (LVLH) position and velocity vectors
r_or = [ p.orb.r; 0; 0 ];
v_or = [ 0; -p.orb.v; 0 ];


Tr_RTN2ECI = Library.RTNtoECI(p.orb, t);

r_eci = Tr_RTN2ECI * r_or;
v_eci = Tr_RTN2ECI * v_or;

% Initial position and velocity vectors of satellite

%Orbit Propagation
orbit = Library.orbit_propagator(p.tspan, [r_eci; v_eci]);   %Edit considering orbital elements
pos = orbit(:,1:3)'; % position vector of satellite
vel = orbit(:,4:6)'; % velocity vector of satellite

% Satellite Dynamics
X = [q;w];

% Iterate for n times orbit time = N : 2*T_orb
for i = 1:p.N
    % Find quaternion and angular velocity
    X(:,i+1) = Library.dynamics(X(:,i), p.sat.I, p.sat.Moment, p.orb, p.dt);    
end

for i = 1:p.N+1
    Euler_ang_true(:,i) = Library.qtoEuler(X(1:4, i));

    % Calculate attitude matrix
    C_BwrtO(i).a = Library.qtoC(X(1:4,i));
    % Sensor measurement
    t = p.tspan(i);
    [w_gyro(:,i), bias(:,i+1)] = Library.gyro_measurement(X(5:7,i), bias(:,i), p.sensor, p.dt);

    [b_o(:,i), b_b(:,i), s_o(:,i), s_b(:,i), n_o(:,i), n_b(:,i)]  =...
        Library.sensor_model(pos(:,i), vel(:,i), p.orb, p.sensor, C_BwrtO(i).a, t);
end

% %%
% figure(1)
% plot(p.tspan, rad2deg(Euler_ang_true(1,:)))
% hold on
% plot(p.tspan, rad2deg(Euler_ang_true(2,:)))
% hold on
% plot(p.tspan, rad2deg(Euler_ang_true(3,:)))
% legend("Roll", "Pitch", "Yaw")
% xticks([0, p.orb.T/2, p.orb.T, 1.5*p.orb.T, 2*T])
%xticklabels({'0','T/2','T','3/2T','2T','2\pi','3\pi'})
% Only Sun and Magnetometer
str = 'Sun-Magnetometer'
for i = 1:p.N+1
    
    sigma = [std(s_b(:,i)) + std(s_o(:,i));
        std(b_b(:,i)) + std(b_o(:,i))];
       
    sgn = sign(X(1:4,i));
    
    b = [s_b(:,i), b_b(:,i)];
    r = [s_o(:,i), b_o(:,i)];
    
    [A_est, P_Triad(i).a] = Library.TRIAD(r,b,sigma);
    [q_TRIAD(:,i)] = Library.Ctoq(A_est, sgn);
    Euler_TRIAD(:,i) = Library.qtoEuler(q_TRIAD(:,i));


    [q_qMethod(:,i), P_qMethod(i).a] = Library.qMethod(r,b,sigma,sgn);
end

for i=1:p.N+1
     Error_q(:,i) = Error_quad(q_TRIAD(:,i), X(1:4,i));
end
% Measurement
z = [q_TRIAD; w_gyro];
%z = [X(1:4,:); w_gyro];

figure(1)
plot(Euler_TRIAD(1,:))
hold on
plot(Euler_ang_true(1,:))
%% Single Frame Aided EKF
tic

% Initial states 
% Start algorithm with some noise on angular velocity and measured
% quaternion from TRIAD
X_est(:,1) = [q_TRIAD(:,1); w_gyro(:,1)];
bias_est = zeros(3,1);

P = [1*eye(3), zeros(3);
    zeros(3), 9.5e-10*eye(3)];  %Initial Cov Matrix

std_Euler_error = [0.005; 0.008; 0.0038];
%std_Euler_error = deg2rad(std_Euler_error.^2);

% Measurement Covariance Matrix with TRIAD 
R = diag(std_Euler_error);

Q = [100*p.sensor.gyro_noise^2*eye(3), zeros(3);
    zeros(3), 100*p.sensor.gyro_bias_sd^2*eye(3)];

for i = 1:1:p.N
    [X_est(:,i), P_upd, bias_est] = Library.MEKF_Update(X_est(:,i), P, z(:,i), bias_est, R);
    [X_est(:,i+1), P] = Library.MEKF_Propagate(X_est(:,i), P_upd, p.dt, Q);

    b_est(:,i) = bias_est;
end

toc

%%

for i=1:p.N+1
Euler_EKF(:,i) = Library.qtoEuler(X_est(:,i));
end


figure(1)
plot(Euler_TRIAD(2,:))
hold on
plot(Euler_ang_true(2,:))
hold on
plot(Euler_EKF(2,:))

legend("TRIAD", "true", "EKF")
%%
w_angular = X(5:7,:);

figure(1)
plot(w_angular(2,:))
hold on
plot(w_gyro(2,:))
legend("true","gyro")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 PLOTTING(n,time,str,rad2deg(Euler_ang_EKF),...
    rad2deg(Euler_ang_TRIAD),...
    rad2deg(Euler_ang_qMethod),...
    rad2deg(Euler_ang_true))
%%
figure(1)
plot(  time(1:n),rad2deg(Euler_ang_EKF(1,1:n)),'.--')
hold on
plot(  time(1:n), rad2deg(Euler_ang_TRIAD(1,1:n)),'.--')
hold on
plot(  time(1:n), rad2deg(Euler_ang_qMethod(1,1:n)),'.--')
hold on
plot(  time(1:n), rad2deg(Euler_ang_true(1,1:n)))
ax=gca;
ax.YGrid = 'on';
ax.XLim = [0   time(n)];
xlabel('$ t $ : Time $(s)$','Interpreter', 'latex ')
ylabel('$ \phi$ : Roll Angle $(^{\circ})$','Interpreter','latex')
lg = legend(' $ EKF $ ',  ' $TRIAD $ ' ,'$ q-Method $ ','True','Interpreter','latex');
title(lg, str)
title(str,'Interpreter', 'latex')

% Parameter definition for code
function params = params()
  params.dt = 0.1;              %Sample time (s)
  params.sat.I = diag([2.1e-3, 2e-3, 1.9e-3]); % atalet matrisi (m^4)
  params.sat.mass = 50;          % kg
  params.sat.Moment = 0;      %The disturbance torque acting on the satellite (N.m)

  params.filter.Q = 1e-6 * eye(6);
  params.filter.R = 1e-4 * eye(3);
  % ... diğer parametreler

% Sensor Parameters

%Parameters for Sensors
params.sensor.gyro_noise = sqrt(10) * 10^(-10);      % Standard deviation gyro noise
params.sensor.gyro_bias_sd = sqrt(10) * 10^(-13);   % Standard deviation bias noise    
params.sensor.magneto_sd = .008;     %The standard deviation of magnetometer error
params.sensor.sun_sd = .002;     %The standard deviation of each Sun sensor noise
params.sensor.nadir_sd = .006;     %The standard deviation of each Horizon sensor noise


%Orbital Elements
h = 500e3;                      %Altitude of Satellite (m)
Re = 6378.140e3;                %Earth radius on Equator (m)

mu = 3.98601e14;        %Earth Gravitational Constant (m^3/s^2)

% Orbital elements,
params.orb.r = ( Re + h );                % Semi-major axis (m)
params.orb.inc = 97.65;                       %Otbit Inclination (deg)
params.orb.RAAN = 171.8428;                %RAAN (deg)
params.orb.argp = 0;                      %Argument of perigee (deg)
params.orb.ecc = 0; % Eccentricity 
params.orb.arlat_0 = 0;    % Initial argument of latitude (deg)
params.orb.wo = sqrt( mu / ( params.orb.r ^ 3 ) );   %Angular Velocity of Orbit (rad/s)
params.orb.v = sqrt(mu / params.orb.r);
params.orb.T = 2 * pi * sqrt(params.orb.r^3 / mu);
%Iteration Parameters   
params.n = 2;   % 2 times orbit time
params.N = params.orb.T * 1/params.dt * params.n;
params.tspan = linspace (0 , params.orb.T*params.n , params.N+1);
end