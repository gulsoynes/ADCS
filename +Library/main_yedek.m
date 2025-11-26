%2022-2023 GRADUATION PROJECT : SINGLE FRAME and KALMAN FILTERING BASE
%METHODS for ATTITUDE DETERMINATION

%This code contains caculations of Earth's magnetic field, sun and nadir
%direction vectors in body and reference frame.
%Second part contains TRIAD method, q method and Kalman Filtering.

%%%Written by Neslihan Gülsoy

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
w = [.0002; .0003; 0.002];

% Initial Perifokal (PQW) position and velocity vectors
r_pf = [ p.orb.r; 0; 0 ];
v_pf = [ 0; p.orb.v; 0 ];

Q = PQWtoECI(p.orb, t);   % transformation matrix from orbit frame to ECI

r_eci = Q * r_pf;
v_eci = Q * v_pf;

% Initial position and velocity vectors of satellite

%Orbit Propagation
orbit = orbit_propagator(p.tspan, [r_eci; v_eci]);   %Edit considering orbital elements
pos = orbit(:,1:3); % position vector of satellite

% Satellite Dynamics
X = [q;w];

for i = 1:p.N
    % Find quaternion and angular velocity
    X(:,i+1) = dynamics(X(:,i), p.sat.I, p.sat.Moment, p.dt);    
end

for i = 1:p.N+1
    Euler_ang_true(:,i) = qtoEuler(X(1:4, i));

    [C(i).a] = qtoC(X(1:4, i));
    
    % Sensor measurement
    t = t + p.dt;
    [b_o(:,i), b_b(:,i), s_o(:,i), s_b(:,i), n_o(:,i), n_b(:,i)]  = sensor_model(pos(i,:)', p.orb, p.sensor, C(i).a, t);

end

%% Only Nadir and Magnetometer
str = 'Nadir-Magnetometer'
for i = 1:p.N+1
    
    sigma = [std(b_b(:,i)) + std(b_o(:,i));
        std(n_b(:,i)) + std(n_o(:,i))];
    
    sgn = sign(X(1:4,i));
     
    b = [n_b(:,i), b_b(:,i)];
    r = [n_o(:,i), b_o(:,i)];
    
    [A_est, P_Triad(i).a] = TRIAD(r,b,sigma);
    [q_TRIAD(:,i)] = Ctoq(A_est, sgn);
   
    [q_qMethod(:,i), P_qMethod(i).a] = qMethod(r,b,sigma,sgn);
    
end

%% Only Sun and Magnetometer
str = 'Sun-Magnetometer'
for i = 1:N+1
    
    sigma = [std(s_b(:,i)) + std(s_o(:,i));
        std(b_b(:,i)) + std(b_o(:,i))];
       
    sgn = sign(X(1:4,i));
    
    b = [s_b(:,i), b_b(:,i)];
    r = [s_o(:,i), b_o(:,i)];
    
    [q_TRIAD(:,i), P_Triad(i).a ] = TRIAD(r,b,sigma,sgn);
    Euler_ang_TRIAD(:,i) = qtoEuler(q_TRIAD(:,i));
 
    [q_qMethod(:,i), P_qMethod(i).a] = qMethod(r,b,sigma,sgn);
    Euler_ang_qMethod(:,i) = qtoEuler(q_qMethod(:,i));
end

%% Only Sun and Nadir Vectors Obs.
str = 'Sun-Nadir'
for i = 1:N+1
    
    sigma = [std(s_b(:,i)) + std(s_o(:,i));
        std(n_b(:,i)) + std(n_o(:,i))];
    sgn = sign(X(1:4,i));
    
    b = [s_b(:,i), n_b(:,i)];
    r = [s_o(:,i), n_o(:,i)];
    
    [q_TRIAD(:,i), P_Triad(i).a] = TRIAD(r,b,sigma,sgn);
    Euler_ang_TRIAD(:,i) = qtoEuler(q_TRIAD(:,i));
    
    [q_qMethod(:,i), P_qMethod(i).a] = qMethod(r,b,sigma,sgn);
    Euler_ang_qMethod(:,i) = qtoEuler(q_qMethod(:,i));
end

%% q Method All Combinations 
str = 'Sun-Magnetometer-Nadir'
for i = 1:N+1
    
    sigma = [std(s_b(:,i)) + std(s_o(:,i));
        std(n_b(:,i)) + std(n_o(:,i));
        std(b_b(:,i)) + std(b_o(:,i))];
       
    sgn = sign(X(1:4,i));
    
    b = [s_b(:,i), b_b(:,i), n_b(:,i)];
    r = [s_o(:,i), b_o(:,i), n_o(:,i)];
 
    [q_qMethod3(:,i)] = qMethod3(r,b,sigma,sgn);
    Euler_ang_qMethod3(:,i) = qtoEuler(q_qMethod3(:,i));
end

RMSE_q_Method3 = sqrt( sum( abs( rad2deg(Euler_ang_true)' - rad2deg(Euler_ang_qMethod3)' ).^2 ) / (length(Euler_ang_true)) )
%% Single Frame Aided EKF
tic
[X,P_EKF] = EKF_plus(q_TRIAD,Euler_ang_TRIAD,X(5:7,:),P_Triad,p.N);
[X2,P_q] = EKF_plus(q_qMethod,Euler_ang_qMethod,X(5:7,:),P_qMethod,pN);
n = length(X);
for i = 1:n
    Euler_ang_EKF(:,i) = qtoEuler(X(1:4,i));
    Euler_ang_EKF2(:,i) = qtoEuler(X2(1:4,i));
end
% Root Mean Square Error (RMSE)
toc
%%
RMSE_TRIADq = sqrt( sum( abs((q_true)' - (q_TRIAD)' ).^2 ) / n )

RMSE_q_Methodq = sqrt( sum( abs( (q_true)' - (q_qMethod)' ).^2 ) / n )

RMSE_TRIAD_EKFq = sqrt( sum( abs( (q_true)' - (X(1:4,:))' ).^2 ) / n )

RMSE_qMethod_EKF2q = sqrt( sum( abs( (q_true)' - (X2(1:4,:))' ).^2 ) / n )

%%
RMSE_TRIAD = sqrt( sum( abs( rad2deg(Euler_ang_true)' - rad2deg(Euler_ang_TRIAD)' ).^2 ) / n )

RMSE_q_Method = sqrt( sum( abs( rad2deg(Euler_ang_true)' - rad2deg(Euler_ang_qMethod)' ).^2 ) / n )

% RMSE_TRIAD_EKF = sqrt( sum( abs( rad2deg(Euler_ang_true)' - rad2deg(Euler_ang_EKF)' ).^2 ) / n )
% 
% RMSE_qMethod_EKF2 = sqrt( sum( abs( rad2deg(Euler_ang_true)' - rad2deg(Euler_ang_EKF2)' ).^2 ) / n )


 %%
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
  params.sat.Moment = 5e-10;      %The disturbance torque acting on the satellite (N.m)
  params.sensor.gyro_noise = 1e-5;
  params.filter.Q = 1e-6 * eye(6);
  params.filter.R = 1e-4 * eye(3);
  % ... diğer parametreler

% Sensor Parameters

%Parameters for Sensors
params.sensor.magneto_sd = .008;     %The standard deviation of magnetometer error
params.sensor.sun_sd = .002;     %The standard deviation of each Sun sensor noise
params.sensor.nadir_sd = .006;     %The standard deviation of each Horizon sensor noise


%Orbital Elements
h = 500e3;                      %Altitude of Satellite (m)
Re = 6378.140e3;                %Earth radius on Equator (m)

mu = 3.98601e14;        %Earth Gravitational Constant (m^3/s^2)

% Orbital elements,
params.orb.r = ( Re + h );                % Semi-major axis (m)
params.orb.inc = 81; %97.65;                       %Otbit Inclination (deg)
params.orb.RAAN = 171.8428;                %RAAN (deg)
params.orb.argp = 0;                      %Argument of perigee (deg)
params.orb.ecc = 0; % Eccentricity 
params.orb.arlat_0 = 0;    % Initial argument of latitude (deg)
params.orb.wo = sqrt( mu / ( params.orb.r ^ 3 ) );   %Angular Velocity of Orbit (rad/s)
params.orb.v = sqrt(mu / params.orb.r);
params.orb.T = 2 * pi * sqrt(params.orb.r^3 / mu);
%Iteration Parameters   
params.N = params.orb.T * 1 * 5;
params.tspan = linspace (0 , params.N*params.dt , params.N+1);
%time = linspace (0 , N*delt , N+1);
end