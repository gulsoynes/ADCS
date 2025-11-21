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
q = [1;0;0;0];
w = [.0002; .0003; 0.002];
% Satellite Dynamics
%time = linspace (0 , N*delt , N+1);
X = [q;w];

time = 
X = dynamics(X, p.sat.I, p.sat.Moment, p.dt);

%%
%Sensor Modelling
%[H_o, H_b] = Magnetometer(p.Ro,p.sat.inc,,sigma_h,time,N);

%%
[S_o, S_b, S_ECI] = SunSensor(omega,Omega,inc,T,C,sigma_s,time,N);
[N_o, N_b] = HorizonSensor(Re,Ro,omega,Omega,T,C,sigma_n,time,N);

for i = 1:length(C)
[s_i(i).a, s_b(i).a] = star_trackers_modelling(RA_deg,DE_deg,C(i).a);
end


%% Only Nadir and Magnetometer
str = 'Nadir-Magnetometer'
for i = 1:N+1
    
    sigma = [std(H_b(:,i)) + std(H_o(:,i));
        std(N_b(:,i)) + std(N_o(:,i))];
    
    sgn = sign(q_true(:,i));
     
    b = [N_b(:,i), H_b(:,i)];
    r = [N_o(:,i), H_o(:,i)];
    
    [q_TRIAD(:,i), P_Triad(i).a ] = TRIAD(r,b,sigma,sgn);
    Euler_ang_TRIAD(:,i) = qtoEuler(q_TRIAD(:,i));
   
    [q_qMethod(:,i), P_qMethod(i).a] = qMethod(r,b,sigma,sgn);
    Euler_ang_qMethod(:,i) = qtoEuler(q_qMethod(:,i));
    
end

%% Only Sun and Magnetometer
str = 'Sun-Magnetometer'
for i = 1:N+1
    
    sigma = [std(S_b(:,i)) + std(S_o(:,i));
        std(H_b(:,i)) + std(H_o(:,i))];
       
    sgn = sign(q_true(:,i));
    
    b = [S_b(:,i), H_b(:,i)];
    r = [S_o(:,i), H_o(:,i)];
    
    [q_TRIAD(:,i), P_Triad(i).a ] = TRIAD(r,b,sigma,sgn);
    Euler_ang_TRIAD(:,i) = qtoEuler(q_TRIAD(:,i));
 
    [q_qMethod(:,i), P_qMethod(i).a] = qMethod(r,b,sigma,sgn);
    Euler_ang_qMethod(:,i) = qtoEuler(q_qMethod(:,i));
end

%% Only Sun and Nadir Vectors Obs.
str = 'Sun-Nadir'
for i = 1:N+1
    
    sigma = [std(S_b(:,i)) + std(S_o(:,i));
        std(N_b(:,i)) + std(N_o(:,i))];
    sgn = sign(q_true(:,i));
    
    b = [S_b(:,i), N_b(:,i)];
    r = [S_o(:,i), N_o(:,i)];
    
    [q_TRIAD(:,i), P_Triad(i).a, A_TRIAD(i).a] = TRIAD(r,b,sigma,sgn);
    Euler_ang_TRIAD(:,i) = qtoEuler(q_TRIAD(:,i));
    
    [q_qMethod(:,i), P_qMethod(i).a] = qMethod(r,b,sigma,sgn);
    Euler_ang_qMethod(:,i) = qtoEuler(q_qMethod(:,i));
end

%% q Method All Combinations 
str = 'Sun-Magnetometer-Nadir'
for i = 1:N+1
    
    sigma = [std(S_b(:,i)) + std(S_o(:,i));
        std(N_b(:,i)) + std(N_o(:,i));
        std(H_b(:,i)) + std(H_o(:,i))];
       
    sgn = sign(q_true(:,i));
    
    b = [S_b(:,i), H_b(:,i), N_b(:,i)];
    r = [S_o(:,i), H_o(:,i), N_o(:,i)];
 
    [q_qMethod3(:,i)] = qMethod3(r,b,sigma,sgn);
    Euler_ang_qMethod3(:,i) = qtoEuler(q_qMethod3(:,i));
end

RMSE_q_Method3 = sqrt( sum( abs( rad2deg(Euler_ang_true)' - rad2deg(Euler_ang_qMethod3)' ).^2 ) / (length(Euler_ang_true)) )
%% Star Trackers
str = 'Star_Trackers'

for i = 1:N+1 
    sigma = [std(s_b(i).a) + std(s_i(i).a)]';
       
    sgn = sign(q_true(:,i));
    
    b = [s_b(i).a];
    r = [s_i(i).a];
 
    [q_qMethod_star(:,i)] = qMethod3(r,b,sigma,sgn);
    Euler_ang_qMethod3(:,i) = qtoEuler(q_qMethod_star(:,i));

end

RMSE_q_Method3 = sqrt( sum( abs( rad2deg(Euler_ang_true)' - rad2deg(Euler_ang_qMethod3)' ).^2 ) / (length(Euler_ang_true)) )
%% Single Frame Aided EKF
tic
[X,P_EKF] = EKF_plus(q_TRIAD,Euler_ang_TRIAD,w_angular,P_Triad,N);
[X2,P_q] = EKF_plus(q_qMethod,Euler_ang_qMethod,w_angular,P_qMethod,N);
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
h = 600e3;                      %Altitude of Satellite (m)
Re = 6378.140e3;                %Earth radius on Equator (m)
params.Ro = ( Re + h );                %Distance Between Satellite and Earth (m)

mu = 3.98601e14;        %Earth Gravitational Constant (m^3/s^2)
params.Wo = sqrt( mu / ( Ro ^ 3 ) );   %Angular Velocity of Orbit (rad/s)

params.sat.inc = 89; %97.65;                       %Otbit Inclination (deg)
params.sat.RAAN = 171.8428;                %RAAN (deg)
params.sat.argp = 0;                      %Argument of perigee (deg)
params.sat.ecc = 0; % Eccentricity 
params.sat.arlat_0 = 10;    % Initial argument of latitude (deg)

%Iteration Parameters   
params.N = 1.5 * (60 * 60) * 2;
%time = linspace (0 , N*delt , N+1);
end