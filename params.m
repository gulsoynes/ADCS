clear; clc;

p = params();

function params = params()
  % params: proje genel ayarları

  params.dt = 0.1;              %Sample time (s)
  params.sat.I = diag([2.1e-3, 2e-3, 1.9e-3]); % atalet matrisi (m^4)
  params.sat.mass = 50;          % kg
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
h = 400e3;                      %Altitude of Satellite (m)
Re = 6378.140e3;                %Earth radius on Equator (m)
Ro = ( Re + h );                %Distance Between Satellite and Earth (m)
params.sat.N_t = 5e-10;      %The disturbance torque acting on the satellite (N.m)

params.sat.inc = 89; %97.65;                       %Otbit Inclination (deg)
params.sat.Omega = 207.4;                     %Right Ascension of the Ascending Node (deg)
params.sat.omega = 10;                      %Argument of perigee (deg)

%Iteration Parameters
delt = .1;          
params.N = 1.5 * (60 * 60) * 2;
%time = linspace (0 , N*delt , N+1);

% w_angular = [Sat.w1_0; Sat.w2_0; Sat.w3_0];




end