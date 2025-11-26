clear; clc;

load('data.mat')

r = r_ECI(:,10);
B = B_ECI(:,10)
% v_eci = [-3384 -4887 4843];
% 
% [r_ECI, r_ECEF] = ECEF2ECI_pos(lat,long, alt, utc);


B_e = DipoleModel(r)


function [r_ECI, r_ECEF] = ECEF2ECI_pos(lat,long, alt, utc)
% Constants
% ECEF to ECI transformation with given latitude, longitude, altitude and
% time
R_a = 6378137;      % Semimajor axis in WGS-84 model
R_p = 6356752.3142;
f = (R_a - R_p) / R_a;
e = sqrt(1-(1-f)^2);
omega_E = 7.292115e-5; % Earth's rotation rate (rad/s)

N = R_a / (sqrt(1 - e^2 * sind(lat)^2));

r_ECEF = [(N + alt) * cosd(lat) * cosd(long);
    (N+alt) * cosd(lat) * sind(long);
    (N * (1-e^2) + alt) * sind(lat)];

% Transform ECEF to ECI using UTC
year = utc(1);  % year
month = utc(2); % month
day = utc(3);   % day
hour = utc(4);
minute = utc(5);
second = utc(6);

JD = JDate(year,month,day,hour,minute,second);

T0 = (JD - 2451545) / (36525);

theta_GMS = 24110.54841 + 8640184.812866 * T0 + 0.093104 * T0^2 ...
    - 6.2e-6 * T0^3 + 1.002737909350795 * (3600 * hour + 60 * minute + second);  % Eq.(2.70)

theta_GMS = mod(theta_GMS, 86400);      % Eq.(2.70)
theta_GMS = theta_GMS / 240;            % Eq.(2.70) (seconds to degree)

R_ECEF2ECI = [cosd(theta_GMS) -sind(theta_GMS) 0 ;
    sind(theta_GMS) cosd(theta_GMS) 0;
    0 0 1];

r_ECI = R_ECEF2ECI * r_ECEF;
v_ECI = cross([0; 0; omega_E], r_ECI);

end