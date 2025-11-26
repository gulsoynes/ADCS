function [b_o, b_b, s_o, s_b, n_o, n_b]  = sensor_model(pos, vel, orb, sensor, C, t)
% ınput : position vector of satellite,
% orbit parameters: radius, inclination and angular velocity
%
% output : magnetic field unit vector, sun vector
% constants
Me = 7.943e24;              % Earth’s magnetic dipole strength (Wb.Gm)
w_E = 7.29e-5;              % Spin Rate of Earth (rad/s)
gamma = deg2rad (11.44);    % Magnetic Dipole Tilt (rad)

% orbit parameters
inc = deg2rad(orb.inc);     % Geodetic inclination
wo = orb.wo;
r = orb.r;          % orbit radius (m)

% Transformation matrix from ECI to LVLH
% Sen sensor gives measurement in ECI frame.
Tr_ECI2Orb = ECItoLVLH(pos, vel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% MAGNETOMETER MODEL %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth's Magnetic Field using Tilted Dipole Model
% Assuming that inc < 78.56 degree; 
% inc_magnetic = inc_geodetic - tilt_of_MF
% Equations from DOI:10.1109/TMECH.2013.2259843
B_y = - Me / r^3 * (cos(inc) * cos(gamma) + sin(inc) * sin(gamma) * cos(w_E * t));

B_x = Me / r^3 * (cos(wo * t) * (sin(inc) * cos(gamma) - ...
    cos(inc) * sin(gamma) * cos(w_E * t)) - sin(wo *t) * sin(gamma));

B_z = 2 * Me / r^3 * (sin(wo * t) * (sin(inc) * cos(gamma) - cos(inc) * sin(gamma) * cos(w_E * t)) ...
    + cos(wo*t) * sin(gamma));

% B_y = - Me / r^3 * (cos(inc + gamma));
% B_x = Me / r^3 * (cos(wo * t) * sin(inc + gamma));
% B_z = 2 * Me / r^3 * (sin(wo * t) * sin(inc + gamma));

b_o = [B_x ; B_y; B_z];   % magnetic field vector in orbit reference frame
b_o = b_o ./norm(b_o);

% Magnetometer model
b_b = C * b_o + sensor.magneto_sd * randn(3,1);
b_b = b_b ./ norm(b_b);    % unit magnetic vector in body frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SUN SENSOR MODEL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference : doi.org/10.1007/978-1-4939-0802-8 pp:420-422
% Sun Heading Vector in ECI frame
JD = JDate(2022,1,1,0,0,0);
T_TDB = (JD - 2451545.0) / 36525;

%mean longitude of the sun (deg)
Lamdba_Msun = 280.460 + 36000.770 * T_TDB;

%mean anomaly of the sun (deg)
M_sun = 357.5277233 + 35999.05034 * T_TDB;

%ecliptic longitude of the sun (deg)
Lamdba_ec = Lamdba_Msun + 1.914666471 * sind(M_sun) ...
    + 0.019994643 * sind(2 * M_sun);

%linear model of the ecliptic of the sun (deg)
e_linear = 23.439291 - 0.0130042 * T_TDB;

% heading from earth to sun
e_E2sun = [cosd(Lamdba_ec);...
    sind(Lamdba_ec) * cosd(e_linear);...
    sind(Lamdba_ec) * sind(e_linear)];

% Distance between Earth and Sun in AU
r_E2sun = 1.000140612 - 0.016708617 * cos(M_sun) - 0.000139589 * cos(2*M_sun);

% Position vector of Earth to Sun is
rv_E2sun = r_E2sun * e_E2sun;

pos_AU = pos * 6.68458712e-12; % satellite position vector in AU

% Position vector of satellite to Sun in AU
pos2sun = rv_E2sun - pos_AU;

% Eclipse might be added.
% Sun position vector wrt ECI frame
s_ECI = pos2sun / norm(pos2sun);

% Sun position vector wrt Orbit frame
s_o = Tr_ECI2Orb * s_ECI;
s_o = s_o ./ norm(s_o);

% Sun sensor model
s_b = C * s_o + sensor.sun_sd * randn(3,1);
s_b = s_b ./ norm(s_b);  % unit sun sensor model in body frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% NADİR SENSOR MODEL %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_ECI = - pos / norm(pos);

n_o = Tr_ECI2Orb * n_ECI;
n_o = n_o./norm(n_o);

% Nadir sensor model
n_b = C * n_o + sensor.nadir_sd * randn(3,1);
n_b = n_b ./ norm(n_b);

end
%%%%%%%%%%%%%%%%%%%%% Additional Functions %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Julian Date Calculator %%%%%%%%%%%%%%%%%%%
function [jd] = JDate(y,m,d,h,min,s)
% Function to compute Julian Date. Matlab built-in function can be also
% used.
%Julian Date Converter
% Input : year, month, day, hour, min and second
% Output: juliandate
jd = 367 * y - fix(7 * (y + fix ((m + 9) / 12) ) /4) + fix(275 * m / 9)...
    + (h + min/60 + s/3600)/24 + d + 1721013.5;
end

%%%%%%%%%%%%%%%%%% ECI to LVLH Transformation %%%%%%%%%%%%%%%%%%%%
function Tr_ECI2Orb = ECItoLVLH(pos, vel)
% funciton gives transformation matrix from orbit frame (LVLH) to ECI
% LVLH frame
% z axis : nadir direction etc center of Earth
% y axis : the opposite of normal vector of orbital plane (cross-track)
% x axis : complete the right handed orthogonal frame (along-track)

% Inputs :
% pos : position vector of satellite
% vel : velocity vector of satellite

% Output :
% Transformation matrix from LVLH to ECI frame
% Note that transpose of this matrix gives transformation from ECI to LVLH

% Equations : z axis is directing Earth
z = - pos ./ norm(pos);

y = - cross(pos,vel) ./ norm(cross(pos,vel));

x = cross(y, z);

Tr_ECI2Orb = [x y z]';

end