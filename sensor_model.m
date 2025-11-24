function [b_ECI, b_b, s_ECI, s_b, n_ECI, n_b]  = sensor_model(pos, orb, sensor, C, t)
% ınput : position vector of satellite, 
% orbit parameters: radius, inclination and angular velocity
% 
% output : magnetic field unit vector, sun vector
% constants
Me = 7.943e24;              % Earth’s magnetic dipole strength (Wb.Gm)
w_E = 7.29e-5;              % Spin Rate of Earth (rad/s)
gamma = deg2rad (11.44);    % Magnetic Dipole Tilt (rad)

% orbit parameters
inc = deg2rad(orb.inc);         % Geodetic inclination
wo = orb.wo;
r = orb.r;          % orbit radius (m)
% Convert to radians
RAAN = deg2rad(orb.RAAN);   % right ascending node
argp = deg2rad(orb.argp);   % argument of perigee (zero for circular orbit)

nu = deg2rad(orb.arlat_0) + wo * t; % true anomaly
nu = mod(nu, 2*pi);

%%%%%%%%%%%%%%% MAGNETOMETER MODEL %%%%%%%%%%%%
% Earth's Magnetic Field using Tilted Dipole Model
% Assuming that inc < 78.56 degree; mag inc = inc - gamma
% Equations from DOI:10.1109/TMECH.2013.2259843
B_y = - Me / r^3 * (cos(inc) * cos(gamma) + sin(inc) * sin(gamma) * cos(w_E * t));

B_x = Me / r^3 * (cos(wo * t) * (sin(inc) * cos(gamma) - ...
    cos(inc) * sin(gamma) * cos(w_E * t)) - sin(wo *t) * sin(gamma));

B_z = 2 * Me / r^3 * (sin(wo * t) * (sin(inc) * cos(gamma) - cos(inc) * sin(gamma) * cos(w_E * t)) ...
    + cos(wo*t) * sin(gamma));

% B_y = - Me / r^3 * (cos(inc + gamma));
% B_x = Me / r^3 * (cos(wo * t) * sin(inc + gamma));
% B_z = 2 * Me / r^3 * (sin(wo * t) * sin(inc + gamma));

B_o = [B_x ; B_y; B_z];   % magnetic field vector in orbit reference frame    

% Orbit frame to ECI transformation
    
    R3_RAAN = [ cos(RAAN)  -sin(RAAN)  0;
                sin(RAAN)  cos(RAAN)  0;
                 0          0          1];

    R1_i = [ 1      0           0;
             0   cos(inc)  -sin(inc);
             0  sin(inc)  cos(inc)];

    R3_w_nu = [ cos(argp+nu)  -sin(argp+nu)   0;
               sin(argp+nu)  cos(argp+nu)   0;
                 0             0             1];

    % PQW → ECI
    Tr = R3_RAAN * R1_i * R3_w_nu;

    b_ECI = Tr * B_o;   
    b_ECI = b_ECI/norm(b_ECI); % unit magnetic vector in ECI

    % Magnetometer model
    b_b = C * b_ECI + sensor.magneto_sd * randn(3,1);
    b_b = b_b/norm(b_b);    % unit magnetic vector in body frame
%%%%%%%%%%%%%%%%%%% SUN SENSOR MODEL %%%%%%%%%%%%%%5
    % Sun Heading Vector in ECI frame
    JD = JDate(2022,1,1,0,0,0); 
    T_TDB = (JD - 2451545.0) / 36525;
    
    Lamdba_Msun = 280.460 + 36000.770 * T_TDB; %mean longitude of the sun (deg)
    M_sun = 357.5277233 + 35999.05034 * T_TDB; %mean anomaly of the sun (deg)
    Lamdba_ec = Lamdba_Msun + 1.914666471 * sind(M_sun) ...
        + 0.019994643 * sind(2 * M_sun); %ecliptic longitude of the sun (deg)
    
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

    pos = pos * 6.68458712e-12; % satellite position vector in AU
    % Position vector of satellite to Sun in AU
    pos2sun = rv_E2sun - pos;

    % Eclipse might be added.
    s_ECI = pos2sun / norm(pos2sun);
    % Sun sensor model
    s_b = C * s_ECI + sensor.sun_sd * randn(3,1);
    s_b = s_b / norm(s_b);  % unit sun sensor model in body frame

    %%%%%%%%%%%%%%%%% NADİR VECTOR %%%%%%%%%%%%%%
    n_ECI = - pos / norm(pos);
    % Nadir sensor model
    n_b = C * n_ECI + sensor.nadir_sd * randn(3,1);

end

% Function to compute Julian Date. Matlab built-in function can be also
% used. 
function [jd] = JDate(y,m,d,h,min,s)
%Julian Date Converter
% Input : year, month, day, hour, min and second
% Output: juliandate
jd = 367 * y - fix(7 * (y + fix ((m + 9) / 12) ) /4) + fix(275 * m / 9)...
    + (h + min/60 + s/3600)/24 + d + 1721013.5; 
end