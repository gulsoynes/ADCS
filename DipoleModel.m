function B = DipoleModel(r_eci)
% Reduced IGRF model is used to define Earth's magnetic field vector
% lon : array
% longitude [deg], positive east, of IGRF calculation
% lat : array
% geodetic latitude [deg] of IGRF calculation
% alt : array
% height [km] above ellipsoid for IGRF calculation
% date : date(s)
% 
% Return
% ------
% Be : array
% Magnetic field [nT] in eastward direction
% Bn : array
% Magnetic field [nT] in northward direction, relative to
% ellipsoid
% Bu : array
% Magnetic field [nT] in upward direction, relative to
% ellipsoid

% IGRF 2005 parameters
theta_m = 169.7;        % coef for magnetic dipole calc. (deg)
alpha_m = 108;        % coef for magnetic dipole calc. (deg)
m_s = 7.77e22;          % Magnetic dipole magnitude (A.m^2)

m = m_s * [sind(theta_m)*cosd(alpha_m);...
    sind(theta_m) * sind(alpha_m);...
    cosd(theta_m)];

% V = a^3/norm(r) * (-29554.63 * sind(lat) + ...
%     (-1669.05 * cosd(lat) * sind(long)) + ...
%     507799 * cosd(lat) * sind(long));

utc = [2020 5 1 0 0 0];
r = eci2ecef(utc, r_eci);

B_ECEF = 3 * (dot(m,r) * r - norm(r)^2 * m) / (norm(r)^5);
% Its in NED frame. Need go to ECI. Need adding lat, long etc
B = ECEF2ECI(utc, B_ECEF);

end


function r_ECI = ECEF2ECI_pos(lat,long, alt, utc)
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