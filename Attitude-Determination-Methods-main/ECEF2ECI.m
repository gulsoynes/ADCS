function v = ECEF2ECI(utc, v)
% Constants
% ECEF to ECI transformation with given latitude, longitude, altitude and
% time

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

v = R_ECEF2ECI * v;

end