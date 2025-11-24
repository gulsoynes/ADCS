function Tr = PQWtoECI(orb, t)
% funciton gives transformation matrix from orbit frame to ECI
% Tr depends on RAAN, inc, argp and true anomaly
wo = orb.wo;
% Convert angles to radians
inc = deg2rad(orb.inc);         % Geodetic inclination
RAAN = deg2rad(orb.RAAN);   % right ascending node
argp = deg2rad(orb.argp);   % argument of perigee (zero for circular orbit)
nu = deg2rad(orb.arlat_0) + wo * t; % true anomaly
nu = mod(nu, 2*pi);
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

end

