function [H_o, H_b, B_ECI] = Magnetometer(Ro,inc,C,sigma_h,time,N)
% Tilted Dipole Model is used to generate magnetic field vector.
% Note that this model operates on Orbit (PQW) Frame.
% Equations are taken from DOI:10.1109/TMECH.2013.2259843
% Input 
Me = 7.943e24;          %Magnetic Dipole Moment of Earth (Wb.Gm)
We = 7.29e-5;           %Spin Rate of Earth (rad/s)
e = deg2rad (11.7);     %Magnetic Dipole Tilt (rad)
mu = 3.98601e14;        %Earth Gravitational Constant (m^3/s^2)
Wo = sqrt( mu / ( Ro ^ 3 ) );   %Angular Velocity of Orbit (rad/s)
%b_c = [.04; .06; .08];      %magnetometer bias vector
for i=1:N+1
    %Magnetic Field Vector in Orbit Frame
    %X Component of Magnetic Field Vector (nT)
    Hx = ( Me / Ro^3 ) * ( cos( Wo*time(i) ) * ( cos(e) * sind(inc) ...
        - sin(e) * cosd(inc) * cos ( We*time(i) ) ) ...
        - sin(e) * sin(Wo*time(i)));

    %Y Component of Magnetic Field Vector (T)
    Hy = ( -Me / Ro^3 ) * ( cos(e) * cosd(inc) ...
        + sin(e) * sind(inc) * cos ( We*time(i) ) );

    %Z Component of Magnetic Field Vector (T)
    Hz = ( 2 * Me / Ro^3 ) * ( sin (Wo * time(i) ) * (cos(e) * sind(inc) ...
        -sin(e)*cosd(inc)*cos(We*time(i))) ...
        + sin(e) * cos(Wo*time(i)));

    H = [Hx; Hy; Hz]; %Magnetic Field Vector

    %Direction Cosine Components of Ho
    Hxo(i,1) = Hx; %/norm(Hx)
    Hyo(i,1) = Hy;
    Hzo(i,1) = Hz;

    %Direction Cosine of Magnetic Field Vector in Orbit Reference
    H_o(:,i) = [Hxo(i,1) ; Hyo(i,1); Hzo(i,1)];

    %Measured Direction Cosine of Magnetic Field Vector in Body Frame
    H_b0 = C(i).a * H_o(:,i) + sigma_h * randn(3,1);
    H_b(:,i) = H_b0;%/norm

    Hxb(i,1) = H_b(1,i);
    Hyb(i,1) = H_b(2,i);
    Hzb(i,1) = H_b(3,i);

end

    Omega = 171.8428;                %RAAN GRACE1 (deg)
    omega = 349.0575;                      %Argument of perigee (deg)
    u = omega;
    RAAN = Omega;
    AoP = omega;    % Argument of perigee
    %Transformation Matrix ECI to Orbit Frame
    Tr = [-cosd(u)*cosd(inc)*sind(Omega)-sind(u)*cosd(Omega),...
        cosd(u)*cosd(inc)*cosd(Omega)-sind(u)*sind(Omega),...
        cosd(u)*sind(inc);
        ...
        -sind(inc)*sind(Omega),...
        sind(inc)*cosd(Omega),...
        -cosd(inc);
        ...
        sind(u)*cosd(inc)*sind(Omega)-cosd(u)*cosd(Omega),...
        -sind(u)*cosd(inc)*cosd(Omega)-cosd(u)*sind(Omega),...
        -sind(u)*sind(inc)];

B_ECI = Tr' * H_o;
end