function [H_o, H_b] = Sensor_Models(sensor_sd)

%sensor_sd : Standard deviation of sensors, magneto, sun, nadir sensor
Ro = 6378.140e3 + h;        %Distance Between Satellite and Earth (m)
Me = 7.943e15;              %Magnetic Dipole Moment of Earth (Wb.m)
We = 7.29e-5;               %Spin Rate of Earth (rad/s)
e = 0;%deg2rad (11.7);      %Magnetic Dipole Tilt (rad)
mu = 3.98601e14;            %Earth Gravitational Constant (m^3/s^2)
Wo = sqrt( mu / ( Ro ^ 3 ) );   %Angular Velocity of Orbit (rad/s)

    %Magnetic Field Vector in Orbit Frame
    %X Component of Magnetic Field Vector (T)
    Hx = ( Me / Ro^3 ) * ( cos( Wo*time(i) ) * ( cos(e) * sind(inc) ...
        - sin(e) * cosd(inc) * cos ( We*time(i) ) ) ...
        - sin ( Wo*time(i)) * sin(e) * sin ( We*time(i)) );
    
    %Y Component of Magnetic Field Vector (T)
    Hy = ( -Me / Ro^3 ) * ( cos(e) * cosd(inc) ...
        + sin(e) * sind(inc) * cos ( We*time(i) ) );
    
    %Z Component of Magnetic Field Vector (T)
    Hz = ( 2 * Me / Ro^3 ) * ( sin (Wo * time(i) ) * (cos(e) * sind(inc) ...
        -sin(e)*cosd(inc)*cos(We*time(i))) ...
        - 2*sin(Wo*time(i))*sin(e)*sin(We*time(i)));
    
    H = [Hx; Hy; Hz]; %Magnetic Field Vector
    
    %Direction Cosine Components of Ho
    Hxo(i,1) = Hx / norm(H) ;
    Hyo(i,1) = Hy / norm(H);
    Hzo(i,1) = Hz / norm(H);
    
    %Direction Cosine of Magnetic Field Vector in Orbit Reference
    H_o(:,i) = [Hxo(i,1) ; Hyo(i,1); Hzo(i,1)];
    
    %Measured Direction Cosine of Magnetic Field Vector in Body Frame
    H_b0 = C(i).a * H_o(:,i) + sigma_h * randn(3,1);
    H_b(:,i) = H_b0 ./ norm(H_b0);
    
    Hxb(i,1) = H_b(1,i);
    Hyb(i,1) = H_b(2,i);
    Hzb(i,1) = H_b(3,i);

end








%Function for Julian Date calc.
function [jd] = JDate(y,m,d,h,min,s)
%Julian Date Converter

jd = 367 * y - fix(7 * (y + fix ((m + 9) / 12) ) /4) + fix(275 * m / 9)...
    + (h + min/60 + s/3600)/24 + d + 1721013.5; 
end