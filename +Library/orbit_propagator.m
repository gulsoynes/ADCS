function S = orbit_propagator(tspan, sys_0)
%Orbit propogation around Earth for circular orbit and J2 perturbation

%Earth Parameters
% u = 3.986004415e5; % km^3/s^2
% R_earth = 6371; % km
J_num = [0;0.0010826269;-0.0000025323;-0.0000016204]; % zonal harmonic coefficients
J2 = J_num(2);

u = 3.98601e14;         %Earth Gravitational Constant (m^3/s^2)
R_earth = 6.371e6;       % m


options = odeset("RelTol",1e-10, "AbsTol",1e-10*ones(1,6));

[t,S] = ode45(@(t,S) OrbitState(t, S, u, R_earth,J2), tspan, sys_0, options); % propagate

end

function [State] = OrbitState(t, S,u,R_earth, J2)
% find new position magnitude at each time step
r = norm(S(1:3));
%r = sqrt(S(1)^2 + S(2)^2 + S(3)^2);
% verify to perturb the orbit or not
%a_dist = Get_Disturbance(J2,r,u,R_earth,S(1:3));
a_dist = [0, 0, 0];
State = [S(4)
    S(5)
    S(6)
    (-u/r^3)*S(1) + a_dist(1)
    (-u/r^3)*S(2) + a_dist(2)
    (-u/r^3)*S(3) + a_dist(3)];
end

function [a_dist,z] = Get_Disturbance(J2,r,u, R_earth,S)

%From textbook p 387
x = S(1); y = S(2); z = S(3);
a_dist = (-3/2 * J2 * (u/r^2) * (R_earth / r)^2) * ...
    [((1-5*(z/r)^2) * x/r); ((1-5*(z/r)^2) * y/r); ((3-5 * (z/r)^2) * z/r)];
end