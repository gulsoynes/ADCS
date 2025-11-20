function [q,w] = EoM(q, w, I, mass, dt)
% This function handles spacecraft dynamics using equation of motion for given paramaters. Equations
% from refence book. And simple integration method is used, can be
% improved. Don't forget the check parameters
%Inputs : Parameters for satellite, and simulation and states
%Outputs : quaternion and angular velocity

w_dot = I\(mass - cross(w, I*w));

w = w + w_dot * dt;

a = q; b = w;
a1 = a(1); a2 = a(2); a3 = a(3);
bv = b(1:3); %[b1;b2;b3];
av = a(1:3); %[a1;a2;a3];

%Check length of vectors and assign scalar part
if length(b) == 4
    b0 = b(4);

elseif length(b) ~= 4
    b0 = 0;
end

if length(a) == 4
    a0 = a(4);

elseif length(a) ~= 4
    a0 = 0;
end

% Mathematical Operations Eq.(2.86) from reference book
crs = [0 -a3 a2;
    a3 0 -a1;
    -a2 a1 0];

out = [a0*eye(3)+crs, av ;
    -transpose(av), a0] * [bv;b0];

q_dot = 1/2 * out; %Eq.(3.20)

%Quaternions as state
q = q + q_dot * dt;
q = q/norm(q);           % Normalizing quaternion
end

