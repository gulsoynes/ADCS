function [X, C] = dynamics(X, J, M, dt)
% This function handles spacecraft dynamics using equation of motion for given paramaters. Equations
% from refence book. And simple integration method is used, can be
% improved. Don't forget the check parameters
%Inputs : Parameters for satellite, and simulation and states
% X : states [q, w]
% J : diagonal nertia matrix
% M : moments actiıng on satellite
% dt : time step
%Outputs : quaternion and angular velocity

J_x = J(1,1); J_y = J(2,2); J_z = J(3,3);
w_x = X(5); w_y = X(6); w_z = X(7);
q1 = X(1); q2 = X(2); q3 = X(3); q0 = X(4); %Scalar part of quaternion

    %Quaternion Rates (1/s)
    
    q1_dot =  .5 * ( q0 * w_x - q3 * w_y + q2 * w_z );
    q2_dot =  .5 * ( q3 * w_x + q0 * w_y - q1 * w_z );
    q3_dot = -.5 * ( q2 * w_x - q1 * w_y - q0 * w_z );
    q0_dot = -.5 * ( q1 * w_x + q2 * w_y + q3 * w_z );
    %Quaternions
    
    q1 = q1 + dt * q1_dot;
    q2 = q2 + dt * q2_dot;
    q3 = q3 + dt * q3_dot;
    q0 = q0 + dt * q0_dot;
    
    %The Angular Velocities (rad/s)
    
    w_x(:,1) = w_x + dt / (J_x) * (J_y - J_z) * w_z * w_y + ...
        dt / (J_x) * M;
    w_y(:,1) = w_y + dt / (J_y) * (J_z - J_x) * w_x * w_z + ...
        dt / (J_y) * M;
    w_z(:,1) = w_z + dt / (J_z) * (J_x - J_y) * w_x * w_y + ...
        dt / (J_z) * M;

    % If J is not diagonal; use below equations 
    % wDot = (J)\(M - cross(w,J*w));
    % w = w + wDot * dt;
    
q = [q1; q2; q3; q0];
q = q ./ norm(q);
w_angular = [w_x; w_y; w_z];
X = [q; w_angular];

% Tranformation matrix from ECI to Body 
C = [(q0^4 + q1^2 - q2^2 - q3^2), 2*(q1*q2 + q3*q0), 2*(q1*q3 - q2*q0);
    2*(q1*q2 - q3*q0), q0^2-q1^2+q2^2-q3^2, (2*(q2*q3 +q1*q0));
    2*(q1*q3 + q2*q0), 2*(q2*q3 - q1*q0), (q0^2 - q1^2 -q2^2 + q3^2)];
end