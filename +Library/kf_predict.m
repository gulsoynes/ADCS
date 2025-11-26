function [X_pred, P_pred] = kf_predict(X_est, w_gyro, P, bias, sensor, dt)
% Inputs 
% X = [quaternion; angular velocities]
% P : covariance matrix
% p : parameters

q_pred = X_est(1:4,:);
% Constants

% Assuming there is no gyro bias but only noise. If there is bias, update G
% and Q matrices with bias 

G = [-eye(3), zeros(3);
    zeros(3), eye(3)]; 

Q = [sensor.gyro_noise^2*eye(3), zeros(3);
    zeros(3), sensor.bias_sd^2 - eye(3)];

% Linearized F matris calculation
omega = X_est(5:7,:);   % Angular velocity of satellite

F = [-cross_product(omega), -eye(3);
    zeros(3), zeros(3)];

% Dynamics propagation
w_pred = w_gyro - bias;
q_pred = quad_propagation(q_pred, w_pred, p.dt);

X_pred = [q_pred; w_pred];
% Covariance propagation
P_dot = F * P + P * transpose(F) + G * Q * transpose(G);

P_pred = P + P_dot * dt;

end

function crs = cross_product(v)

crs = [0 -v(3) v(2);
    v(3) 0 -v(1);
    -v(2) v(1) 0];
end

function q = quad_propagation(q, w, dt)
% This function handles spacecraft dynamics using equation of motion for given paramaters. Equations
% from refence book. And simple integration method is used, can be
% improved. Don't forget the check parameters
%Inputs : Parameters for satellite, and simulation and states
% dt : time step
%Outputs : quaternion and angular velocity

q1 = q(1); q2 = q(2); q3 = q(3); q0 = q(4); %Scalar part of quaternion
w_x = w(1); w_y = w(2); w_z = w(3);
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
   
q = [q1; q2; q3; q0];
q = q ./ norm(q);

end