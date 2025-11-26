function [X_upd,P_upd, bias_upd] = kf_update(X_pred, P_pred, z, R)
% Inputs 
% z_quad : Measured quaternion from TRIAD or q-Method

% Constants
H = [eye(3) zeros(3);
    zeros(3) zeros(3)];     % Measurement matrix

% Kalman Gain Calculation
K = P_pred * transpose(H) / (H * P_pred * transpose(H) + R);

% Update Covariance matrix
P_upd = (eye(6) - K * H) * P_pred;

q_pred = X_pred(1:4,:);
q_meas = z(1:4,:);
w_meas = z(5:7,:);
% Update error-state

% Error state update formula from  Optimal Estimation of Dynamic Systems page 456
Ksi = Ksi_fc(q_pred);
% burası düzeltilecek. 
delta_x = K * [(2 * transpose(Ksi) * q_meas); w_meas];

q_upd = q_pred + 1/2 * Ksi * delta_x(1:3,:);
bias_pred = zeros(3,1);
w_upd = X_pred(5:7,:);
bias_upd = bias_pred + delta_x(4:6,:);

X_upd = [q_upd(1:4,:); w_upd];

end

function Ksi = Ksi_fc(q)

Ksi = [q(4) -q(3) q(2);
    q(3) q(4) -q(1);
    -q(2) q(1) q(4);
    -q(1) -q(2) -q(3)];
end