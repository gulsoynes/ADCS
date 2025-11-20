% Residuals

load("HIP.mat");
% s -- star tracker vector in inertial frame

star_ID = [32349, 30438, 71683, 69673]; %Sirius, Canopus, Rigil Kent, Arcturus

RA_deg = HIP(star_ID,2);
DE_deg = HIP(star_ID,3);
V_Mag = HIP(:,4);