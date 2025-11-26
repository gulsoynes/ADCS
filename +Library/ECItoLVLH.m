function Tr_ECI2Orb = ECItoLVLH(pos, vel)
% funciton gives transformation matrix from orbit frame (LVLH) to ECI
% LVLH frame
% z axis : nadir direction etc center of Earth
% y axis : the opposite of normal vector of orbital plane (cross-track)
% x axis : complete the right handed orthogonal frame (along-track)

% Inputs :
% pos : position vector of satellite
% vel : velocity vector of satellite

% Output :
% Transformation matrix from LVLH to ECI frame
% Note that transpose of this matrix gives transformation from ECI to LVLH

% Equations : z axis is directing Earth
z = - pos ./ norm(pos);

y = - cross(pos,vel) ./ norm(cross(pos,vel));

x = cross(y, z);

Tr_ECI2Orb = [x y z];

end