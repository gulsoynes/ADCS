function q = Eul_to_q(eul)

%It use different sequence from qtoEuler and qtoC

yaw = eul(3); pitch = eul(2); roll = eul(1);

q0 = cos(roll/2) * cos(pitch/2) * cos(yaw/2) - ...
    sin(roll/2) * cos(pitch/2) * sin(yaw/2);

q1 = cos(roll/2) * sin(pitch/2) * cos(yaw/2) + ...
    sin(roll/2) * sin(pitch/2) * sin(yaw/2);

q2 = cos(roll/2) * sin(pitch/2) * sin(yaw/2) - ...
    sin(roll/2) * sin(pitch/2) * cos(yaw/2);

q3 = sin(roll/2) * cos(pitch/2) * cos(yaw/2) + ...
    cos(roll/2) * cos(pitch/2) * sin(yaw/2);

q = [q1; q2; q3; q0];

end

