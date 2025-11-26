function  Euler = CtoEuler(C) % Euler kısmı [roll, pitch, yaw] idi dikkat et!!
%3-2-1 sequence
yaw = atan2(C(1,2) , C(1,1));

%pitch = - acos(C(3,3));
pitch = -asin(C(1,3));
%roll = atan2(C(3,1) , -C(3,2));
roll = atan2(C(2,3), C(3,3));
Euler = [roll, pitch, yaw];

end

