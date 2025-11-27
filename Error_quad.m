function Error_quad = Error_quad(q,qp)
% Error between 2 quaternion q and p
% Quaternion multip
% Inputs : 
% q : Measured
% qp : True
Ksi = Ksi_fc(q);
Error_quad = 2 * transpose(Ksi) * qp;
end

function Ksi = Ksi_fc(q)

Ksi = [q(4) -q(3) q(2);
    q(3) q(4) -q(1);
    -q(2) q(1) q(4);
    -q(1) -q(2) -q(3)];
end