function B = Dipole_Model(r,i)

mu_B = 7.812e6;     %Earth's magnetic constant
omega_0 = sqrt(mu/(r^3));

u_0 = 0;
u = omega_0 * t + u_0;

%Magnetic Field of direct dipole in Orbit Frame :

B = mu_B/r^3 * [cos(u) * sin(i);
                -cos(i);
                2 * sin(u) * sin(i)];

end

