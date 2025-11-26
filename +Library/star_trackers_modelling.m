function [s_i, s_b] = star_trackers_modelling(RA_deg,DE_deg,C)
sigma_sensor1 = [5; 5; 40] * 0.0002777777777777778;  %arcsec to radians
sigma_sensor = [5; 5; 40] * 10^-4;

for i = 1:4
    s_i(1:3, i) = [cosd(RA_deg(i)) * cosd(DE_deg(i));
        sind(RA_deg(i)) * cosd(DE_deg(i));
        sind(DE_deg(i))];
end
s_i = s_i ./ vecnorm(s_i);

for i = 1:4
    s_b(1:3, i) = C * s_i(1:3,i) + sigma_sensor .* randn(3,1);
end
end