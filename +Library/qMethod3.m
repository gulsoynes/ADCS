function [q_opt] = qMethod3(r, b, sigma, sgn)
b1 = b(:,1);
b2 = b(:,2);
b3 = b(:,3);
sigma1 = sigma(1,:); %More accuarate measurement standart deviation
sigma2 = sigma(2,:); %Less accuarate
sigma3 = sigma(3,:);

[~, N] = size(r);       %Number of Sensors
var = 1./ (sigma).^2;   %Weight of Sensors Measurement
B = zeros(3,3);         %Attitude Profile Matrix

for i = 1:N
    
    B = B + ((var(i,1) * b(:,i)) * r(:,i)');
    
end

z = [B(2,3) - B(3,2);
    B(3,1) - B(1,3);
    B(1,2) - B(2,1)];

%Symmetric Traceless Matrix
K = [B + B' - (trace(B) * eye(3)),  z ;
    z',  trace(B)];

%Find Normalized Eigenvectors and Eigenvalues of K Matrix
[evec, eval] = eig(K);

%Find Eigenvector of Maximum Eigenvalue of K Matrix
[~, index] = max(real(diag(eval)));

%Set Optimum Quaternion as Eigenvector of Maximum Eigenvalue of K Matrix
q_opt = evec(:,index);

%Sign check
index = find(sign(q_opt) ~= sgn);

if ~isempty(index)
    q_opt(index) = - q_opt(index);
end

%Cov = 1/4 * ((sigma1^2 * eye(3)) + ...
   % ((((sigma2^2 - sigma1^2) * (b1 * b1')) + ...
   % (sigma1^2 * (dot(b1, b2) * ((b1* b2') + (b2 * b1'))))) / ...
    %((norm(cross(b1, b2)))^2)));
end