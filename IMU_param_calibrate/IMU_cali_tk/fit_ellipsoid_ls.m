function [ellipsoid_params, v] = fit_ellipsoid_ls(stable_means, k)


M = size(stable_means, 1);
x = stable_means(:,1); y = stable_means(:,2); z = stable_means(:,3);  % 单位：g

D = [x.^2, y.^2, z.^2, ...
     2*y.*z, 2*x.*z, 2*x.*y, ...
     2*x, 2*y, 2*z, ...
     ones(M,1)];


C1 = [-1,        k/2-1,    k/2-1, 0, 0, 0;
      k/2-1,    -1,        k/2-1, 0, 0, 0;
      k/2-1,    k/2-1,    -1,    0, 0, 0;
      0,        0,        0,    -k, 0, 0;
      0,        0,        0,     0, -k, 0;
      0,        0,        0,     0, 0, -k];
C = blkdiag(C1, zeros(4,4)); 


DDT = D' * D;
S11 = DDT(1:6, 1:6);  
S12 = DDT(1:6, 7:10); 
S22 = DDT(7:10, 7:10);


if rank(S22) < 4
    S22_inv = pinv(S22);
else
    S22_inv = inv(S22);
end


S_schur = S11 - S12 * S22_inv * S12';  
C1_inv = inv(C1);
[V1, Lambda] = eig(C1_inv * S_schur);


eigenvalues = diag(Lambda);
positive_idx = find(eigenvalues > 1e-6);
if length(positive_idx) ~= 1
    [~, max_idx] = max(eigenvalues);
    positive_idx = max_idx;
end
v1 = V1(:, positive_idx);
v2 = -S22_inv * S12' * v1;

v = [v1; v2];
constraint_val = v' * C * v;
v = v / sqrt(constraint_val);

ellipsoid_params = v';
end
