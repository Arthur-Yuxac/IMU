function [is_ellipsoid, A,B,C] = is_valid_ellipsoid(ellipsoid_params, eps)

if nargin < 2
    eps = 1e-6;
end

a = ellipsoid_params(1); b = ellipsoid_params(2); c = ellipsoid_params(3);
f = ellipsoid_params(4); g = ellipsoid_params(5); h = ellipsoid_params(6);

M = [a, h, g;
     h, b, f;
     g, f, c];

eig_vals = eig(M);
is_ellipsoid = all(eig_vals > eps);


if is_ellipsoid
    [V, Lambda] = eig(M);
    lambda = diag(Lambda);
    A = lambda(1); B = lambda(2); C = lambda(3);
else
    A = 0; B = 0; C = 0;
end
end

