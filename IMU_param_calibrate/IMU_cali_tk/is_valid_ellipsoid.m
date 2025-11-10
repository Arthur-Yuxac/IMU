function [is_ellipsoid, A,B,C] = is_valid_ellipsoid(ellipsoid_params, eps)
% 功能：根据文献判定拟合结果是否为椭球（二次项矩阵正定）
% 输入：ellipsoid_params-椭球参数，eps-最小正系数阈值（默认1e-6）
% 输出：is_ellipsoid-是否为椭球（逻辑值），A/B/C-标准形式系数
if nargin < 2
    eps = 1e-6;
end

% 提取二次项系数（文献方程1）
a = ellipsoid_params(1); b = ellipsoid_params(2); c = ellipsoid_params(3);
f = ellipsoid_params(4); g = ellipsoid_params(5); h = ellipsoid_params(6);

% 二次项矩阵（文献2.1节，需正定）
M = [a, h, g;
     h, b, f;
     g, f, c];
% 特征值判定（正定矩阵所有特征值>0）
eig_vals = eig(M);
is_ellipsoid = all(eig_vals > eps);

% 计算标准形式系数（用于后续验证）
if is_ellipsoid
    [V, Lambda] = eig(M);
    lambda = diag(Lambda);
    A = lambda(1); B = lambda(2); C = lambda(3);
else
    A = 0; B = 0; C = 0;
end
end

