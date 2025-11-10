function [optimal_k, ellipsoid_params] = search_optimal_k(stable_means, k_init, k_min, tol)
% 功能：迭代搜索最优k值（文献4.1节算法）
% 输入：stable_means-稳定段均值（g），k_init-初始大k值，k_min-最小k值，tol-收敛阈值
% 输出：optimal_k-最优k值，ellipsoid_params-对应椭球参数

% 初始搜索：从大k减小，找到拟合为椭球的k上界
current_k = k_init;
found = false;
while current_k > k_min && ~found
    [ellipsoid_params, ~] = fit_ellipsoid_ls(stable_means, current_k);
    [is_ellipsoid, ~, ~] = is_valid_ellipsoid(ellipsoid_params);
    if is_ellipsoid
        found = true;
        k_upper = current_k;
        k_lower = current_k / 2;
    else
        current_k = current_k / 2;  % 按文献逐步减半k
    end
end
if ~found
    error('文献方法：未找到有效k值，建议增加稳定段数量或优化数据');
end

% 二分法细化最优k（文献4.1节）
while (k_upper - k_lower) > tol
    k_mid = (k_upper + k_lower) / 2;
    [ellipsoid_params_mid, ~] = fit_ellipsoid_ls(stable_means, k_mid);
    [is_ellipsoid_mid, ~, ~] = is_valid_ellipsoid(ellipsoid_params_mid);
    if is_ellipsoid_mid
        k_upper = k_mid;
    else
        k_lower = k_mid;
    end
end
optimal_k = k_upper;
[ellipsoid_params, ~] = fit_ellipsoid_ls(stable_means, optimal_k);
fprintf('文献方法迭代搜索得到最优k值：%.2f\n', optimal_k);
end
