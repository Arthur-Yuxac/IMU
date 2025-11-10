close all; clear; clc;


%% 2. 迭代搜索最优k值（文献4.1节）
k_init = 5;    % 初始大k值（文献建议）
k_min = 3;       % 最小k值（文献中k>3保证约束有效性）
k_tol = 0.01;    % k值收敛阈值
[optimal_k, ellipsoid_params] = search_optimal_k(stable_means, k_init, k_min, k_tol);
% optimal_k = 4;

% 输出椭球参数（对应文献方程1）
fprintf('最优k值下的椭球参数：\n');
fprintf('a=%.4e, b=%.4e, c=%.4e, f=%.4e, g=%.4e, h=%.4e, p=%.4e, q=%.4e, r=%.4e, d=%.4e\n', ...
    ellipsoid_params(1:10));

%% 3. 反解加速度计误差参数（单位：g）
% 提取椭球参数（文献方程1系数）
a = ellipsoid_params(1); b = ellipsoid_params(2); c = ellipsoid_params(3);
f = ellipsoid_params(4); g = ellipsoid_params(5); h = ellipsoid_params(6);
p = ellipsoid_params(7); q = ellipsoid_params(8); r = ellipsoid_params(9);
d = ellipsoid_params(10);
gravity = 1;  % 重力加速度（单位：g，修正后模长目标）

% 3.1 计算零偏（椭球中心映射，文献2.2节）
b_x = -p / a;  % X轴零偏（g）
b_y = -q / b;  % Y轴零偏（g）
b_z = -r / c;  % Z轴零偏（g）

% 3.2 计算比例因子与耦合系数（文献3.2节）
M_ellipsoid = [a, h, g; h, b, f; g, f, c];  % 二次项矩阵
[V, Lambda] = eig(M_ellipsoid);% 差系数
lambda = diag(Lambda);

% 中心代入后的常数项（文献推导）
D_prime = a*b_x^2 + b*b_y^2 + c*b_z^2 + 2*f*b_y*b_z + 2*g*b_x*b_z + 2*h*b_x*b_y + ...
          2*p*b_x + 2*q*b_y + 2*r*b_z + d;

% 比例因子（无单位，修正比例）
S_xx = sqrt(gravity^2 * lambda(1) / (-D_prime));
S_yy = sqrt(gravity^2 * lambda(2) / (-D_prime));
S_zz = sqrt(gravity^2 * lambda(3) / (-D_prime));

% 轴间耦合系数（主轴方向矩阵非对角元素，文献3.2节）
S_xy = V(1,2); S_xz = V(1,3);
S_yx = V(2,1); S_yz = V(2,3);
S_zx = V(3,1); S_zy = V(3,2);

% 输出误差参数
fprintf('\n=== 加速度计标定结果（单位：g） ===\n');
fprintf('零偏：b_x=%.4f, b_y=%.4f, b_z=%.4f\n', b_x, b_y, b_z);
fprintf('比例因子：S_xx=%.4f, S_yy=%.4f, S_zz=%.4f\n', S_xx, S_yy, S_zz);
fprintf('轴间耦合系数：\nS_xy=%.4f, S_xz=%.4f\nS_yx=%.4f, S_yz=%.4f\nS_zx=%.4f, S_zy=%.4f\n', ...
    S_xy, S_xz, S_yx, S_yz, S_zx, S_zy);

%% 4. 误差补偿（修正原始数据，单位：g）
% 原始数据（mg）转换为g
accelRaw_g = accelRaw / 1000;  % accelRaw为原始加速度数据（mg）

% 补偿公式：修正值 = S^{-1}*(原始值 - 零偏)
b_vec = [b_x, b_y, b_z];  % 零偏向量（g）
b_vec_fix = [0.0021297, -0.0104362, 0.0200418];
S_diag = diag([S_xx, S_yy, S_zz]);  % 比例因子对角矩阵
S = S_diag + [0 S_xy S_xz;
              S_yx 0 S_yz;
              S_zx S_zy 0;];
accel_corrected = (accelRaw_g - ones(size(accelRaw_g,1),1)*b_vec) / S;  % 修正后数据（g）
accel_corrected_fix = (accelRaw_g - ones(size(accelRaw_g, 1), 1) * b_vec_fix) / S;

%% 5. 标定效果验证（目标：模长均值≈1g，文献4.2节实验标准）
% 5.1 稳定段均值修正后模长验证
mean_norm_raw = mean(sqrt(sum(stable_means.^2, 2)));
norm_error_raw = abs(mean_norm_raw - gravity) / gravity * 100;
stable_means_corrected = (stable_means - ones(size(stable_means,1),1)*b_vec) / S_diag;
mean_norm = mean(sqrt(sum(stable_means_corrected.^2, 2)));
norm_error = abs(mean_norm - gravity) / gravity * 100;  % 相对误差（%）
stable_means_corrected_fix = (stable_means - ones(size(stable_means,1),1)*b_vec_fix) / S_diag;
mean_norm_fix = mean(sqrt(sum(stable_means_corrected_fix.^2, 2)));
norm_error_fix = abs(mean_norm_fix - gravity) / gravity * 100;  % 相对误差（%）
fprintf('\n=== 标定效果验证 ===\n');
fprintf('原始平均段均值模长均值： %.4f g(目标：1g) \n', mean_norm_raw);
fprintf('模长相对误差：%.2f %%（文献实验优秀标准：<1%%）\n', norm_error_raw);
fprintf('修正后稳定段均值模长均值：%.4f g（目标：1 g）\n', mean_norm);
fprintf('模长相对误差：%.2f %%（文献实验优秀标准：<1%%）\n', norm_error);
fprintf('替换零偏修正后稳定段均值模长均值：%.4f g（目标：1 g）\n', mean_norm_fix);
fprintf('替换零偏后模长相对误差：%.2f %%（文献实验优秀标准：<1%%）\n', norm_error_fix);

% 5.2 三维可视化（文献图2-7风格）
figure('Position', [100, 100, 800, 700]);
% 原始稳定段均值（g）
scatter3(stable_means(:,1), stable_means(:,2), stable_means(:,3), 100, 'r', 'filled', 'DisplayName', '原始稳定段均值（g）');
hold on;
% 拟合椭球（文献方法生成）
[x_ell, y_ell, z_ell] = generate_ellipsoid_surface(ellipsoid_params, 50);  % 50为网格密度
mesh(x_ell, y_ell, z_ell, 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'DisplayName', sprintf('拟合椭球（k=%.2f）', optimal_k));
% 修正后稳定段均值（g）
scatter3(stable_means_corrected(:,1), stable_means_corrected(:,2), stable_means_corrected(:,3), 100, 'g', 'filled', 'DisplayName', '修正后稳定段均值（g）');
% 理想重力球面（半径1g）
[xs, ys, zs] = sphere(50);
xs = xs * gravity; ys = ys * gravity; zs = zs * gravity;
mesh(xs, ys, zs, 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'DisplayName', '理想重力球面（1g）');

xlabel('X轴加速度（g）'); ylabel('Y轴加速度（g）'); zlabel('Z轴加速度（g）');
title('基于《Least_squares_ellipsoid_specific_fitting.pdf》的IMU标定效果');
legend('Location', 'best');
grid on; view(3);


IMU_params = struct();
IMU_params.accel.bias = b_vec';% g
IMU_params.accel.scale_mat = [S_xx S_xy S_xz;
                              S_yx S_yy S_yz;
                              S_zx S_zy S_zz];
IMU_params.gyro.bias = [0.7345;   -1.5181;   -0.5765];% mds
IMU_params.gyro.scale_mat = [1 0 0;
                             0 1 0;
                             0 0 1];
save('IMU_params.mat', 'IMU_params');


%% 辅助函数：生成椭球表面（文献方法）
function [x, y, z] = generate_ellipsoid_surface(ellipsoid_params, n)
% 输入：ellipsoid_params-椭球参数，n-网格密度
% 输出：椭球表面坐标（单位：g）
a = ellipsoid_params(1); b = ellipsoid_params(2); c = ellipsoid_params(3);
f = ellipsoid_params(4); g = ellipsoid_params(5); h = ellipsoid_params(6);
p = ellipsoid_params(7); q = ellipsoid_params(8); r = ellipsoid_params(9);
d = ellipsoid_params(10);

% 椭球中心（文献2.1节）
x0 = -p/a; y0 = -q/b; z0 = -r/c;
% 二次项矩阵对角化（文献3.1节）
M = [a, h, g; h, b, f; g, f, c];
[V, Lambda] = eig(M);
lambda = diag(Lambda);
% 椭球半轴长度（文献推导）
a_semi = sqrt(-d / lambda(1));
b_semi = sqrt(-d / lambda(2));
c_semi = sqrt(-d / lambda(3));

% 球面网格映射为椭球（文献图生成逻辑）
[xs, ys, zs] = sphere(n);
xs = a_semi * xs; ys = b_semi * ys; zs = c_semi * zs;
xyz = [xs(:), ys(:), zs(:)] * V';  % 主轴变换
x = reshape(xyz(:,1), n+1, n+1) + x0;
y = reshape(xyz(:,2), n+1, n+1) + y0;
z = reshape(xyz(:,3), n+1, n+1) + z0;
end
