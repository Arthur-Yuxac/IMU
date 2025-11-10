%% 椭球法加速度计标定（基于静态区间检测结果）
close all; clear; clc;

%% 1. 加载数据和静态检测结果
% 加载原始加速度数据（单位：g，格式：[时间戳; x; y; z]）
load('sensor_data_with_errors.mat', 'accel_with_error_time', 'time_stamp');
accel_raw = accel_with_error_time(2:4, :)';  % 提取3轴数据，转为[N×3]矩阵（每行[x,y,z]）
time = time_stamp';  % 时间戳（N×1）

% 加载静态区间检测结果（s_filter：1=静态，0=运动）
load('Static_Detection_Result.mat', 's_filter');
s_filter = s_filter';  % 转为[N×1]向量，与数据长度匹配

% 检查数据长度一致性
% if length(s_filter) ~= size(accel_raw, 1)
%     error('静态标记与加速度数据长度不匹配！');
% end


%% 2. 从静态区间提取有效数据（用于椭球拟合）
% 提取静态区间索引
static_indices = find(s_filter == 1);
if isempty(static_indices)
    error('未检测到静态区间，请检查静态检测结果！');
end
fprintf('共检测到 %d 个静态样本点\n', length(static_indices));

% 提取静态数据并去重（保留不同姿态的稳定段均值）
% 方法：按时间分段，每段连续静态数据取均值（避免同一姿态重复采样）
static_data = accel_raw(static_indices, :);  % 静态样本点
static_time = time(static_indices);          % 静态点时间戳

% 分段取均值（相邻静态点时间间隔超过0.5秒视为新姿态）
segment_threshold = 0.5;  % 姿态切换时间阈值（秒）
diff_time = diff(static_time);
segment_flags = [1; diff_time > segment_threshold];  % 分段标记（1=新分段起点）
segment_ids = cumsum(segment_flags);  % 分段编号
num_segments = max(segment_ids);      % 总姿态数

% 计算每个静态段的均值（每个姿态一个均值点）
static_means = zeros(num_segments, 3);
for i = 1:num_segments
    seg_idx = segment_ids == i;
    static_means(i, :) = mean(static_data(seg_idx, :), 1);
end
fprintf('提取不同姿态的静态段均值：%d 组（用于椭球拟合）\n', num_segments);

if num_segments < 6
    warning('静态姿态数量较少（建议≥6），可能影响标定精度！');
end


%% 3. 椭球方程拟合（核心步骤）
% 静态均值点满足椭球方程：Ax² + By² + Cz² + 2Fyz + 2Gxz + 2Hxy + 2Px + 2Qy + 2Rz + D = 0
x = static_means(:, 1);
y = static_means(:, 2);
z = static_means(:, 3);
N = size(static_means, 1);

k_init = 1000;
k_min = 3;
k_tol = 0.01;
[optimal_k, ellipsoid_params] = search_optimal_k(static_means, k_init, k_min, k_tol);

a = ellipsoid_params(1); b = ellipsoid_params(2); c = ellipsoid_params(3);
f = ellipsoid_params(4); g = ellipsoid_params(5); h = ellipsoid_params(6);
p = ellipsoid_params(7); q = ellipsoid_params(8); r = ellipsoid_params(9);
d = ellipsoid_params(10);



%% 4. 从椭球参数提取加速度计误差参数
% 误差模型：测量值 = inv(T) * inv(diag(K)) * (理想值 + b) + 噪声
% 理想值满足：理想值·理想值 = 1（单位球，g=1）

% 4.1 计算椭球中心（测量值空间的零偏映射）
center_x = -p / a;      % 椭球中心x坐标
center_y = -q / b;      % 椭球中心y坐标
center_z = -r / c;  % 椭球中心z坐标
z0 = [center_x; center_y; center_z];  % 中心向量

% 4.2 二次项矩阵特征分解（主轴方向与缩放）
M = [a, h, g; h, b, f; g, f, c];
% M = [a, h, g; 0, b, f; 0, 0, c];
[V, Lambda] = eig(M);  % V：特征向量矩阵（主轴方向），Lambda：特征值对角阵
lambda = diag(Lambda); % 特征值

% 4.3 计算比例因子K和交叉耦合矩阵T
D_prime = z0' * M * z0 + d;  % 中心代入后的常数项
scale_factor = sqrt(-lambda / D_prime);  % 比例因子（K的对角元素）
coupling_matrix = V;                     % 交叉耦合矩阵（主轴方向）

% 4.4 计算物理零偏b（转换至理想值空间）
K = diag(scale_factor);
b = K * coupling_matrix * z0;  % 零偏向量（g）


%% 5. 输出标定结果
fprintf('\n=== 加速度计椭球法标定结果 ===\n');
fprintf('零偏（b）：\n');
fprintf('b_x = %.6f g\n', b(1));
fprintf('b_y = %.6f g\n', b(2));
fprintf('b_z = %.6f g\n', b(3));

fprintf('\n比例因子（K）：\n');
fprintf('K_x = %.6f\n', scale_factor(1));
fprintf('K_y = %.6f\n', scale_factor(2));
fprintf('K_z = %.6f\n', scale_factor(3));

fprintf('\n交叉耦合矩阵（T）：\n');
fprintf('[%+.6f  %+.6f  %+.6f]\n', coupling_matrix(1, :));
fprintf('[%+.6f  %+.6f  %+.6f]\n', coupling_matrix(2, :));
fprintf('[%+.6f  %+.6f  %+.6f]\n', coupling_matrix(3, :));


%% 6. 标定效果验证
% 6.1 修正静态均值数据
z0_mat = repmat(z0', num_segments, 1);  % 椭球中心矩阵
static_corrected = (static_means - z0_mat) * coupling_matrix' * K;  % 修正为理想值

% 6.2 计算修正前后的模长（理想值应为1g）
raw_norm_mean = mean(sqrt(sum(static_means.^2, 2)));
corr_norm_mean = mean(sqrt(sum(static_corrected.^2, 2)));

fprintf('\n=== 标定效果验证 ===\n');
fprintf('静态均值原始模长均值：%.6f g\n', raw_norm_mean);
fprintf('修正后模长均值：%.6f g（目标：1g）\n', corr_norm_mean);
fprintf('相对误差：%.4f%%\n', abs(corr_norm_mean - 1) / 1 * 100);


%% 7. 可视化标定效果
figure('Position', [100, 100, 800, 700]);
hold on; grid on; axis equal;

% 原始静态均值点
scatter3(static_means(:, 1), static_means(:, 2), static_means(:, 3), ...
         80, 'r', 'filled', 'DisplayName', '原始静态均值');

% 拟合椭球表面
[x_ell, y_ell, z_ell] = generate_ellipsoid(ellipsoid_params, 50);
mesh(x_ell, y_ell, z_ell, 'FaceAlpha', 0.3, 'EdgeColor', 'b', ...
     'DisplayName', '拟合椭球');

% 修正后静态均值点
scatter3(static_corrected(:, 1), static_corrected(:, 2), static_corrected(:, 3), ...
         80, 'g', 'filled', 'DisplayName', '修正后静态均值');

% 理想单位球（1g）
[xs, ys, zs] = sphere(50);
surf(xs, ys, zs, 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'DisplayName', '理想单位球');

xlabel('X (g)'); ylabel('Y (g)'); zlabel('Z (g)');
title('基于静态区间的椭球法标定效果');
legend('Location', 'best');
view(3);


%% 辅助函数：生成椭球表面
function [x, y, z] = generate_ellipsoid(params, n)
    A = params(1); B = params(2); C = params(3);
    F = params(4); G = params(5); H = params(6);
    P = params(7); Q = params(8); R = params(9); D = params(10);
    
    % 椭球中心
    x0 = -P / A;
    y0 = -Q / B;
    z0 = -R / C;
    
    % 二次项矩阵对角化
    M = [A, H, G; H, B, F; G, F, C];
    [V, Lambda] = eig(M);
    lambda = diag(Lambda);
    
    % 半轴长度
    a_semi = sqrt(-D / lambda(1));
    b_semi = sqrt(-D / lambda(2));
    c_semi = sqrt(-D / lambda(3));
    
    % 球面映射为椭球
    [xs, ys, zs] = sphere(n);
    xs = a_semi * xs; ys = b_semi * ys; zs = c_semi * zs;
    xyz = [xs(:), ys(:), zs(:)] * V';  % 主轴变换
    x = reshape(xyz(:, 1), n+1, n+1) + x0;
    y = reshape(xyz(:, 2), n+1, n+1) + y0;
    z = reshape(xyz(:, 3), n+1, n+1) + z0;
end