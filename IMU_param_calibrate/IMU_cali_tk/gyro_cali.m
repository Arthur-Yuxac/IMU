% 脚本4：陀螺仪标定（指定代价函数+求解方法，与加速度标定数据结构一致）
% 核心：采用您提供的gyroCostFunctLSQNONLINUsingOnlyTheFilter代价函数，lsqnonlin直接调用求解
% 数据结构：完全对齐加速度标定（参数向量、矩阵格式、结果存储1:1匹配）

clear; clc; close all;

%% 1. 加载前期数据（与加速度标定共用数据源，变量名完全一致）
% -------------------------- 脚本1：预处理数据 --------------------------
load('Preprocessed_Data.mat'); 
total_sample = length(time_gyro);
% 陀螺仪数据（去硬件零偏后，对应加速度计a_xp/a_yp/a_zp）
omega_x = omega_x(1:total_sample); 
omega_y = omega_y(1:total_sample); 
omega_z = omega_z(1:total_sample);
time_gyro = time_gyro(1:total_sample);
time_interval = 0.01; % 采样周期（与您代价函数默认值一致，对应100Hz）

% -------------------------- 脚本2：静态检测结果 --------------------------
load('Static_Detection_Result.mat'); 
s_filter = s_filter(1:total_sample); % 静态标记（与加速度标定共用）
var_3D = var_3D;
time_acc = time_acc(1:total_sample);

% -------------------------- 脚本3：加速度计校准结果（参考基准） --------------------------
load('Accelerometer_Calib_Result.mat'); 
% 加速度计校准参数（保留原始结构，与陀螺仪标定参数格式对齐）
acc_calib_params = acc_calib_result.acc_calib_params; % 加速度计9参数向量
T_a = acc_calib_result.T_a;       % 轴对准矩阵（T^a）
K_a = acc_calib_result.K_a;       % 比例因子矩阵（K^a）
b_a = acc_calib_result.b_a;       % 偏置向量（b^a）
calib_acc = acc_calib_result.calib_acc(:, 1:total_sample); % 校准后加速度（a^O）


%% 2. 步骤1：陀螺仪偏置去除（类比加速度计偏置估计逻辑）
% 提取第一个完整长静态段（与您之前的脚本逻辑一致）
init_long_qs_start = 0;
init_long_qs_end = 0;
flag_is_first_long_static = 1;

for i = 1:total_sample
    if s_filter(i) == 0 && flag_is_first_long_static == 1
    elseif s_filter(i) == 1 && flag_is_first_long_static == 1
        init_long_qs_start = i;
        flag_is_first_long_static = 2;
    elseif s_filter(i) == 1 && flag_is_first_long_static == 2
    elseif s_filter(i) == 0 && flag_is_first_long_static == 2
        init_long_qs_end = i - 1;
        break;
    end
end

% 估计偏置（对应加速度计b_a，变量名b_g对齐b_a）
b_g_x = mean(omega_x(init_long_qs_start:init_long_qs_end));
b_g_y = mean(omega_y(init_long_qs_start:init_long_qs_end));
b_g_z = mean(omega_z(init_long_qs_start:init_long_qs_end));
b_g = [b_g_x; b_g_y; b_g_z]; % 陀螺仪偏置向量（3×1，与b_a格式一致）

% 去偏后数据（对应您代价函数中的omega_x_hat/omega_y_hat/omega_z_hat）
omega_x_hat = omega_x - b_g_x;
omega_y_hat = omega_y - b_g_y;
omega_z_hat = omega_z - b_g_z;

fprintf('陀螺仪偏置估计完成：b_g = [%.6f, %.6f, %.6f]\n', b_g_x, b_g_y, b_g_z);


%% 3. 步骤2：构建静态区间信息矩阵（QS_time_interval_calib_info_matrix，与您指定格式一致）
% 矩阵格式：6行N列（[起始,结束,样本数, u_a_x, u_a_y, u_a_z]），完全匹配您代价函数的输入要求
QS_time_interval_calib_info_matrix = zeros(6, 1); 
l = 1;
flag_static = 0;
static_start = 0;
static_samples = 0;
num_samples_qs = ceil(1 / time_interval); % 每个静态段最小样本数（1秒）

for i = 1:total_sample
    if flag_static == 0 && s_filter(i) == 1
        static_start = i;
        static_samples = 1;
        flag_static = 1;
    elseif flag_static == 1 && s_filter(i) == 1
        static_samples = static_samples + 1;
    elseif flag_static == 1 && s_filter(i) == 0
        if static_samples >= num_samples_qs
            static_end = i - 1;
            % 前3行：静态段起始、结束、样本数
            QS_time_interval_calib_info_matrix(1:3, l) = [static_start; static_end; static_samples];
            % 后3行：静态段平均重力向量（校准后加速度归一化）
            avg_gravity = mean(calib_acc(:, static_start:static_end), 2);
            QS_time_interval_calib_info_matrix(4:6, l) = avg_gravity;
            l = l + 1;
        end
        flag_static = 0;
        static_samples = 0;
    end
end

% 筛选有效区间（确保≥2个，满足代价函数循环要求）
QS_time_interval_calib_info_matrix(:, all(QS_time_interval_calib_info_matrix == 0)) = [];
num_static = size(QS_time_interval_calib_info_matrix, 2);
if num_static < 2
    error('有效静态区间仅%d个，无法满足代价函数计算要求', num_static);
end
fprintf('构建静态区间信息矩阵：%d列（与代价函数输入格式一致）\n', num_static);


%% 4. 步骤3：陀螺仪参数优化（采用您指定的代价函数与求解方法）
% 4.1 待估参数向量E（复合参数：尺度+轴对准，对应您代价函数的E，与加速度计参数维度一致）
% E = [s_x, γ_yz, γ_zy, γ_xz, s_y, γ_zx, γ_xy, γ_yx, s_z]（9个参数，对齐acc_calib_params）
theta_pr_gyro = [1, 0, 0, 0, 1, 0, 0, 0, 1]; % 初始值（理想无误差，与您脚本初始值逻辑一致）

% 4.2 优化配置（与您指定的option参数完全一致）
% option = optimset('TolX', 1e-7, 'TolFun', 1e-6, 'MaxFunEvals', 400);
% 4.2 优化配置（添加Display字段显示迭代过程）
option = optimset(...
    'TolX', 1e-7, ...
    'TolFun', 1e-6, ...
    'MaxFunEvals', 400, ...
    'Display', 'iter'); % 新增：显示迭代过程（iter=显示每一步迭代信息）

% 4.3 调用lsqnonlin求解（完全采用您指定的调用格式）
% 目标函数：@(theta_pr_gyro) gyroCostFunctLSQNONLINUsingOnlyTheFilter(...)
theta_pr_gyro = lsqnonlin(...
    @(theta_pr_gyro) gyroCostFunctLSQNONLINUsingOnlyTheFilter(theta_pr_gyro, QS_time_interval_calib_info_matrix, omega_x_hat, omega_y_hat, omega_z_hat, time_interval), ...
    theta_pr_gyro, ...
    [], ...
    [], ...
    option);

fprintf('陀螺仪参数优化完成，复合参数E = [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f]\n', ...
    theta_pr_gyro(1), theta_pr_gyro(2), theta_pr_gyro(3), theta_pr_gyro(4), theta_pr_gyro(5), ...
    theta_pr_gyro(6), theta_pr_gyro(7), theta_pr_gyro(8), theta_pr_gyro(9));


%% 5. 步骤4：提取校准矩阵（与加速度计T_a/K_a格式完全对齐）
% 从复合参数E中拆分轴对准矩阵与比例因子矩阵（对应您代价函数的矩阵定义）
% 5.1 轴对准矩阵misal_matrix（对应加速度计T_a，3×3格式一致）
misal_matrix = [1, theta_pr_gyro(2), theta_pr_gyro(3);
                theta_pr_gyro(4), 1, theta_pr_gyro(6);
                theta_pr_gyro(7), theta_pr_gyro(8), 1]; % 与您代价函数中misalignmentMatrix定义一致

% 5.2 比例因子矩阵scale_matrix（对应加速度计K_a，3×3对角矩阵格式一致）
scale_matrix = diag([theta_pr_gyro(1), theta_pr_gyro(5), theta_pr_gyro(9)]); % 与您代价函数中scalingMatrix定义一致

% 5.3 校准后角速度（对应加速度计calib_acc，计算逻辑与您代价函数中omega_bar一致）
omega_hat = [omega_x_hat; omega_y_hat; omega_z_hat];
calib_omega = misal_matrix * scale_matrix * omega_hat; % 3×N，与calib_acc维度一致


%% 6. 步骤5：精度验证（复用您代价函数的残差计算逻辑）
% 调用代价函数计算残差（重力方向误差）
residuals = gyroCostFunctLSQNONLINUsingOnlyTheFilter(theta_pr_gyro, QS_time_interval_calib_info_matrix, omega_x_hat, omega_y_hat, omega_z_hat, time_interval);
avg_error_deg = mean(residuals) * (180/pi); % 平均角度误差（度）
max_error_deg = max(residuals) * (180/pi); % 最大角度误差（度）

fprintf('\n==================== 陀螺仪标定精度验证 ====================\n');
fprintf('重力方向平均误差：%.3f 度\n', avg_error_deg);
fprintf('重力方向最大误差：%.3f 度\n', max_error_deg);


%% 7. 步骤6：结果保存（结构体字段与加速度标定完全对齐）
% 陀螺仪标定结果结构体（gyro_calib_result对应acc_calib_result，字段1:1匹配）
gyro_calib_result = struct(...
    'gyro_calib_params', theta_pr_gyro, ... % 对应acc_calib_params（9参数向量）
    'T_g', misal_matrix, ...                 % 对应T_a（轴对准矩阵）
    'K_g', scale_matrix, ...                 % 对应K_a（比例因子矩阵）
    'b_g', b_g, ...                           % 对应b_a（偏置向量）
    'calib_omega', calib_omega, ...           % 对应calib_acc（校准后数据）
    'residuals', residuals, ...               % 对应加速度计static_mag（验证残差）
    'QS_time_interval_calib_info_matrix', QS_time_interval_calib_info_matrix, ... % 静态区间信息
    'time_interval', time_interval);          % 对应加速度计sample_period

% 保存结果（文件名与加速度标定风格一致）
save('Gyroscope_Calib_Result.mat', 'gyro_calib_result', 'time_gyro');
fprintf('\n陀螺仪标定结果保存完成，文件：Gyroscope_Calib_Result.mat\n');


%% 8. 步骤7：可视化（与加速度标定绘图风格一致）
figure('Name', '陀螺仪标定结果（指定代价函数）', 'Position', [100, 100, 1000, 600]);

% 子图1：校准前后陀螺仪x轴对比（类比加速度计x轴对比）
subplot(2,1,1);
plot(time_gyro, omega_x_hat, 'b-', 'LineWidth', 1); hold on;
plot(time_gyro, calib_omega(1, :), 'r-', 'LineWidth', 1.2);
% 标记静态区间（灰色背景，与加速度标定一致）
for k = 1:num_static
    start_t = time_gyro(QS_time_interval_calib_info_matrix(1, k));
    end_t = time_gyro(QS_time_interval_calib_info_matrix(2, k));
    area([start_t, end_t], [min(omega_x_hat), max(omega_x_hat)], ...
         'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.3);
end
xlabel('时间 (s)'); ylabel('陀螺仪x轴角速度（单位与原始数据一致）');
legend('校准前（去偏后）', '校准后', '静态区间');
title('陀螺仪x轴校准前后对比（采用指定代价函数）');
grid on; grid minor;

% 子图2：残差验证（类比加速度计静态模值误差）
subplot(2,1,2);
res_idx = 1:length(residuals);
plot(res_idx, residuals*(180/pi), 'b.-', 'MarkerSize', 10, 'LineWidth', 1);
xlabel('静态区间对序号'); ylabel('重力方向误差（度）');
title('陀螺仪标定残差验证（代价函数输出结果）');
grid on; grid minor;


%% 代价函数
function [res_vector] = gyroCostFunctLSQNONLINUsingOnlyTheFilter(E, QS_time_interval_info_matrix, omega_x_hat, omega_y_hat, omega_z_hat, time_interval)
if(nargin<6)
    time_interval=0.01;
end
omega_hat = [omega_x_hat; omega_y_hat; omega_z_hat];

misalignmentMatrix = [1, E(2), E(3); E(4), 1, E(6); E(7), E(8), 1];
scalingMatrix = diag([E(1), E(5), E(9)]);

omega_bar = misalignmentMatrix*scalingMatrix*omega_hat; % Huai, because the bias is removed beforehand, to obtain the calibrated omega

omega_x = omega_bar(1,:);
omega_y = omega_bar(2,:);
omega_z = omega_bar(3,:);

vector = zeros(3,5);

for pr = 1:size(QS_time_interval_info_matrix, 2) - 1
    
    vector((pr-1)*3 + 1:(pr)*3, 1) = QS_time_interval_info_matrix(4:6,pr); % the previous gravity vector 
    vector((pr-1)*3 + 1:(pr)*3, 5) = QS_time_interval_info_matrix(4:6,pr + 1); % the current gravity vector
    % extract the dynamic angular rate data
    gyroUnbiasUncalibratedValues = [omega_x(QS_time_interval_info_matrix(2,pr) + 1:QS_time_interval_info_matrix(1,pr + 1) - 1); omega_y(QS_time_interval_info_matrix(2,pr) + 1:QS_time_interval_info_matrix(1,pr + 1) - 1); omega_z(QS_time_interval_info_matrix(2,pr) + 1:QS_time_interval_info_matrix(1,pr + 1) - 1)];
    R = rotationRK4(gyroUnbiasUncalibratedValues, time_interval); 

    vector((pr-1)*3 + 1:(pr)*3, 2:4) = R; % the rotation from the previous to current gravity vector
    
end

residuals = zeros(length(vector(:,1))/3, 1);

for i = 1:length(vector(:,1))/3    
    v = vector((i-1)*3 + 1:(i)*3, 5)/(vector((i-1)*3 + 1, 5)^2 + vector((i-1)*3 + 2, 5)^2 + vector((i)*3, 5)^2)^(1/2) -...
    vector((i-1)*3 + 1:(i)*3, 2:4)*vector((i-1)*3 + 1:(i)*3, 1)/(vector((i-1)*3 + 1, 5)^2 + vector((i-1)*3 + 2, 5)^2 + vector((i)*3, 5)^2)^(1/2);
    v = v';
    residuals(i,1) = (v(1)^2 + v(2)^2 + v(3)^2)^(1/2);    
end
res_vector = residuals;
end


%% 配套的rotationRK4函数（与您代价函数调用逻辑一致）
function [R] = rotationRK4(gyro_data, time_interval)
% 四阶龙格-库塔积分计算旋转矩阵（与您代价函数中的调用要求匹配）
q = [1; 0; 0; 0]; % 初始四元数（无旋转）
num_samples = size(gyro_data, 2);

for i = 1:num_samples
    omega = gyro_data(:, i);
    % 角速度反对称矩阵（匹配您代价函数中的旋转逻辑）
    Omega = [0, -omega(1), -omega(2), -omega(3);
             omega(1), 0, omega(3), -omega(2);
             omega(2), -omega(3), 0, omega(1);
             omega(3), omega(2), -omega(1), 0];
    
    % RK4积分步骤
    k1 = 0.5 * Omega * q;
    k2 = 0.5 * Omega * (q + k1 * time_interval / 2);
    k3 = 0.5 * Omega * (q + k2 * time_interval / 2);
    k4 = 0.5 * Omega * (q + k3 * time_interval);
    
    % 更新并归一化四元数
    q = q + (k1 + 2*k2 + 2*k3 + k4) * time_interval / 6;
    q = q / norm(q);
end

% 四元数转旋转矩阵（3×3，匹配您代价函数中R的输出格式）
w = q(1); x = q(2); y = q(3); z = q(4);
R = [1 - 2*y^2 - 2*z^2, 2*x*y - 2*z*w, 2*x*z + 2*y*w;
     2*x*y + 2*z*w, 1 - 2*x^2 - 2*z^2, 2*y*z - 2*x*w;
     2*x*z - 2*y*w, 2*y*z + 2*x*w, 1 - 2*x^2 - 2*y^2];
end