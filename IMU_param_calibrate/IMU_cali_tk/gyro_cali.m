clear; clc; close all;

load('Preprocessed_Data.mat');
total_sample = length(time_gyro);
omega_x = omega_x(1:total_sample) ;
omega_y = omega_y(1:total_sample);
omega_z = omega_z(1:total_sample);
time_gyro = time_gyro(1:total_sample);
dt = 1 / 250;

load('Static_Detection_Result.mat');
s_filter = s_filter(1:total_sample);
var_3D = var_3D;
time_acc = time_acc(1:total_sample);

load('Accelerometer_Calib_Result.mat');
acc_calib_params = acc_calib_result.acc_calib_params;
T_a = acc_calib_result.T_a;
K_a = acc_calib_result.K_a;
b_a = acc_calib_result.b_a;
calib_acc = acc_calib_result.calib_acc(:, 1:total_sample);

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

b_g_x = mean(omega_x(init_long_qs_start:init_long_qs_end));
b_g_y = mean(omega_y(init_long_qs_start:init_long_qs_end));
b_g_z = mean(omega_z(init_long_qs_start:init_long_qs_end));
b_g = [b_g_x; b_g_y; b_g_z];

omega_x_hat = omega_x - b_g_x;
omega_y_hat = omega_y - b_g_y;
omega_z_hat = omega_z - b_g_z;

fprintf('陀螺仪偏置估计完成：b_g = [%.6f, %.6f, %.6f]\n', b_g_x, b_g_y, b_g_z);

QS_dt_calib_info_matrix = zeros(6, 1);
l = 1;
flag_static = 0;
static_start = 0;
static_samples = 0;
num_samples_qs = ceil(1 / dt);

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
            QS_dt_calib_info_matrix(1:3, l) = [static_start; static_end; static_samples];
            avg_gravity = mean(calib_acc(:, static_start:static_end), 2);
            QS_dt_calib_info_matrix(4:6, l) = avg_gravity;
            l = l + 1;
        end
        flag_static = 0;
        static_samples = 0;
    end
end

QS_dt_calib_info_matrix(:, all(QS_dt_calib_info_matrix == 0)) = [];
num_static = size(QS_dt_calib_info_matrix, 2);
if num_static < 2
    error('有效静态区间仅%d个，无法满足代价函数计算要求', num_static);
end
fprintf('构建静态区间信息矩阵：%d列（与代价函数输入格式一致）\n', num_static);

theta_pr_gyro = [1, 0, 0, 0, 1, 0, 0, 0, 1];
lb = [0.95, -0.02, -0.02, -0.02, 0.95, -0.02, -0.02, -0.02, 0.95];
ub = [1.05,  0.02,  0.02,  0.02, 1.05,  0.02,  0.02,  0.02, 1.05];

option = optimset(...
    'TolX', 1e-7, ...
    'TolFun', 1e-6, ...
    'MaxFunEvals', 400, ...
    'Display', 'iter');


theta_pr_gyro = lsqnonlin(...
    @(theta_pr_gyro) gyroCostFunctLSQNONLINUsingOnlyTheFilter(theta_pr_gyro, QS_dt_calib_info_matrix, omega_x_hat, omega_y_hat, omega_z_hat, dt), ...
    theta_pr_gyro, ...
    [], ...
    [], ...
    option);

fprintf('陀螺仪参数优化完成，复合参数E = [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f]\n', ...
    theta_pr_gyro(1), theta_pr_gyro(2), theta_pr_gyro(3), theta_pr_gyro(4), theta_pr_gyro(5), ...
    theta_pr_gyro(6), theta_pr_gyro(7), theta_pr_gyro(8), theta_pr_gyro(9));

T_g = [1, theta_pr_gyro(2), theta_pr_gyro(3);
                theta_pr_gyro(4), 1, theta_pr_gyro(6);
                theta_pr_gyro(7), theta_pr_gyro(8), 1];

K_g = diag([theta_pr_gyro(1), theta_pr_gyro(5), theta_pr_gyro(9)]);

omega_hat = [omega_x_hat; omega_y_hat; omega_z_hat] * pi / 180;
calib_omega = T_g * K_g * omega_hat;

residuals = gyroCostFunctLSQNONLINUsingOnlyTheFilter(theta_pr_gyro, QS_dt_calib_info_matrix, omega_x_hat, omega_y_hat, omega_z_hat, dt);
avg_error_deg = mean(residuals);
max_error_deg = max(residuals);

fprintf('\n==================== 陀螺仪标定精度验证 ====================\n');
fprintf('重力方向平均误差：%.3f \n', avg_error_deg);
fprintf('重力方向最大误差：%.3f \n', max_error_deg);

gyro_calib_result = struct(...
    'gyro_calib_params', theta_pr_gyro, ...
    'T_g', T_g, ...
    'K_g', K_g, ...
    'b_g', b_g, ...
    'calib_omega', calib_omega, ...
    'residuals', residuals, ...
    'QS_dt_calib_info_matrix', QS_dt_calib_info_matrix, ...
    'dt', dt);

save('Gyroscope_Calib_Result.mat', 'gyro_calib_result', 'time_gyro');
fprintf('\n陀螺仪标定结果保存完成，文件：Gyroscope_Calib_Result.mat\n');

figure('Name', '陀螺仪标定结果（指定代价函数）', 'Position', [100, 100, 1000, 600]);

subplot(2,1,1);
plot(time_gyro, omega_x_hat, 'b-', 'LineWidth', 1); hold on;
plot(time_gyro, calib_omega(1, :), 'r-', 'LineWidth', 1.2);
for k = 1:num_static
    start_t = time_gyro(QS_dt_calib_info_matrix(1, k));
    end_t = time_gyro(QS_dt_calib_info_matrix(2, k));
    area([start_t, end_t], [min(omega_x_hat), max(omega_x_hat)], ...
         'FaceColor', [0.8, 0.8, 0.8], 'FaceAlpha', 0.3);
end
xlabel('时间 (s)'); ylabel('陀螺仪x轴角速度（单位与原始数据一致）');
legend('校准前（去偏后）', '校准后', '静态区间');
title('陀螺仪x轴校准前后对比（采用指定代价函数）');
grid on; grid minor;

subplot(2,1,2);
res_idx = 1:length(residuals);
plot(res_idx, residuals, 'b.-', 'MarkerSize', 10, 'LineWidth', 1);
xlabel('静态区间对序号'); ylabel('重力方向误差（度）');
title('陀螺仪标定残差验证（代价函数输出结果）');
grid on; grid minor;

function [res_vector] = gyroCostFunctLSQNONLINUsingOnlyTheFilter(E, QS_dt_info_matrix, omega_x_hat, omega_y_hat, omega_z_hat, dt)

omega_hat = [omega_x_hat; omega_y_hat; omega_z_hat];

misalignmentMatrix = [1, E(2), E(3); E(4), 1, E(6); E(7), E(8), 1];
scalingMatrix = diag([E(1), E(5), E(9)]);

omega_bar = misalignmentMatrix*scalingMatrix*omega_hat;

omega_x = omega_bar(1,:);
omega_y = omega_bar(2,:);
omega_z = omega_bar(3,:);

vector = zeros(3,5);

for pr = 1:size(QS_dt_info_matrix, 2) - 1
    
    vector((pr-1)*3 + 1:(pr)*3, 1) = QS_dt_info_matrix(4:6,pr);
    vector((pr-1)*3 + 1:(pr)*3, 5) = QS_dt_info_matrix(4:6,pr + 1);
    gyroUnbiasUncalibratedValues = [omega_x(QS_dt_info_matrix(2,pr) + 1:QS_dt_info_matrix(1,pr + 1) - 1); omega_y(QS_dt_info_matrix(2,pr) + 1:QS_dt_info_matrix(1,pr + 1) - 1); omega_z(QS_dt_info_matrix(2,pr) + 1:QS_dt_info_matrix(1,pr + 1) - 1)];
    R = rotationRK4(gyroUnbiasUncalibratedValues, dt); 

    vector((pr-1)*3 + 1:(pr)*3, 2:4) = R;
    
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



function [R] = rotationRK4(gyro_data, dt)
q = [1; 0; 0; 0];
num_samples = size(gyro_data, 2);

    for i = 1:num_samples
        w = gyro_data(:, i);
        Omega = [0, -w(1), -w(2), -w(3);
                w(1), 0, w(3), -w(2);
                w(2), -w(3), 0, w(1);
                w(3), w(2), -w(1), 0];
    
        k1 = 0.5 * Omega * q;
        k2 = 0.5 * Omega * (q + k1 * dt / 2);
        k3 = 0.5 * Omega * (q + k2 * dt / 2);
        k4 = 0.5 * Omega * (q + k3 * dt);
    
        q = q + (k1 + 2*k2 + 2*k3 + k4) * dt / 6;
        q = q / norm(q);
    end

w = q(1); x = q(2); y = q(3); z = q(4);
R = [1 - 2*y^2 - 2*z^2, 2*x*y - 2*z*w, 2*x*z + 2*y*w;
     2*x*y + 2*z*w, 1 - 2*x^2 - 2*z^2, 2*y*z - 2*x*w;
     2*x*z - 2*y*w, 2*y*z + 2*x*w, 1 - 2*x^2 - 2*y^2];
end

