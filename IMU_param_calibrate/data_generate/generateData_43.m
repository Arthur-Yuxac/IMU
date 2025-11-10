% IMU绕中心点连续旋转模拟 - 含误差模型与时间戳的数据生成
clear; clc;

% 参数设置
g = 9.81;               % 重力加速度 (m/s²)
fs = 888;               % 采样频率 (Hz)
dt = 1/fs;              % 采样间隔 (s)
static_duration = 10;   % 每个静态点停留时间 (s)
static_points = 20;     % 静态点总数
w_rot = 0.5*pi;           % 转动角速率 (rad/s) - 1转/秒

% 计算每个阶段的采样点数
points_per_static = round(static_duration * fs);  % 每个静态点的采样数
rotation_points = 1.5 * fs;                            % 假设每次旋转0.5s

% 总采样点数计算
total_static_points = static_points * points_per_static;
num_rotations = static_points - 1;               % 静态点之间的旋转次数
total_rotation_points = num_rotations * rotation_points;
total_points = total_static_points + total_rotation_points;  % 总采样点数
time_stamp = (0 : total_points-1) * dt;          % 时间戳

% 初始方向（四元数表示）
q = [1 0 0 0]';         % 初始四元数 [w x y z]，表示无旋转

% 预分配数据数组
accel_ideal = zeros(3, total_points);   % 加速度计理想数据(g)
gyro_ideal = zeros(3, total_points);    % 陀螺仪理想数据(rad/s)
accel_with_error = zeros(3, total_points);  % 带误差的加速度数据
gyro_with_error = zeros(3, total_points);   % 带误差的陀螺仪数据

% 生成均匀分布在球面上的转动轴（静态点之间的旋转轴）
axes = zeros(3, num_rotations);
for i = 1:num_rotations
    % 球坐标生成均匀分布的转动轴
    theta = 2 * pi * rand();         % 方位角 [0, 2π]
    phi = acos(2 * rand() - 1);      % 极角 [0, π]
    % 转换为笛卡尔坐标
    axes(1, i) = sin(phi) * cos(theta);
    axes(2, i) = sin(phi) * sin(theta);
    axes(3, i) = cos(phi);
end

% 误差模型参数（参考规范设置）
% 加速度计误差参数
accel_bias = [0.02; 0.015; 0.01];          % 零偏(g)
accel_K = [1.01; 0.98; 1.012];             % 比例因子误差
accel_T = [1 0.005 0.003; 
           0 1     0.002; 
           0 0     1];  % 交叉耦合矩阵
accel_noise_std = 0.005;                    % 噪声标准差(g)

% 陀螺仪误差参数
gyro_bias = [0.05; 0.04; 0.06];            % 零偏(rad/s)
gyro_K = [1.02; 0.99; 1.015];              % 比例因子误差
gyro_T = [1 0.008 0.004; 
          0.003 1 0.005; 
          0.002 0.001 1];  % 交叉耦合矩阵
gyro_noise_std = 0.001;                    % 噪声标准差(rad/s)
gyro_drift_rate = 0.001;                   % 零漂率(rad/s²)

% 生成陀螺仪零漂（随时间累积）
gyro_drift = gyro_drift_rate * cumsum(dt * ones(3, total_points));

% 生成数据
current_idx = 1;  % 当前数据索引

for static_idx = 1:static_points
    % 静态阶段
    fprintf('开始第 %d 个静态阶段，持续 %d 秒...\n', static_idx, static_duration);
    
    % 获取当前姿态的旋转矩阵
    R = quat2rotm(q');
    % 计算理想重力加速度（转换为g单位）
    accel_g_ideal = (R * [0; 0; -1]);  % 单位：g（已除以g）
    
    % 填充静态阶段数据
    for i = 1:points_per_static
        % 加速度计理想数据（单位：g）
        accel_ideal(:, current_idx) = accel_g_ideal;
        % 陀螺仪理想数据（静态时为0）
        gyro_ideal(:, current_idx) = [0; 0; 0];
        
        current_idx = current_idx + 1;
        if current_idx > total_points, break; end
    end
    
    % 最后一个静态点后不需要旋转
    if static_idx == static_points, break; end
    
    % 旋转阶段
    fprintf('开始第 %d 个旋转阶段...\n', static_idx);
    current_axis = axes(:, static_idx);  % 当前旋转轴
    
    for i = 1:rotation_points
        % 计算旋转角度
        theta = w_rot * dt;
        
        % 计算旋转四元数
        q_rot = [cos(theta/2); current_axis * sin(theta/2)];
        
        % 更新姿态四元数
        q = quatmultiply(q', q_rot')';
        q = q / norm(q);  % 归一化
        
        % 计算旋转矩阵
        R = quat2rotm(q');
        
        % 加速度计理想数据（单位：g）
        accel_ideal(:, current_idx) = (R * [0; 0; -1]);  % 单位：g
        % 陀螺仪理想数据（旋转轴×角速度）
        gyro_ideal(:, current_idx) = current_axis * w_rot;
        
        current_idx = current_idx + 1;
        if current_idx > total_points, break; end
    end
end

% 应用误差模型
% 加速度计误差计算
accel_noise = accel_noise_std * randn(3, total_points);
accel_with_error = inv(accel_T) * inv(diag(accel_K)) * (accel_ideal + accel_bias) + accel_noise;

% 陀螺仪误差计算
gyro_noise = gyro_noise_std * randn(3, total_points);
gyro_with_error = inv(gyro_T) * inv(diag(gyro_K)) * (gyro_ideal + gyro_bias + gyro_drift) + gyro_noise;

% 数据整合（时间戳+数据）
accel_ideal_time = [time_stamp; accel_ideal];
accel_with_error_time = [time_stamp; accel_with_error];
gyro_ideal_time = [time_stamp; gyro_ideal];
gyro_with_error_time = [time_stamp; gyro_with_error];

% 可视化
% 3D轨迹对比
figure('Color', 'w', 'Position', [100 100 900 600]);
hold on; grid on; axis equal;
plot3(accel_ideal(1,:), accel_ideal(2,:), accel_ideal(3,:), 'b-', 'LineWidth', 0.8);
plot3(accel_with_error(1,:), accel_with_error(2,:), accel_with_error(3,:), 'r.', 'MarkerSize', 2);
xlabel('X (g)'); ylabel('Y (g)'); zlabel('Z (g)');
title('加速度计轨迹(含10秒停留固定点)');
view(3);

% 加速度计各轴时序
figure('Color', 'w', 'Position', [100 200 1000 600]);
subplot(3,1,1);
plot(time_stamp, accel_ideal(1,:), 'b-', 'LineWidth', 0.8);
hold on; plot(time_stamp, accel_with_error(1,:), 'r-', 'LineWidth', 0.5);
xlabel('时间(s)'); ylabel('X轴加速度(g)'); grid on; legend('理想', '带误差');
subplot(3,1,2);
plot(time_stamp, accel_ideal(2,:), 'b-', 'LineWidth', 0.8);
hold on; plot(time_stamp, accel_with_error(2,:), 'r-', 'LineWidth', 0.5);
xlabel('时间(s)'); ylabel('Y轴加速度(g)'); grid on;
subplot(3,1,3);
plot(time_stamp, accel_ideal(3,:), 'b-', 'LineWidth', 0.8);
hold on; plot(time_stamp, accel_with_error(3,:), 'r-', 'LineWidth', 0.5);
xlabel('时间(s)'); ylabel('Z轴加速度(g)'); grid on;
sgtitle('加速度计各轴时序(888Hz，固定点停留10秒)');

% 陀螺仪各轴时序对比
figure('Color', 'w', 'Position', [100 300 1000 600]);
subplot(3,1,1);
plot(time_stamp, gyro_ideal(1,:), 'b-', 'LineWidth', 0.8);
hold on; plot(time_stamp, gyro_with_error(1,:), 'r-', 'LineWidth', 0.5);
xlabel('时间(s)'); ylabel('X轴角速度(rad/s)'); grid on; legend('理想', '带误差');
subplot(3,1,2);
plot(time_stamp, gyro_ideal(2,:), 'g-', 'LineWidth', 0.8);
hold on; plot(time_stamp, gyro_with_error(2,:), 'r-', 'LineWidth', 0.5);
xlabel('时间(s)'); ylabel('Y轴角速度(rad/s)'); grid on;
subplot(3,1,3);
plot(time_stamp, gyro_ideal(3,:), 'b-', 'LineWidth', 0.8);
hold on; plot(time_stamp, gyro_with_error(3,:), 'r-', 'LineWidth', 0.5);
xlabel('时间(s)'); ylabel('Z轴角速度(rad/s)'); grid on;
sgtitle('陀螺仪各轴时序对比(含零漂和噪声)');

% 数据保存
save('imu_rotation_data_with_errors.mat', ...
     'accel_ideal_time', 'accel_with_error_time', ...
     'gyro_ideal_time', 'gyro_with_error_time', ...
     'time_stamp', 'fs', 'static_duration');
fprintf('数据保存完成！总时长: %.2f秒，总点数: %d\n', max(time_stamp), total_points);

% 辅助函数：四元数乘法
function q = quatmultiply(q1, q2)
    w1 = q1(1); x1 = q1(2); y1 = q1(3); z1 = q1(4);
    w2 = q2(1); x2 = q2(2); y2 = q2(3); z2 = q2(4);
    
    q = [w1*w2 - x1*x2 - y1*y2 - z1*z2;
         w1*x2 + x1*w2 + y1*z2 - z1*y2;
         w1*y2 - x1*z2 + y1*w2 + z1*x2;
         w1*z2 + x1*y2 - y1*x2 + z1*w2]';
end

% 辅助函数：四元数转旋转矩阵
function R = quat2rotm(q)
    w = q(1); x = q(2); y = q(3); z = q(4);
    
    R = [1-2*y^2-2*z^2,   2*x*y-2*z*w,   2*x*z+2*y*w;
         2*x*y+2*z*w,   1-2*x^2-2*z^2,   2*y*z-2*x*w;
         2*x*z-2*y*w,     2*y*z+2*x*w, 1-2*x^2-2*y^2];
end