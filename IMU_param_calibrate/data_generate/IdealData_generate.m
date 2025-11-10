%% 初始化
close all; clear; clc;

%% 核心参数设置
% 采样频率与时间参数
fs = 888;                  % 采样频率(Hz)
dt = 1 / fs;               % 采样间隔(秒)
stable_hold_time = 30;     % 固定点停留时间(秒)
num_stable_hold = round(stable_hold_time * fs);  % 停留点数(确保30秒)
num_process_segments = 3*fs; % 两点间过渡段点数

% 椭球参数(单位球)
a = 1; b = 1; c = 1;

% 经线设置(方位角)
num_meridians = 8;
theta_deg = [0, 45, 90, 135, 180, 225, 270, 315];
theta = theta_deg * pi / 180;

% 固定点纬度(弧度)
phi_stable_deg = [0, 45, 90, 135, 180];
phi_stable = phi_stable_deg * pi / 180;


%% 1. 生成带时间戳的加速度计理想数据
% 预计算总数据点数
points_per_meridian = length(phi_stable)*num_stable_hold + ...
                      (length(phi_stable)-1)*num_process_segments;
total_points = num_meridians * points_per_meridian;
time_stamp = (0 : total_points-1) * dt;  % 全局时间戳

% 生成加速度计理想轨迹
accel_ideal = zeros(3, total_points);
current_idx = 1;  % 数据索引指针

for i = 1:num_meridians
    theta_i = theta(i);  % 当前经线方位角
    
    for p = 1:length(phi_stable)
        % 计算固定点坐标
        sin_phi = sin(phi_stable(p));
        cos_phi = cos(phi_stable(p));
        x_stable = b * sin_phi * cos(theta_i);
        y_stable = c * sin_phi * sin(theta_i);
        z_stable = a * cos_phi;
        stable_point = [x_stable; y_stable; z_stable];
        
        % 写入固定点数据(停留30秒)
        end_idx = current_idx + num_stable_hold - 1;
        accel_ideal(:, current_idx:end_idx) = repmat(stable_point, 1, num_stable_hold);
        current_idx = end_idx + 1;
        
        % 生成过渡段数据
        if p < length(phi_stable)
            phi_start = phi_stable(p);
            phi_end = phi_stable(p+1);
            phi_process = linspace(phi_start, phi_end, num_process_segments + 1);
            phi_process = phi_process(2:end);  % 排除起点
            
            for k = 1:length(phi_process)
                sin_phi_p = sin(phi_process(k));
                cos_phi_p = cos(phi_process(k));
                x_p = b * sin_phi_p * cos(theta_i);
                y_p = c * sin_phi_p * sin(theta_i);
                z_p = a * cos_phi_p;
                accel_ideal(:, current_idx) = [x_p; y_p; z_p];
                current_idx = current_idx + 1;
            end
        end
    end
end


%% 2. 生成带时间戳的陀螺仪理想数据
gyro_ideal = zeros(3, total_points);  % 3轴角速度(rad/s)
base_omega = 15;  % 基础旋转速率

for k = 2:total_points
    pos_curr = accel_ideal(:, k);
    pos_prev = accel_ideal(:, k-1);
    delta_pos = pos_curr - pos_prev;
    
    if norm(delta_pos) > 1e-6  % 运动段
        omega_dir = cross(pos_prev, pos_curr);
        omega_dir = omega_dir / norm(omega_dir);
        omega_mag = base_omega * norm(delta_pos) / dt;
        gyro_ideal(:, k) = omega_dir * omega_mag;
    else  % 停留段(角速度为0)
        gyro_ideal(:, k) = [0; 0; 0];
    end
end
gyro_ideal(:, 1) = gyro_ideal(:, 2);  % 第一个点与第二个点一致


%% 3. 生成带误差的传感器数据
% 3.1 加速度计误差参数
accel_bias = [0.02; 0.015; 0.01];          % 零偏(g)
accel_K = [1.01; 0.98; 1.012];             % 比例因子误差
accel_T = [1 0.005 0.003; 
           0 1     0.002; 
           0 0     1];  % 交叉耦合矩阵
accel_noise_std = 0.005;                    % 噪声标准差(g)

% 3.2 陀螺仪误差参数（符合MEMS陀螺典型特性）
gyro_bias = [0.05; 0.04; 0.06];            % 零偏(rad/s)
gyro_K = [1.02; 0.99; 1.015];              % 比例因子误差
gyro_T = [1 0.008 0.004; 
          0.003 1 0.005; 
          0.002 0.001 1];  % 交叉耦合矩阵
gyro_noise_std = 0.001;                    % 噪声标准差(rad/s)，比加速度计小一个量级
gyro_drift_rate = 0.001;                   % 零漂率(rad/s²)

% 3.3 计算带误差数据
accel_with_error = inv(accel_T) * inv(diag(accel_K)) *(accel_ideal + accel_bias);
%                    accel_noise_std * randn(3, total_points);

% 陀螺仪添加零漂（随时间缓慢变化的零偏）
gyro_drift = gyro_drift_rate * cumsum(dt * ones(3, total_points));  % 累积漂移
gyro_with_error = inv(gyro_T) * inv(diag(gyro_K))* (gyro_ideal + gyro_bias + gyro_drift);
%                   gyro_noise_std * randn(3, total_points);


%% 4. 数据整合(时间戳+数据)
accel_ideal_time = [time_stamp; accel_ideal];
accel_with_error_time = [time_stamp; accel_with_error];
gyro_ideal_time = [time_stamp; gyro_ideal];
gyro_with_error_time = [time_stamp; gyro_with_error];


%% 5. 可视化
% 5.1 3D轨迹对比
figure('Color', 'w', 'Position', [100 100 900 600]);
hold on; grid on; axis equal;
plot3(accel_ideal(1,:), accel_ideal(2,:), accel_ideal(3,:), 'b-', 'LineWidth', 0.8);
plot3(accel_with_error(1,:), accel_with_error(2,:), accel_with_error(3,:), 'r.', 'MarkerSize', 2);
% 标记固定点
for i = 1:num_meridians
    theta_i = theta(i);
    for p = 1:length(phi_stable)
        sin_phi = sin(phi_stable(p));
        cos_phi = cos(phi_stable(p));
        x_stable = b * sin_phi * cos(theta_i);
        y_stable = c * sin_phi * sin(theta_i);
        z_stable = a * cos_phi;
        plot3(x_stable, y_stable, z_stable, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    end
end
xlabel('X (g)'); ylabel('Y (g)'); zlabel('Z (g)');
title('加速度计轨迹(含30秒停留固定点)');
view(3);

% 5.2 加速度计各轴时序
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
sgtitle('加速度计各轴时序(888Hz，固定点停留30秒)');

% 5.3 陀螺仪各轴时序对比
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


%% 6. 数据保存
save('sensor_data_with_errors.mat', ...
     'accel_ideal_time', 'accel_with_error_time', ...
     'gyro_ideal_time', 'gyro_with_error_time', ...
     'time_stamp', 'fs', 'stable_hold_time');
fprintf('数据保存完成！总时长: %.2f秒，总点数: %d\n', max(time_stamp), total_points);