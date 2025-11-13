clear; clc;

g = 9.81;
fs = 800;
dt = 1/fs;
static_duration = 10;
static_points = 20;
w_rot = 2*pi;

points_per_static = round(static_duration * fs);
rotation_points = 0.2*fs;

total_static_points = static_points * points_per_static;
num_rotations = static_points - 1;
total_rotation_points = num_rotations * rotation_points;
total_points = total_static_points + total_rotation_points;
time_stamp = (0 : total_points-1) * dt;

q = [1 0 0 0]';

accel_ideal = zeros(3, total_points);
gyro_ideal = zeros(3, total_points);
accel_with_error = zeros(3, total_points);
gyro_with_error = zeros(3, total_points);

axes = zeros(3, num_rotations);
for i = 1:num_rotations
    theta = 2 * pi * rand();
    phi = acos(2 * rand() - 1);
    axes(1, i) = sin(phi) * cos(theta);
    axes(2, i) = sin(phi) * sin(theta);
    axes(3, i) = cos(phi);
end

accel_bias = [0.02; 0.015; 0.01];
accel_K = [1.01; 0.98; 1.012];
accel_T = [1 0.005 0.003; 
           0 1     0.002; 
           0 0     1];
accel_noise_std = 0.005;

gyro_bias = [0.05; 0.04; 0.06];
gyro_K = [1.02; 0.99; 1.015];
gyro_T = [1 0.008 0.004; 
          0.003 1 0.005; 
          0.002 0.001 1];
gyro_noise_std = 0.001;
gyro_drift_rate = 0.00001;

gyro_drift = gyro_drift_rate * cumsum(dt * ones(3, total_points));

current_idx = 1;

for static_idx = 1:static_points
    fprintf('开始第 %d 个静态阶段，持续 %d 秒...\n', static_idx, static_duration);
    
    R = quat2rotm(q');
    accel_g_ideal = (R' * [0; 0; -1]);
    
    for i = 1:points_per_static
        accel_ideal(:, current_idx) = accel_g_ideal;
        gyro_ideal(:, current_idx) = [0; 0; 0];
        
        current_idx = current_idx + 1;
        if current_idx > total_points, break; end
    end
    
    if static_idx == static_points, break; end
    
    fprintf('开始第 %d 个旋转阶段...\n', static_idx);
    current_axis = axes(:, static_idx);
    
    for i = 1:rotation_points
        theta = w_rot * dt;
        
        q_rot = [cos(theta/2); current_axis * sin(theta/2)];
        
        q = quatmultiply(q', q_rot')';
        q = q / norm(q);
        
        R = quat2rotm(q');
        
        accel_ideal(:, current_idx) = (R' * [0; 0; -1]);
        gyro_ideal(:, current_idx) = current_axis * w_rot * 180 / pi;
        
        current_idx = current_idx + 1;
        if current_idx > total_points, break; end
    end
end

accel_noise = accel_noise_std * randn(3, total_points);
accel_with_error = inv(accel_T) * inv(diag(accel_K)) * accel_ideal + accel_bias + accel_noise;

% --------------------------修改处--------------------------
% 按文献模型修正陀螺仪误差施加顺序：ω^S = (T^g K^g)⁻¹ ω^O - (b^g + 漂移) + 噪声
gyro_noise = gyro_noise_std * randn(3, total_points);
gyro_error_matrix = gyro_T * diag(gyro_K);
gyro_with_error = inv(gyro_error_matrix) * gyro_ideal + (gyro_bias + gyro_drift) + gyro_noise;
% ----------------------------------------------------------

accel_ideal_time = [time_stamp; accel_ideal];
accel_with_error_time = [time_stamp; accel_with_error];
gyro_ideal_time = [time_stamp; gyro_ideal];
gyro_with_error_time = [time_stamp; gyro_with_error];

figure('Color', 'w', 'Position', [100 100 900 600]);
hold on; grid on; axis equal;
plot3(accel_ideal(1,:), accel_ideal(2,:), accel_ideal(3,:), 'b-', 'LineWidth', 0.8);
plot3(accel_with_error(1,:), accel_with_error(2,:), accel_with_error(3,:), 'r.', 'MarkerSize', 2);
xlabel('X (g)'); ylabel('Y (g)'); zlabel('Z (g)');
title('加速度计轨迹(含10秒停留固定点)');
view(3);

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

figure('Color', 'w', 'Position', [100 300 1000 600]);
subplot(3,1,1);
plot(time_stamp, gyro_ideal(1,:), 'b-', 'LineWidth', 0.8);
hold on; plot(time_stamp, gyro_with_error(1,:), 'r-', 'LineWidth', 0.5);
xlabel('时间(s)'); ylabel('X轴角速度(dps)'); grid on; legend('理想', '带误差');
subplot(3,1,2);
plot(time_stamp, gyro_ideal(2,:), 'g-', 'LineWidth', 0.8);
hold on; plot(time_stamp, gyro_with_error(2,:), 'r-', 'LineWidth', 0.5);
xlabel('时间(s)'); ylabel('Y轴角速度(dps)'); grid on;
subplot(3,1,3);
plot(time_stamp, gyro_ideal(3,:), 'b-', 'LineWidth', 0.8);
hold on; plot(time_stamp, gyro_with_error(3,:), 'r-', 'LineWidth', 0.5);
xlabel('时间(s)'); ylabel('Z轴角速度(dps)'); grid on;
sgtitle('陀螺仪各轴时序对比(含零漂和噪声)');

save('imu_rotation_data_with_errors.mat', ...
     'accel_ideal_time', 'accel_with_error_time', ...
     'gyro_ideal_time', 'gyro_with_error_time', ...
     'time_stamp', 'fs', 'static_duration');
fprintf('数据保存完成！总时长: %.2f秒，总点数: %d\n', max(time_stamp), total_points);

function q = quatmultiply(q1, q2)
    w1 = q1(1); x1 = q1(2); y1 = q1(3); z1 = q1(4);
    w2 = q2(1); x2 = q2(2); y2 = q2(3); z2 = q2(4);
    
    q = [w1*w2 - x1*x2 - y1*y2 - z1*z2;
         w1*x2 + x1*w2 + y1*z2 - z1*y2;
         w1*y2 - x1*z2 + y1*w2 + z1*x2;
         w1*z2 + x1*y2 - y1*x2 + z1*w2]';
end

function R = quat2rotm(q)
    w = q(1); x = q(2); y = q(3); z = q(4);
    
    R = [1-2*y^2-2*z^2,   2*x*y-2*z*w,   2*x*z+2*y*w;
         2*x*y+2*z*w,   1-2*x^2-2*z^2,   2*y*z-2*x*w;
         2*x*z-2*y*w,     2*y*z+2*x*w, 1-2*x^2-2*y^2];
end