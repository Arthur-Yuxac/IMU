% 脚本1：IMU数据预处理（去除硬件零偏）
% 输入：accel_with_error_time(4*N) - 加速度计数据（行1：时间，行2-4：x/y/z原始值）
%       gyro_with_error_time(4*N) - 陀螺仪数据（行1：时间，行2-4：x/y/z原始值）
% 输出：a_xp/a_yp/a_zp - 去零偏后的加速度计信号
%       omega_x/omega_y/omega_z - 去零偏后的陀螺仪信号
%       time_acc/time_gyro - 对应的时间序列

% clear; clc; close all;
% load('sensor_data_with_errors');
load('imu_rotation_data_with_errors.mat');
%% 1. 输入数据提取
% 加速度计数据（行1：时间，行2-4：x/y/z原始输出）
time_acc = accel_with_error_time(1, :);
accel_raw_x = accel_with_error_time(2, :);
accel_raw_y = accel_with_error_time(3, :);
accel_raw_z = accel_with_error_time(4, :);

% 陀螺仪数据（行1：时间，行2-4：x/y/z原始输出）
time_gyro = gyro_with_error_time(1, :);
gyro_raw_x = gyro_with_error_time(2, :);
gyro_raw_y = gyro_with_error_time(3, :);
gyro_raw_z = gyro_with_error_time(4, :);

% 确保时间序列长度一致（论文假设数据同步采集）
total_sample = min(length(time_acc), length(time_gyro));
time_acc = time_acc(1:total_sample);
time_gyro = time_gyro(1:total_sample);
accel_raw_x = accel_raw_x(1:total_sample);
accel_raw_y = accel_raw_y(1:total_sample);
accel_raw_z = accel_raw_z(1:total_sample);
gyro_raw_x = gyro_raw_x(1:total_sample);
gyro_raw_y = gyro_raw_y(1:total_sample);
gyro_raw_z = gyro_raw_z(1:total_sample);

%% 2. 设定硬件零偏（需根据传感器手册或预实验测量，论文中为传感器出厂默认零位）
% 示例值（需根据实际传感器调整，单位：g, dps）
offset_acc_x = 0.02; % 加速度计x轴零偏
offset_acc_y = 0.015; % 加速度计y轴零偏
offset_acc_z = 0.01; % 加速度计z轴零偏
offset_gyro_x = 0.05; % 陀螺仪x轴零偏
offset_gyro_y = 0.04; % 陀螺仪y轴零偏
offset_gyro_z = 0.06; % 陀螺仪z轴零偏

%% 3. 去除硬件零偏（对应论文中a^S、ω^S的定义，即去零偏后的原始信号）
% a_xp = accel_raw_x - offset_acc_x;
% a_yp = accel_raw_y - offset_acc_y;
% a_zp = accel_raw_z - offset_acc_z;
% omega_x = gyro_raw_x - offset_gyro_x;
% omega_y = gyro_raw_y - offset_gyro_y;
% omega_z = gyro_raw_z - offset_gyro_z;
a_xp = accel_raw_x;
a_yp = accel_raw_y;
a_zp = accel_raw_z;
omega_x = gyro_raw_x;
omega_y = gyro_raw_y;
omega_z = gyro_raw_z;

%% 4. 结果保存与可视化
save('Preprocessed_Data.mat', 'a_xp', 'a_yp', 'a_zp', 'omega_x', 'omega_y', 'omega_z', 'time_acc', 'time_gyro');

