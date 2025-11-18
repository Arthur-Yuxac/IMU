close all; clear; clc;

% 加载数据
load('IMU_1107_1549_Recorder1_20251107141152.mat');
% load('D:/C607/传感器/data/IMU_1031_1610_Recorder1_20251031160728/IMU_1031_1610_Recorder1_20251031160728.mat');

% 提取原始数据（保持原始结构）
accel_x = CANFDMessage__acceleration_mg_0(:, 2);  % 列向量
accel_y = CANFDMessage__acceleration_mg_1(:, 2);
accel_z = CANFDMessage__acceleration_mg_2(:, 2);

gyro_x = CANFDMessage__angular_rate_mdps_0(:, 2);
gyro_y = CANFDMessage__angular_rate_mdps_1(:, 2);
gyro_z = CANFDMessage__angular_rate_mdps_2(:, 2);

% 提取原始时间戳（转换为秒，列向量）
time_original = CANFDMessage__AMS330_Timestamp(:, 2) * 25e-6;

% --------------------------
% 处理时间戳异常（重复/非递增/间隔过大）
% --------------------------
% 1. 去重并确保时间递增（保留第一个出现的时间点，按时间排序）
[~, unique_idx] = unique(time_original, 'first');
sorted_idx = sort(unique_idx);  % 确保时间单调递增
time_unique = time_original(sorted_idx);  % 去重后时间（列向量）
% 同步提取对应的数据（列向量）
accel_x_u = accel_x(sorted_idx);
accel_y_u = accel_y(sorted_idx);
accel_z_u = accel_z(sorted_idx);
gyro_x_u = gyro_x(sorted_idx);
gyro_y_u = gyro_y(sorted_idx);
gyro_z_u = gyro_z(sorted_idx);

% 2. 计算时间间隔，过滤有效间隔（必须>0且不超过阈值）
dt = diff(time_unique);  % 相邻时间差（列向量）
threshold = 0.02;  % 异常间隔阈值
valid_dt = dt(dt > 0 & dt <= threshold);  % 有效间隔（>0且<=阈值）
abnormal_flag = ~isempty(dt) && (any(dt <= 0) || any(dt > threshold));  % 异常标志

% 3. 估计正常采样频率（基于有效间隔）
if isempty(dt)
    error('数据点不足（少于2个），无法处理');
end
if isempty(valid_dt)
    % 所有间隔无效时，用非负间隔估算（避免时间倒序导致的错误）
    non_neg_dt = dt(dt >= 0);
    if isempty(non_neg_dt)
        error('时间戳严重错乱（所有间隔为负），无法处理');
    end
    fs_estimated = 1 / mean(non_neg_dt);
else
    fs_estimated = 1 / mean(valid_dt);
end

% --------------------------
% 插值补充缺失数据（保持维度一致）
% --------------------------
if abnormal_flag
    % 存在异常：创建均匀时间轴（行向量，与原始时间范围一致）
    t_start = time_unique(1);
    t_end = time_unique(end);
    t_full = t_start:1/fs_estimated:t_end;  % 均匀时间（行向量，1×M）
    
    % 插值传感器数据（原始数据为列向量，需转置为行向量后插值）
    accel_x_interp = interp1(time_unique, accel_x_u, t_full, 'linear');  % 1×M
    accel_y_interp = interp1(time_unique, accel_y_u, t_full, 'linear');
    accel_z_interp = interp1(time_unique, accel_z_u, t_full, 'linear');
    gyro_x_interp = interp1(time_unique, gyro_x_u, t_full, 'linear');
    gyro_y_interp = interp1(time_unique, gyro_y_u, t_full, 'linear');
    gyro_z_interp = interp1(time_unique, gyro_z_u, t_full, 'linear');
else
    % 无异常：直接使用去重后的数据（转置为行向量，保持1×N）
    t_full = time_unique';  % 列向量转置为行向量
    accel_x_interp = accel_x_u';
    accel_y_interp = accel_y_u';
    accel_z_interp = accel_z_u';
    gyro_x_interp = gyro_x_u';
    gyro_y_interp = gyro_y_u';
    gyro_z_interp = gyro_z_u';
end

% --------------------------
% 单位转换与数据整合（保持4×N结构）
% --------------------------
% 时间起点归零（行向量，1×M）
time = t_full - t_full(1);

% 传感器数据单位转换（mg→g，mdps→dps），并组合为3×M矩阵
accel = [accel_x_interp; accel_y_interp; accel_z_interp] * 1e-3;  % 3×M
gyro = [gyro_x_interp; gyro_y_interp; gyro_z_interp] * 1e-3;      % 3×M

% 整合为4×N矩阵（第1行为时间，后3行为三个轴数据）
accel_with_error_time = [time; accel];  % 4×M
gyro_with_error_time = [time; gyro];    % 4×M

% --------------------------
% 输出处理信息
% --------------------------
fprintf('数据处理完成！\n');
fprintf('原始数据点数: %d\n', length(time_original));
fprintf('去重后数据点数: %d\n', length(time_unique));
fprintf('处理后数据点数: %d\n', size(accel_with_error_time, 2));  % 列数即N
fprintf('采样频率: %.2f Hz\n', fs_estimated);