close all;clear;clc;
load('IMU_1107_1549_Recorder1_20251107141152.mat');
accel_x = CANFDMessage__acceleration_mg_0(:, 2);
accel_y = CANFDMessage__acceleration_mg_1(:, 2);
accel_z = CANFDMessage__acceleration_mg_2(:, 2);


gyro_x = CANFDMessage__angular_rate_mdps_0(:, 2);
gyro_y = CANFDMessage__angular_rate_mdps_1(:, 2);
gyro_z = CANFDMessage__angular_rate_mdps_2(:, 2);


time = CANFDMessage__AMS330_Timestamp(:, 2) * 25e-6;
dt_origin = diff(time);
N = length(time);
abnormal_idx = find(dt_origin > 0.02) ; 

for idx = 2 : N
    if ismember(idx, abnormal_idx)
        time(idx) = (time(idx-1) + time(idx+1)) / 2;
        accel_x(idx) = (accel_x(idx-1) + accel_x(idx+1)) / 2;
        accel_y(idx) = (accel_y(idx-1) + accel_y(idx+1)) / 2;
        accel_z(idx) = (accel_z(idx-1) + accel_z(idx+1)) / 2;
        gyro_x(idx) = (gyro_x(idx-1) + gyro_x(idx+1)) / 2;
        gyro_y(idx) = (gyro_y(idx-1) + gyro_y(idx+1)) / 2;
        gyro_z(idx) = (gyro_z(idx-1) + gyro_z(idx+1)) / 2;
    end
end

accel = [accel_x'; accel_y'; accel_z'] * 1e-3;
gyro = [gyro_x'; gyro_y'; gyro_z'] * 1e-3;

time = time(2:end) - time(2);
gyro = gyro(:, 2:end);
accel = accel(:, 2:end);
fs = 1 / mean(diff(time));
accel_with_error_time = [time'; accel];
gyro_with_error_time = [time'; gyro];




