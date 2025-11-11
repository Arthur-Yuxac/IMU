clear; clc; close all;

load('Preprocessed_Data.mat');
load('Static_Detection_Result.mat');
total_sample = length(time_acc);
a_xp = a_xp(1:total_sample);
a_yp = a_yp(1:total_sample);
a_zp = a_zp(1:total_sample);

static_indices = find(s_filter == 1);
if length(static_indices) < 100
    error('静态样本不足，无法完成加速度计校准，请重新采集数据或调整静态检测参数');
end
selectedAccData = [a_xp(static_indices); a_yp(static_indices); a_zp(static_indices)];
M = size(selectedAccData, 2);

theta0 = zeros(1, 9);
theta0(4:6) = 1e-6;

ObjectiveFunction = @(theta) accCalibCost(theta, selectedAccData);

options = optimset(...
    'MaxFunEvals', 150000, ...
    'MaxIter', 6000, ...
    'TolFun', 1e-10, ...
    'Display', 'iter');
[acc_calib_params, resnorm] = lsqnonlin(ObjectiveFunction, theta0, [], [], options);

alpha_yz = acc_calib_params(1);
alpha_zy = acc_calib_params(2);
alpha_zx = acc_calib_params(3);
s_x_a = acc_calib_params(4);
s_y_a = acc_calib_params(5);
s_z_a = acc_calib_params(6);
b_x_a = acc_calib_params(7);
b_y_a = acc_calib_params(8);
b_z_a = acc_calib_params(9);

T_a = [1, -alpha_yz, alpha_zy;
       0, 1, -alpha_zx;
       0, 0, 1];
K_a = diag([s_x_a, s_y_a, s_z_a]);
b_a = [b_x_a; b_y_a; b_z_a];

a_S_matrix = [a_xp; a_yp; a_zp];
calib_acc = T_a * K_a * (a_S_matrix - b_a * ones(1, total_sample));

static_calib_acc = calib_acc(:, static_indices);
static_mag = sqrt(sum(static_calib_acc.^2, 1));
avg_static_mag = mean(static_mag);
fprintf('静态时校准后加速度平均magnitude：%.4f g（理论值：1）\n', avg_static_mag);

acc_calib_result = struct(...
    'acc_calib_params', acc_calib_params, ...
    'T_a', T_a, ...
    'K_a', K_a, ...
    'b_a', b_a, ...
    'calib_acc', calib_acc, ...
    'static_mag', static_mag);
save('Accelerometer_Calib_Result.mat', 'acc_calib_result', 'time_acc');

disp('脚本3：加速度计校准完成，结果已保存至 Accelerometer_Calib_Result.mat');

function [residual] = accCalibCost(theta, selectedAccData)
    g_mag = 9.81;
    M = size(selectedAccData, 2);
    residual = zeros(M, 1);

    T_a = [1, -theta(1), theta(2); 
           0, 1, -theta(3); 
           0, 0, 1];
    K_a = diag(theta(4:6));
    b_a = theta(7:9)';

    for i = 1:M
        a_S = selectedAccData(:, i);
        a_O = T_a * K_a * (a_S - b_a);
        residual(i) = 1^2 - norm(a_O)^2;
    end
end

