% 脚本3：加速度计校准（基于论文多位置最小二乘）
% 输入：Preprocessed_Data.mat（脚本1）、Static_Detection_Result.mat（脚本2）
% 输出：acc_calib_params - 加速度计校准参数（θ^acc）
%       T_a - 加速度计对准矩阵（论文式(6)）
%       K_a - 加速度计比例因子矩阵（论文式(6)）
%       b_a - 加速度计偏置向量（论文式(6)）
%       calib_acc - 校准后的加速度计数据

clear; clc; close all;

%% 1. 加载前期数据
load('Preprocessed_Data.mat');
load('Static_Detection_Result.mat');
total_sample = length(time_acc);
a_xp = a_xp(1:total_sample);
a_yp = a_yp(1:total_sample);
a_zp = a_zp(1:total_sample);

%% 2. 提取静态区间数据（论文要求至少9个不同静态位置，Section IV-D）
static_indices = find(s_filter == 1);
if length(static_indices) < 100 % 确保有足够静态样本（建议>100）
    error('静态样本不足，无法完成加速度计校准，请重新采集数据或调整静态检测参数');
end
selectedAccData = [a_xp(static_indices); a_yp(static_indices); a_zp(static_indices)];
M = size(selectedAccData, 2); % 静态样本数

%% 3. 定义加速度计校准参数与代价函数（论文式(8)、式(10)）
% 待估计参数θ^acc：[α_yz, α_zy, α_zx, s_x^a, s_y^a, s_z^a, b_x^a, b_y^a, b_z^a]
theta0 = zeros(1, 9); % 初始值（论文建议从理想值开始）
theta0(4:6) = 1e-6; % 比例因子初始值（需根据传感器量程调整，示例：1e-6 m/(s²·AD计数)）

% 代价函数：最小化论文式(10) L(θ^acc) = sum(||g||² - ||a^O||²)²
ObjectiveFunction = @(theta) accCalibCost(theta, selectedAccData);

%% 4. 最小二乘优化（论文用Levenberg-Marquardt算法，对应lsqnonlin）
options = optimset(...
    'MaxFunEvals', 150000, ...
    'MaxIter', 6000, ...
    'TolFun', 1e-10, ...
    'Display', 'iter'); % 显示迭代过程（可选）
[acc_calib_params, resnorm] = lsqnonlin(ObjectiveFunction, theta0, [], [], options);

%% 5. 提取加速度计校准矩阵（论文式(6)）
alpha_yz = acc_calib_params(1);
alpha_zy = acc_calib_params(2);
alpha_zx = acc_calib_params(3);
s_x_a = acc_calib_params(4);
s_y_a = acc_calib_params(5);
s_z_a = acc_calib_params(6);
b_x_a = acc_calib_params(7);
b_y_a = acc_calib_params(8);
b_z_a = acc_calib_params(9);

% 对准矩阵T^a（论文式(6)）
T_a = [1, -alpha_yz, alpha_zy;
       0, 1, -alpha_zx;
       0, 0, 1];
% 比例因子矩阵K^a（论文式(6)）
K_a = diag([s_x_a, s_y_a, s_z_a]);
% 偏置向量b_a（论文式(6)）
b_a = [b_x_a; b_y_a; b_z_a];

%% 6. 计算校准后的加速度计数据（论文式(6)：a^O = T^a*K^a*(a^S + b_a)）
a_S_matrix = [a_xp; a_yp; a_zp]; % 3*N 去零偏原始信号
calib_acc = T_a * K_a * (a_S_matrix + b_a * ones(1, total_sample)); % 校准后加速度（单位：m/s²）

%% 7. 结果验证（静态时校准后加速度magnitude应接近9.81）
static_calib_acc = calib_acc(:, static_indices);
static_mag = sqrt(sum(static_calib_acc.^2, 1));
avg_static_mag = mean(static_mag);
fprintf('静态时校准后加速度平均magnitude：%.4f m/s²（理论值：9.81 m/s²）\n', avg_static_mag);

%% 8. 结果保存与可视化
acc_calib_result = struct(...
    'acc_calib_params', acc_calib_params, ...
    'T_a', T_a, ...
    'K_a', K_a, ...
    'b_a', b_a, ...
    'calib_acc', calib_acc, ...
    'static_mag', static_mag);
save('Accelerometer_Calib_Result.mat', 'acc_calib_result', 'time_acc');

% 绘制校准前后静态加速度magnitude对比
% figure('Name', '加速度计校准结果');
% subplot(2,1,1);
% plot(time_acc(static_indices), static_mag, 'b.');
% hold on; plot(time_acc(static_indices), ones(1, length(static_indices)), 'r--');
% xlabel('时间 (s)'); ylabel('加速度magnitude (m/s²)');
% legend('校准后静态加速度', '理论重力值(9.81)'); title('静态加速度校准验证');
% subplot(2,1,2);
% plot(time_acc, calib_acc(1,:), 'b-');
% xlabel('时间 (s)'); ylabel('校准后加速度x轴 (m/s²)');
% title('校准后加速度计x轴数据');

disp('脚本3：加速度计校准完成，结果已保存至 Accelerometer_Calib_Result.mat');

%% 辅助函数：加速度计校准代价函数（论文式(10)）
function [residual] = accCalibCost(theta, selectedAccData)
    % theta: 加速度计参数θ^acc
    % selectedAccData: 3*M 静态加速度原始数据
    g_mag = 9.81; % 当地重力加速度magnitude（论文中||g||）
    M = size(selectedAccData, 2);
    residual = zeros(M, 1);

    % 提取参数并构建矩阵（论文式(6)）
    T_a = [1, -theta(1), theta(2); 
           0, 1, -theta(3); 
           0, 0, 1];
    K_a = diag(theta(4:6));
    b_a = theta(7:9)';

    % 计算每个静态样本的残差（论文式(10)）
    for i = 1:M
        a_S = selectedAccData(:, i);
        a_O = T_a * K_a * (a_S - b_a); % 修正后的加速度a^O
        residual(i) = 1^2 - norm(a_O)^2; % 残差项：||g||² - ||a^O||²
    end
end
