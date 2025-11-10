% 脚本2：静态区间检测（基于论文方差阈值法）
% 输入：Preprocessed_Data.mat - 脚本1输出的去零偏数据
% 输出：s_filter - 静态区间标记（1=静态，0=运动）
%       var_3D - 初始静态段3D方差（阈值基准）
%       threshold_opt - 最优静态检测阈值

clear; clc; close all;

%% 1. 加载预处理数据
load('Preprocessed_Data.mat');
total_sample = length(time_acc);
a_xp = a_xp(1:total_sample);
a_yp = a_yp(1:total_sample);
a_zp = a_zp(1:total_sample);
% load('sensor_data_with_errors.mat');
% a_xp = accel_with_error_time(2, 1:total_sample);
% a_yp = accel_with_error_time(3, 1:total_sample);
% a_zp = accel_with_error_time(4, 1:total_sample);

%% 2. 设定静态检测参数（对应论文Section IV-A）
w_d = 889; % 滑动窗口大小（建议为奇数，论文中t_w=2s，需根据采样率调整，此处假设采样率100Hz，窗口=101个样本≈1s）
half_w_d = floor(w_d / 2);
max_times_the_var = 10; % 阈值倍数遍历范围（论文中k=1~10，自动选最优）
init_static_samples = 3000; % 初始静态段样本数（用于计算var_3D，论文中T_init）

%% 3. 计算初始静态段3D方差（论文中的ς_init²）
var_x_init = var(a_xp(1:init_static_samples));
var_y_init = var(a_yp(1:init_static_samples));
var_z_init = var(a_zp(1:init_static_samples));
var_3D = var_x_init^2 + var_y_init^2 + var_z_init^2; % 3D方差magnitude（论文式(14)简化）

%% 4. 滑动窗口计算每个时刻的加速度方差
normal_x = zeros(1, total_sample); % x轴窗口方差
normal_y = zeros(1, total_sample); % y轴窗口方差
normal_z = zeros(1, total_sample); % z轴窗口方差

for i = (half_w_d + 1) : (total_sample - half_w_d)
    % 窗口内数据（中心在i）
    win_x = a_xp(i - half_w_d : i + half_w_d);
    win_y = a_yp(i - half_w_d : i + half_w_d);
    win_z = a_zp(i - half_w_d : i + half_w_d);
    % 计算窗口方差
    normal_x(i) = var(win_x);
    normal_y(i) = var(win_y);
    normal_z(i) = var(win_z);
end

% 计算3D方差magnitude（论文中的ς(t)²）
s_square = normal_x.^2 + normal_y.^2 + normal_z.^2;

%% 5. 遍历最优阈值（论文中自动选择最小残差对应的k）
res_norm_vector = zeros(11, max_times_the_var); % 存储参数与残差（前9为acc参数，10为残差，11为阈值）

for times_the_var = 1 : max_times_the_var
    % 临时静态标记（当前阈值下）
    s_filter_temp = zeros(1, total_sample);
    for i = half_w_d : (total_sample - half_w_d)
        if s_square(i) < times_the_var * var_3D
            s_filter_temp(i) = 1; % 标记为静态
        end
    end

    % 提取静态区间数据（用于加速度计校准残差计算，简化版论文式(10)）
    static_indices = find(s_filter_temp == 1);
    if length(static_indices) < 9 % 论文要求至少9个静态位置（Section IV-D）
        res_norm_vector(10, times_the_var) = Inf; % 静态数据不足，残差设为无穷大
        continue;
    end
    selectedAccData = [a_xp(static_indices); a_yp(static_indices); a_zp(static_indices)];

    % 加速度计参数初步优化（计算残差，论文式(10)的lsqnonlin简化）
    theta_pr_init = zeros(1, 9); % 加速度计参数θ^acc：[α_yz,α_zy,α_zx,s_x^a,s_y^a,s_z^a,b_x^a,b_y^a,b_z^a]
    ObjectiveFunction = @(theta) accCostFunct(theta, selectedAccData);
    [theta_pr, rsnorm] = lsqnonlin(ObjectiveFunction, theta_pr_init, [], [], optimset('MaxIter', 1000));

    % 存储结果
    res_norm_vector(1:9, times_the_var) = theta_pr;
    res_norm_vector(10, times_the_var) = rsnorm;
    res_norm_vector(11, times_the_var) = times_the_var * var_3D;
end

%% 6. 选择最优阈值（最小残差对应的阈值）
[min_res, min_idx] = min(res_norm_vector(10, :));
threshold_opt = 3 * res_norm_vector(11, min_idx);
s_filter = zeros(1, total_sample);
for i = half_w_d : (total_sample - half_w_d)
    if s_square(i) < threshold_opt
        s_filter(i) = 1;
    end
end

%% 7. 结果保存与可视化
save('Static_Detection_Result.mat', 's_filter', 'var_3D', 'threshold_opt', 's_square', 'time_acc');

% 绘制静态检测结果
figure('Name', '静态区间检测结果');
subplot(2,1,1);
plot(time_acc, s_square, 'b-');
hold on; plot(time_acc, threshold_opt * ones(1, total_sample), 'r--');
xlabel('时间 (s)'); ylabel('3D方差 magnitude');
legend('实时3D方差', '最优阈值'); title('方差与阈值对比');
subplot(2,1,2);
plot(time_acc, a_xp, 'b-'); hold on;
plot(time_acc, s_filter, 'r-'); % 放大静态标记以便观察
xlabel('时间 (s)'); ylabel('加速度计x轴值 / 静态标记');
legend('加速度计x轴', '静态区间（放大500倍）'); title('静态区间标记');

disp('脚本2：静态区间检测完成，结果已保存至 Static_Detection_Result.mat');

%% 辅助函数：加速度计代价函数（论文式(10)简化）
function [residual] = accCostFunct(theta, selectedAccData)
    % theta: [α_yz, α_zy, α_zx, s_x^a, s_y^a, s_z^a, b_x^a, b_y^a, b_z^a]
    % selectedAccData: 3*M 静态加速度数据（M为静态样本数）
    g_mag = 9.81; % 当地重力加速度magnitude（论文中||g||）
    M = size(selectedAccData, 2);
    residual = zeros(M, 1);

    % 论文式(6)：a^O = T^a * K^a * (a^S + b^a)（简化噪声项）
    T_a = [1, -theta(1), theta(2); 
           0, 1, -theta(3); 
           0, 0, 1]; % 加速度计对准矩阵T^a
    K_a = diag(theta(4:6)); % 加速度计比例因子矩阵K^a
    b_a = theta(7:9)'; % 加速度计偏置向量b^a

    for i = 1:M
        a_S = selectedAccData(:, i);
        a_O = T_a * K_a * (a_S + b_a); % 修正后的加速度（论文中a^O）
        residual(i) = g_mag^2 - norm(a_O)^2; % 论文式(10)的残差项
    end
end
