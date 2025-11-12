clear; clc; close all;

load('Preprocessed_Data.mat');
total_sample = length(time_acc);
a_xp = a_xp(1:total_sample);
a_yp = a_yp(1:total_sample);
a_zp = a_zp(1:total_sample);

w_d = 2*(fs+1);
half_w_d = floor(w_d / 2);
max_times_the_var = 10;
init_static_samples = fs * 3;

var_x_init = var(a_xp(1:init_static_samples));
var_y_init = var(a_yp(1:init_static_samples));
var_z_init = var(a_zp(1:init_static_samples));
var_3D = var_x_init^2 + var_y_init^2 + var_z_init^2;

normal_x = zeros(1, total_sample);
normal_y = zeros(1, total_sample);
normal_z = zeros(1, total_sample);

for i = (half_w_d + 1) : (total_sample - half_w_d)
    win_x = a_xp(i - half_w_d : i + half_w_d);
    win_y = a_yp(i - half_w_d : i + half_w_d);
    win_z = a_zp(i - half_w_d : i + half_w_d);
    normal_x(i) = var(win_x);
    normal_y(i) = var(win_y);
    normal_z(i) = var(win_z);
end

s_square = normal_x.^2 + normal_y.^2 + normal_z.^2;

res_norm_vector = zeros(11, max_times_the_var);

for times_the_var = 1 : max_times_the_var
    s_filter_temp = zeros(1, total_sample);
    for i = half_w_d : (total_sample - half_w_d)
        if s_square(i) < times_the_var * var_3D
            s_filter_temp(i) = 1;
        end
    end

    static_indices = find(s_filter_temp == 1);
    if length(static_indices) < 9
        res_norm_vector(10, times_the_var) = Inf;
        continue;
    end
    selectedAccData = [a_xp(static_indices); a_yp(static_indices); a_zp(static_indices)];

    theta_pr_init = zeros(1, 9);
    ObjectiveFunction = @(theta) accCostFunct(theta, selectedAccData);
    [theta_pr, rsnorm] = lsqnonlin(ObjectiveFunction, theta_pr_init, [], [], optimset('MaxIter', 1000));

    res_norm_vector(1:9, times_the_var) = theta_pr;
    res_norm_vector(10, times_the_var) = rsnorm;
    res_norm_vector(11, times_the_var) = times_the_var * var_3D;
end

[min_res, min_idx] = min(res_norm_vector(10, :));
threshold_opt = 5 * res_norm_vector(11, min_idx);  %判断静态条件，可放缩
s_filter = zeros(1, total_sample);
for i = half_w_d : (total_sample - half_w_d)
    if s_square(i) < threshold_opt
        s_filter(i) = 1;
    end
end

save('Static_Detection_Result.mat', 's_filter', 'var_3D', 'threshold_opt', 's_square', 'time_acc');

figure('Name', '静态区间检测结果');
subplot(2,1,1);
plot(time_acc, s_square, 'b-');
hold on; plot(time_acc, threshold_opt * ones(1, total_sample), 'r--');
xlabel('时间 (s)'); ylabel('3D方差 magnitude');
legend('实时3D方差', '最优阈值'); title('方差与阈值对比');
subplot(2,1,2);
plot(time_acc, a_xp, 'b-'); hold on;
plot(time_acc, s_filter, 'r-');
xlabel('时间 (s)'); ylabel('加速度计x轴值 / 静态标记');
legend('加速度计x轴', '静态区间（放大500倍）'); title('静态区间标记');

disp('脚本2：静态区间检测完成，结果已保存至 Static_Detection_Result.mat');

function [residual] = accCostFunct(theta, selectedAccData)
%     g_mag = 9.81;
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