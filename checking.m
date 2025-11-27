close all;clear;clc;

load('euler_ideal.mat');
% load('euler_est_ukf.mat');
% load('euler_est_eskf.mat');
load('euler_est_ekf.mat');
% 6. 与理想姿态对比分析
fprintf('\n===== EKF估计与理想姿态对比 =====\n');

% 计算误差
euler_error = euler_est - euler_ideal;

% 转换为角度制便于分析
euler_est_deg = rad2deg(euler_est);
euler_ideal_deg = rad2deg(euler_ideal);
euler_error_deg = rad2deg(euler_error);

% 处理角度环绕问题，将误差映射到[-180, 180]度范围内
euler_error_deg_wrapped = wrapTo180(euler_error_deg);

% 剔除误差大于180度的数据点（实际上通过wrapTo180已经处理，这里做额外检查）
valid_indices = true(1, length(time));
large_error_threshold = 180; % 定义大误差阈值

for i = 1:length(time)
    if any(abs(euler_error_deg_wrapped(:,i)) > large_error_threshold)
        valid_indices(i) = false;
    end
end

% 应用有效索引
time_valid = time(valid_indices);
euler_est_deg_valid = euler_est_deg(:, valid_indices);
euler_ideal_deg_valid = euler_ideal_deg(:, valid_indices);
euler_error_deg_valid = euler_error_deg_wrapped(:, valid_indices);

fprintf('原始数据点数: %d\n', length(time));
fprintf('有效数据点数: %d\n', sum(valid_indices));
fprintf('剔除数据点数: %d (%.2f%%)\n', sum(~valid_indices), sum(~valid_indices)/length(time)*100);

% 统计误差指标（使用有效数据）
rmse_roll = sqrt(mean(euler_error_deg_valid(1,:).^2));
rmse_pitch = sqrt(mean(euler_error_deg_valid(2,:).^2));
rmse_yaw = sqrt(mean(euler_error_deg_valid(3,:).^2));

mean_error_roll = mean(abs(euler_error_deg_valid(1,:)));
mean_error_pitch = mean(abs(euler_error_deg_valid(2,:)));
mean_error_yaw = mean(abs(euler_error_deg_valid(3,:)));

max_error_roll = max(abs(euler_error_deg_valid(1,:)));
max_error_pitch = max(abs(euler_error_deg_valid(2,:)));
max_error_yaw = max(abs(euler_error_deg_valid(3,:)));

std_error_roll = std(euler_error_deg_valid(1,:));
std_error_pitch = std(euler_error_deg_valid(2,:));
std_error_yaw = std(euler_error_deg_valid(3,:));

fprintf('滚转角误差: RMSE=%.4f°, 平均绝对误差=%.4f°, 最大误差=%.4f°, 标准差=%.4f°\n', ...
        rmse_roll, mean_error_roll, max_error_roll, std_error_roll);
fprintf('俯仰角误差: RMSE=%.4f°, 平均绝对误差=%.4f°, 最大误差=%.4f°, 标准差=%.4f°\n', ...
        rmse_pitch, mean_error_pitch, max_error_pitch, std_error_pitch);
fprintf('偏航角误差: RMSE=%.4f°, 平均绝对误差=%.4f°, 最大误差=%.4f°, 标准差=%.4f°\n', ...
        rmse_yaw, mean_error_yaw, max_error_yaw, std_error_yaw);

% 可视化对比
figure('Position', [100, 100, 1200, 800]);

% 滚转角对比
subplot(3,2,1);
plot(time_valid, euler_ideal_deg_valid(1,:), 'b-', 'LineWidth', 2, 'DisplayName', '理想值');
hold on;
plot(time_valid, euler_est_deg_valid(1,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'EKF估计');
title('滚转角对比 '); ylabel('角度 (deg)'); grid on;
legend('Location', 'best');

% 俯仰角对比
subplot(3,2,3);
plot(time_valid, euler_ideal_deg_valid(2,:), 'b-', 'LineWidth', 2, 'DisplayName', '理想值');
hold on;
plot(time_valid, euler_est_deg_valid(2,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'EKF估计');
title('俯仰角对比 '); ylabel('角度 (deg)'); grid on;
legend('Location', 'best');

% 偏航角对比
subplot(3,2,5);
plot(time_valid, euler_ideal_deg_valid(3,:), 'b-', 'LineWidth', 2, 'DisplayName', '理想值');
hold on;
plot(time_valid, euler_est_deg_valid(3,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'EKF估计');
title('偏航角对比 '); xlabel('时间 (s)'); ylabel('角度 (deg)'); grid on;
legend('Location', 'best');

% 误差分析
subplot(3,2,2);
plot(time_valid, euler_error_deg_valid(1,:), 'k-', 'LineWidth', 1.5);
title('滚转角误差 '); ylabel('误差 (deg)'); grid on;
ylim([-max_error_roll*1.1, max_error_roll*1.1]);

subplot(3,2,4);
plot(time_valid, euler_error_deg_valid(2,:), 'k-', 'LineWidth', 1.5);
title('俯仰角误差 '); ylabel('误差 (deg)'); grid on;
ylim([-max_error_pitch*1.1, max_error_pitch*1.1]);

subplot(3,2,6);
plot(time_valid, euler_error_deg_valid(3,:), 'k-', 'LineWidth', 1.5);
title('偏航角误差 '); xlabel('时间 (s)'); ylabel('误差 (deg)'); grid on;
ylim([-max_error_yaw*1.1, max_error_yaw*1.1]);

sgtitle('EKF姿态估计与理想值对比分析');

% 误差分布直方图
figure('Position', [100, 100, 1000, 600]);

subplot(2,3,1);
histogram(euler_error_deg_valid(1,:), 50, 'FaceColor', 'blue', 'EdgeColor', 'black');
title('滚转角误差分布'); xlabel('误差 (deg)'); ylabel('频数'); grid on;
xlim([-max_error_roll*1.1, max_error_roll*1.1]);

subplot(2,3,2);
histogram(euler_error_deg_valid(2,:), 50, 'FaceColor', 'green', 'EdgeColor', 'black');
title('俯仰角误差分布'); xlabel('误差 (deg)'); ylabel('频数'); grid on;
xlim([-max_error_pitch*1.1, max_error_pitch*1.1]);

subplot(2,3,3);
histogram(euler_error_deg_valid(3,:), 50, 'FaceColor', 'red', 'EdgeColor', 'black');
title('偏航角误差分布'); xlabel('误差 (deg)'); ylabel('频数'); grid on;
xlim([-max_error_yaw*1.1, max_error_yaw*1.1]);

% 累积分布函数
subplot(2,3,4);
ecdf(abs(euler_error_deg_valid(1,:)));
title('滚转角误差CDF'); xlabel('绝对误差 (deg)'); ylabel('累积概率'); grid on;

subplot(2,3,5);
ecdf(abs(euler_error_deg_valid(2,:)));
title('俯仰角误差CDF'); xlabel('绝对误差 (deg)'); ylabel('累积概率'); grid on;

subplot(2,3,6);
ecdf(abs(euler_error_deg_valid(3,:)));
title('偏航角误差CDF'); xlabel('绝对误差 (deg)'); ylabel('累积概率'); grid on;

sgtitle('姿态估计误差统计分析 (剔除大误差数据)');

% 显示被剔除的数据点（可选）
if any(~valid_indices)
    figure('Position', [100, 100, 800, 600]);
    
    removed_time = time(~valid_indices);
    removed_errors = euler_error_deg_wrapped(:, ~valid_indices);
    
    plot3(removed_errors(1,:), removed_errors(2,:), removed_errors(3,:), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    hold on;
    plot3(0, 0, 0, 'g+', 'MarkerSize', 15, 'LineWidth', 3);
    grid on;
    xlabel('滚转误差 (deg)'); ylabel('俯仰误差 (deg)'); zlabel('偏航误差 (deg)');
    title(sprintf('被剔除的数据点 (%d个)', sum(~valid_indices)));
    legend('被剔除点', '理想点', 'Location', 'best');
end

fprintf('\n对比分析完成！\n');

% 如果没有wrapTo180函数，使用以下替代函数
function wrapped_angles = wrapTo180(angles)
    % 将角度映射到[-180, 180]度范围
    wrapped_angles = mod(angles + 180, 360) - 180;
end