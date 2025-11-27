load('euler_ideal.mat');load('euler_est.mat');
% 6. 与理想姿态对比分析
fprintf('\n===== EKF估计与理想姿态对比 =====\n');

% 计算误差
euler_error = euler_est - euler_ideal;

% 转换为角度制便于分析
euler_est_deg = rad2deg(euler_est);
euler_ideal_deg = rad2deg(euler_ideal);
euler_error_deg = rad2deg(euler_error);

% 统计误差指标
rmse_roll = sqrt(mean(euler_error_deg(1,:).^2));
rmse_pitch = sqrt(mean(euler_error_deg(2,:).^2));
rmse_yaw = sqrt(mean(euler_error_deg(3,:).^2));

mean_error_roll = mean(abs(euler_error_deg(1,:)));
mean_error_pitch = mean(abs(euler_error_deg(2,:)));
mean_error_yaw = mean(abs(euler_error_deg(3,:)));

max_error_roll = max(abs(euler_error_deg(1,:)));
max_error_pitch = max(abs(euler_error_deg(2,:)));
max_error_yaw = max(abs(euler_error_deg(3,:)));

fprintf('滚转角误差: RMSE=%.4f°, 平均绝对误差=%.4f°, 最大误差=%.4f°\n', ...
        rmse_roll, mean_error_roll, max_error_roll);
fprintf('俯仰角误差: RMSE=%.4f°, 平均绝对误差=%.4f°, 最大误差=%.4f°\n', ...
        rmse_pitch, mean_error_pitch, max_error_pitch);
fprintf('偏航角误差: RMSE=%.4f°, 平均绝对误差=%.4f°, 最大误差=%.4f°\n', ...
        rmse_yaw, mean_error_yaw, max_error_yaw);

% 可视化对比
figure('Position', [100, 100, 1200, 800]);

% 滚转角对比
subplot(3,2,1);
plot(time, euler_ideal_deg(1,:), 'b-', 'LineWidth', 2, 'DisplayName', '理想值');
hold on;
plot(time, euler_est_deg(1,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'EKF估计');
title('滚转角对比'); ylabel('角度 (deg)'); grid on;
legend('Location', 'best');

% 俯仰角对比
subplot(3,2,3);
plot(time, euler_ideal_deg(2,:), 'b-', 'LineWidth', 2, 'DisplayName', '理想值');
hold on;
plot(time, euler_est_deg(2,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'EKF估计');
title('俯仰角对比'); ylabel('角度 (deg)'); grid on;
legend('Location', 'best');

% 偏航角对比
subplot(3,2,5);
plot(time, euler_ideal_deg(3,:), 'b-', 'LineWidth', 2, 'DisplayName', '理想值');
hold on;
plot(time, euler_est_deg(3,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'EKF估计');
title('偏航角对比'); xlabel('时间 (s)'); ylabel('角度 (deg)'); grid on;
legend('Location', 'best');

% 误差分析
subplot(3,2,2);
plot(time, euler_error_deg(1,:), 'k-', 'LineWidth', 1.5);
title('滚转角误差'); ylabel('误差 (deg)'); grid on;
ylim([-max_error_roll*1.1, max_error_roll*1.1]);

subplot(3,2,4);
plot(time, euler_error_deg(2,:), 'k-', 'LineWidth', 1.5);
title('俯仰角误差'); ylabel('误差 (deg)'); grid on;
ylim([-max_error_pitch*1.1, max_error_pitch*1.1]);

subplot(3,2,6);
plot(time, euler_error_deg(3,:), 'k-', 'LineWidth', 1.5);
title('偏航角误差'); xlabel('时间 (s)'); ylabel('误差 (deg)'); grid on;
ylim([-max_error_yaw*1.1, max_error_yaw*1.1]);

sgtitle('EKF姿态估计与理想值对比分析');

% 误差分布直方图
figure('Position', [100, 100, 1000, 600]);

subplot(2,3,1);
histogram(euler_error_deg(1,:), 50, 'FaceColor', 'blue', 'EdgeColor', 'black');
title('滚转角误差分布'); xlabel('误差 (deg)'); ylabel('频数'); grid on;

subplot(2,3,2);
histogram(euler_error_deg(2,:), 50, 'FaceColor', 'green', 'EdgeColor', 'black');
title('俯仰角误差分布'); xlabel('误差 (deg)'); ylabel('频数'); grid on;

subplot(2,3,3);
histogram(euler_error_deg(3,:), 50, 'FaceColor', 'red', 'EdgeColor', 'black');
title('偏航角误差分布'); xlabel('误差 (deg)'); ylabel('频数'); grid on;

% 累积分布函数
subplot(2,3,4);
ecdf(abs(euler_error_deg(1,:)));
title('滚转角误差CDF'); xlabel('绝对误差 (deg)'); ylabel('累积概率'); grid on;

subplot(2,3,5);
ecdf(abs(euler_error_deg(2,:)));
title('俯仰角误差CDF'); xlabel('绝对误差 (deg)'); ylabel('累积概率'); grid on;

subplot(2,3,6);
ecdf(abs(euler_error_deg(3,:)));
title('偏航角误差CDF'); xlabel('绝对误差 (deg)'); ylabel('累积概率'); grid on;

sgtitle('姿态估计误差统计分析');

fprintf('\n对比分析完成！\n');