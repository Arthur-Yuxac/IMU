% 1. 加载数据与参数
% load('IMU_1015_1549_Recorder1_20251015165305.mat');  % 加载IMU原始数据
% data_generate;
generateData_43;
% load('IMU_1107_1549_Recorder1_20251107141152.mat');

% load('./IMU_param_calibrate/IMU_cali_tk/Accelerometer_Calib_Result.mat');
% load('./IMU_param_calibrate/IMU_cali_tk/Gyroscope_Calib_Result.mat');
% T_a = acc_calib_result.T_a;   K_a = acc_calib_result.K_a;   b_a = acc_calib_result.b_a;
% T_g = gyro_calib_result.T_g;   K_g = gyro_calib_result.K_g;   b_g = gyro_calib_result.b_g;

% 时间间隔计算（使用实际间隔而非均值）
N = length(time);
dt = zeros(1, N);
dt(1) = 0;  % 第一个时刻无间隔
for i = 2:N
    dt(i) = time(i) - time(i-1);
end
avg_freq = 1/mean(dt(2:end));  % 平均采样频率

% 2. IMU误差补偿
accel_compensated = zeros(3, N);
gyro_compensated = zeros(3, N);
% 预计算校准矩阵逆
calib_mat_a = T_a * K_a;
calib_mat_g = T_g * K_g;

for i = 1:N
    % 加速度计补偿：去零偏→刻度与轴间耦合校正
    accel_compensated(:,i) = calib_mat_a * (accel(:,i) - b_a);
    % 陀螺仪补偿：去零偏→刻度与轴间耦合校正
    gyro_compensated(:,i) = calib_mat_g * (gyro(:,i) - b_g)*pi/180;
end

% 3. UKF参数配置（优化初始状态与噪声参数）
% 初始状态优化：利用静止初始加速度估计姿态
x = [1;0;0;0; 0;0;0];  % [四元数; 陀螺零偏]
% if N >= 5  % 用前5个点平滑初始加速度
%     a0 = mean(accel_compensated(:,1:5), 2);
%     if norm(a0) > 0.5  % 确保有效重力测量
%         roll0 = atan2(a0(2), a0(3));          % 滚转角
%         pitch0 = atan2(-a0(1), norm(a0(2:3))); % 俯仰角
%         yaw0 = 0;                             % 偏航角初始为0
%         x(1:4) = euler2quat(roll0, pitch0, yaw0);  % 初始四元数
%     end
% end

% 协方差矩阵初始化（基于传感器手册参数）
P = eye(7);
P(1:4,1:4) = 1e-6 * eye(4);         % 四元数初始协方差（小角度误差）
P(5:7,5:7) = (10*pi/180/3600)^2 * eye(3);  % 陀螺零偏协方差（根据器件噪声）

% UKF参数
n = 7; % 状态维度
alpha = 1e-4; % 缩放参数
beta = 2;     % 分布参数
kappa = 0;    % 二阶缩放参数
lambda = alpha^2 * (n + kappa) - n;

% 过程噪声（动态调整框架）
Q_base = zeros(7);
gyro_noise_density = 0.01 * pi/180;  % 陀螺噪声密度 (rad/s/√Hz)
gyro_random_walk = 0.001 * pi/180;   % 陀螺随机游走 (rad/s²/√Hz)
Q_base(1:4,1:4) = (gyro_noise_density^2) * eye(4);
Q_base(5:7,5:7) = (gyro_random_walk^2) * eye(3);

% 测量噪声（加速度计，动态调整）
acc_noise_density = 0.01;     % 加速度计噪声密度 (m/s²/√Hz)
R_base = (acc_noise_density^2) * eye(3);

% 4. UKF主循环
q_est = zeros(4, N);     % 估计四元数
euler_est = zeros(3, N); % 估计欧拉角(rad)
bg_est = zeros(3, N);    % 估计陀螺零偏
g = 9.81;                % 重力加速度

% 定义四元数导数函数
quat_deriv = @(q, w) 0.5 * [
    0, -w(1), -w(2), -w(3);
    w(1), 0, w(3), -w(2);
    w(2), -w(3), 0, w(1);
    w(3), w(2), -w(1), 0
] * q;  % 四元数微分方程：q_dot = 0.5 * Omega(w) * q

% UKF权重计算
Wm = zeros(2*n+1, 1);
Wc = zeros(2*n+1, 1);
Wm(1) = lambda / (n + lambda);
Wc(1) = Wm(1) + (1 - alpha^2 + beta);
for i = 2:2*n+1
    Wm(i) = 1 / (2*(n + lambda));
    Wc(i) = Wm(i);
end

% 添加数值稳定性参数
epsilon = 1e-10;

for i = 1:N
    a = accel_compensated(:,i); 
    w = gyro_compensated(:,i);
    current_dt = dt(i);
    if current_dt <= 0 
        current_dt = mean(dt(2:end));
    end

    q = x(1:4);   % 当前四元数
    bg = x(5:7);  % 当前零偏

    % 陀螺温度补偿（同前）
    % temp = temp_imu(i);
    % temp_ref = 25;  
    % bg_temp_coeff = [1e-6; 1e-6; 1e-6];  % 需实际校准
    % bg_temp = bg_temp_coeff * (temp - temp_ref);
    % w_clean = w - bg - bg_temp;  % 去零偏+温度补偿后的角速度
    w_clean = w - bg;

    % UKF预测步骤 - 生成sigma点（添加数值稳定性）
    sigma_points = zeros(n, 2*n+1);
    
    % 确保协方差矩阵正定
    P_ensure_pos_def = P + epsilon * eye(n);
    [sqrt_P, flag] = chol((n + lambda) * P_ensure_pos_def, 'lower');
    if flag ~= 0
        % 如果Cholesky分解仍然失败，使用对角矩阵
        sqrt_P = sqrt((n + lambda)) * diag(sqrt(max(diag(P_ensure_pos_def), epsilon)));
    end
    
    sigma_points(:,1) = x;
    for j = 1:n
        sigma_points(:,j+1) = x + sqrt_P(:,j);
        sigma_points(:,j+n+1) = x - sqrt_P(:,j);
    end

    % 传播sigma点
    sigma_points_pred = zeros(n, 2*n+1);
    for j = 1:2*n+1
        q_sigma = sigma_points(1:4,j);
        bg_sigma = sigma_points(5:7,j);
        
        % 使用RK4积分传播四元数
        w_clean_sigma = w - bg_sigma;
        
        k1 = quat_deriv(q_sigma, w_clean_sigma); 
        k2 = quat_deriv(q_sigma + k1 * current_dt/2, w_clean_sigma);  
        k3 = quat_deriv(q_sigma + k2 * current_dt/2, w_clean_sigma); 
        k4 = quat_deriv(q_sigma + k3 * current_dt, w_clean_sigma); 
        
        q_pred_sigma = q_sigma + (k1 + 2*k2 + 2*k3 + k4) * current_dt / 6;
        q_pred_sigma = q_pred_sigma / norm(q_pred_sigma);
        
        % 零偏保持不变
        bg_pred_sigma = bg_sigma;
        
        sigma_points_pred(:,j) = [q_pred_sigma; bg_pred_sigma];
    end

    % 计算预测状态和协方差
    x_pred = zeros(n,1);
    for j = 1:2*n+1
        x_pred = x_pred + Wm(j) * sigma_points_pred(:,j);
    end
    
    % 归一化四元数部分
    x_pred(1:4) = x_pred(1:4) / norm(x_pred(1:4));
    
    P_pred = zeros(n,n);
    for j = 1:2*n+1
        dx = sigma_points_pred(:,j) - x_pred;
        % 对四元数部分进行特殊处理，确保在流形上
        if j > 1
            dq = quat_mult(sigma_points_pred(1:4,j), quat_conj(x_pred(1:4)));
            if dq(1) < 0  % 确保四元数符号一致
                dx(1:4) = -dx(1:4);
            end
        end
        P_pred = P_pred + Wc(j) * (dx * dx');
    end
    
    % 添加过程噪声并确保对称性
    Q = Q_base;
    Q(1:4,1:4) = Q(1:4,1:4) * current_dt;
    Q(5:7,5:7) = Q(5:7,5:7) * current_dt;
    P_pred = P_pred + Q;
    P_pred = 0.5 * (P_pred + P_pred');  % 确保对称性
    
    % 更新步骤
    a_norm = norm(a);
    g_norm = g;
    dynamic_ratio = abs(a_norm - g_norm) / g_norm;
    update_enable = (a_norm > 0.5) && (dynamic_ratio < 0.2);

    if update_enable
        % 生成新的sigma点用于更新（添加数值稳定性）
        P_pred_ensure_pos_def = P_pred + epsilon * eye(n);
        [sqrt_P_pred, flag] = chol((n + lambda) * P_pred_ensure_pos_def, 'lower');
        if flag ~= 0
            % 如果Cholesky分解失败，使用对角矩阵
            sqrt_P_pred = sqrt((n + lambda)) * diag(sqrt(max(diag(P_pred_ensure_pos_def), epsilon)));
        end
        
        sigma_points_pred_new = zeros(n, 2*n+1);
        sigma_points_pred_new(:,1) = x_pred;
        for j = 1:n
            sigma_points_pred_new(:,j+1) = x_pred + sqrt_P_pred(:,j);
            sigma_points_pred_new(:,j+n+1) = x_pred - sqrt_P_pred(:,j);
        end

        % 观测模型：预测的重力向量
        z_sigma = zeros(3, 2*n+1);
        for j = 1:2*n+1
            q_sigma = sigma_points_pred_new(1:4,j);
            C_pred = quat2rotmat(q_sigma);
            g_hat_body = C_pred' * [0; 0; -g];
            z_sigma(:,j) = g_hat_body / norm(g_hat_body);
        end

        % 计算预测观测值
        z_pred = zeros(3,1);
        for j = 1:2*n+1
            z_pred = z_pred + Wm(j) * z_sigma(:,j);
        end

        % 计算协方差
        Pzz = zeros(3,3);
        Pxz = zeros(n,3);
        for j = 1:2*n+1
            dz = z_sigma(:,j) - z_pred;
            dx = sigma_points_pred_new(:,j) - x_pred;
            % 对四元数部分进行特殊处理
            if j > 1
                dq = quat_mult(sigma_points_pred_new(1:4,j), quat_conj(x_pred(1:4)));
                if dq(1) < 0
                    dx(1:4) = -dx(1:4);
                end
            end
            
            Pzz = Pzz + Wc(j) * (dz * dz');
            Pxz = Pxz + Wc(j) * (dx * dz');
        end

        % 动态调整测量噪声R
        R = R_base * (1 + 10*dynamic_ratio);

        % 实际观测值（归一化加速度）
        z = a / a_norm;

        % 卡尔曼增益与状态更新
        S = Pzz + R;
        S = 0.5 * (S + S');  % 确保对称性
        
        % 检查S是否可逆
        [~, pos_def] = chol(S);
        if pos_def == 0
            K = Pxz / S;
            
            y = z - z_pred;
            x = x_pred + K * y;
            x(1:4) = x(1:4) / norm(x(1:4));  % 归一化四元数
            
            % 协方差更新（添加数值稳定性）
            P = P_pred - K * S * K';
            P = 0.5 * (P + P');  % 确保对称性
            P = P + epsilon * eye(n);  % 添加小的正定扰动
        else
            % 如果S不可逆，跳过更新
            x = x_pred;
            P = P_pred;
        end
    else
        x = x_pred;
        P = P_pred;
    end

    % 结果存储
    q_est(:,i) = x(1:4);
    euler_est(:,i) = quat2euler(x(1:4));
    bg_est(:,i) = x(5:7);
end

% 5. 结果输出与可视化
fprintf('\n===== UKF姿态解算结果 =====\n');
fprintf('平均采样频率: %.2f Hz | 数据点: %d\n', avg_freq, N);
fprintf('最终姿态角(deg): 滚转=%.2f, 俯仰=%.2f, 横摆=%.2f\n', ...
        rad2deg(euler_est(1,end)), rad2deg(euler_est(2,end)), rad2deg(euler_est(3,end)));

% 可视化
figure;
subplot(3,1,1); plot(time, rad2deg(euler_est(1,:))); title('滚转角(deg)'); grid on;
subplot(3,1,2); plot(time, rad2deg(euler_est(2,:))); title('俯仰角(deg)'); grid on;
subplot(3,1,3); plot(time, rad2deg(euler_est(3,:))); title('偏航角(deg)'); grid on;
sgtitle('UKF姿态估计');

save('euler_est_ukf.mat', 'euler_est');

% 辅助函数定义
function q_out = quat_mult(q, p)
    % 四元数乘法
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    p0 = p(1); p1 = p(2); p2 = p(3); p3 = p(4);
    
    q_out = [q0*p0 - q1*p1 - q2*p2 - q3*p3;
             q0*p1 + q1*p0 + q2*p3 - q3*p2;
             q0*p2 - q1*p3 + q2*p0 + q3*p1;
             q0*p3 + q1*p2 - q2*p1 + q3*p0];
end

function q_conj = quat_conj(q)
    % 四元数共轭
    q_conj = [q(1); -q(2); -q(3); -q(4)];
end