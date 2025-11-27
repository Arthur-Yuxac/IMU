close all;clear;clc;
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


% 3. EKF参数配置（优化初始状态与噪声参数）
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

% 过程噪声（动态调整框架）
Q_base = zeros(7);
gyro_noise_density = 0.01 * pi/180;  % 陀螺噪声密度 (rad/s/√Hz)
gyro_random_walk = 0.001 * pi/180;   % 陀螺随机游走 (rad/s²/√Hz)
Q_base(1:4,1:4) = (gyro_noise_density^2) * eye(4);
Q_base(5:7,5:7) = (gyro_random_walk^2) * eye(3);

% 测量噪声（加速度计，动态调整）
acc_noise_density = 0.01;     % 加速度计噪声密度 (m/s²/√Hz)
R_base = (acc_noise_density^2) * eye(3);


% 4. EKF主循环
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

    k1 = quat_deriv(q, w_clean); 
    k2 = quat_deriv(q + k1 * current_dt/2, w_clean);  
    k3 = quat_deriv(q + k2 * current_dt/2, w_clean); 
    k4 = quat_deriv(q + k3 * current_dt, w_clean); 
    
    q_pred = q + (k1 + 2*k2 + 2*k3 + k4) * current_dt / 6;
    q_pred = q_pred / norm(q_pred); 

    % 零偏预测
    bg_pred = bg;

    % 状态转移矩阵F
    F = eye(7);
    Omega = [0, -w_clean(1), -w_clean(2), -w_clean(3);
             w_clean(1), 0, w_clean(3), -w_clean(2);
             w_clean(2), -w_clean(3), 0, w_clean(1);
             w_clean(3), w_clean(2), -w_clean(1), 0];
    F(1:4,1:4) = eye(4) + 0.5 * Omega * current_dt; 
    % 四元数对零偏的导数
    F(1:4,5:7) = -0.5 * [0, -current_dt, 0;
                         current_dt, 0, 0;
                         0, 0, -current_dt;
                         0, 0, current_dt];

    % 过程噪声Q
    Q = Q_base;
    Q(1:4,1:4) = Q(1:4,1:4) * current_dt;
    Q(5:7,5:7) = Q(5:7,5:7) * current_dt;

    % 预测协方差
    x_pred = [q_pred; bg_pred];
    P_pred = F * P * F' + Q;


    % 更新步骤
    a_norm = norm(a);
    g_norm = g;
    dynamic_ratio = abs(a_norm - g_norm) / g_norm;
    update_enable = (a_norm > 0.5) && (dynamic_ratio < 0.2);

    if update_enable
        % 计算预测重力向量
        C_pred = quat2rotmat(q_pred);
        g_hat_body = C_pred' * [0; 0; -g];
        
        % 观测残差
        y = (a / a_norm) - (g_hat_body / norm(g_hat_body));

        % 观测矩阵H
        q0 = q_pred(1); q1 = q_pred(2); q2 = q_pred(3); q3 = q_pred(4);
        Hq = [ -2*q2*g,   2*q3*g,  -2*q0*g,   2*q1*g;
                2*q1*g,   2*q0*g,   2*q3*g,   2*q2*g;
                2*q0*g,  -2*q1*g,  -2*q2*g,   2*q3*g ];
        g_hat_norm = norm(g_hat_body);
        if g_hat_norm > 0
            Hq = Hq / g_hat_norm;
        end
        Hg = zeros(3,3);
        H = [Hq, Hg];

        % 动态调整测量噪声R
        R = R_base * (1 + 10*dynamic_ratio);

        % 卡尔曼增益与状态更新
        S = H * P_pred * H' + R;
        K = P_pred * H' / S;
        x = x_pred + K * y;
        x(1:4) = x(1:4) / norm(x(1:4));  % 归一化
        P = (eye(7) - K * H) * P_pred;
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
fprintf('\n===== 优化后EKF姿态解算结果 =====\n');
fprintf('平均采样频率: %.2f Hz | 数据点: %d\n', avg_freq, N);
fprintf('最终姿态角(deg): 滚转=%.2f, 俯仰=%.2f, 横摆=%.2f\n', ...
        rad2deg(euler_est(1,end)), rad2deg(euler_est(2,end)), rad2deg(euler_est(3,end)));


% 可视化
figure;
subplot(3,1,1); plot(time, rad2deg(euler_est(1,:))); title('滚转角(deg)'); grid on;
subplot(3,1,2); plot(time, rad2deg(euler_est(2,:))); title('俯仰角(deg)'); grid on;
subplot(3,1,3); plot(time, rad2deg(euler_est(3,:))); title('偏航角(deg)'); grid on;
sgtitle('EKF姿态估计');

save('euler_est_ekf.mat', 'euler_est');