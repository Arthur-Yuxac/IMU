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

% 3. 误差状态卡尔曼滤波(ESKF)参数配置
% 初始状态优化：利用静止初始加速度估计姿态
x = [1;0;0;0; 0;0;0];  % [四元数; 陀螺零偏]

% 误差状态：姿态误差(3x1) + 陀螺零偏误差(3x1)
dx = zeros(6,1);

% 协方差矩阵初始化（基于传感器手册参数）
P = eye(6);
P(1:3,1:3) = 1e-6 * eye(3);         % 姿态误差初始协方差
P(4:6,4:6) = (10*pi/180/3600)^2 * eye(3);  % 陀螺零偏误差协方差

% 过程噪声（动态调整框架）
Q_base = zeros(6);
gyro_noise_density = 0.01 * pi/180;  % 陀螺噪声密度 (rad/s/√Hz)
gyro_random_walk = 0.001 * pi/180;   % 陀螺随机游走 (rad/s²/√Hz)
Q_base(1:3,1:3) = (gyro_noise_density^2) * eye(3);
Q_base(4:6,4:6) = (gyro_random_walk^2) * eye(3);

% 测量噪声（加速度计，动态调整）
acc_noise_density = 0.01;     % 加速度计噪声密度 (m/s²/√Hz)
R_base = (acc_noise_density^2) * eye(3);

% 4. ESKF主循环
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
    
    w_clean = w - bg;

    % 名义状态预测（使用四阶龙格库塔法）
    k1 = quat_deriv(q, w_clean); 
    k2 = quat_deriv(q + k1 * current_dt/2, w_clean);  
    k3 = quat_deriv(q + k2 * current_dt/2, w_clean); 
    k4 = quat_deriv(q + k3 * current_dt, w_clean); 
    
    q_pred = q + (k1 + 2*k2 + 2*k3 + k4) * current_dt / 6;
    q_pred = q_pred / norm(q_pred); 

    % 零偏预测
    bg_pred = bg;

    % 误差状态转移矩阵F
    F = eye(6);
    % 从四元数到旋转矩阵
    C_nb = quat2rotmat(q_pred);
    
    % 姿态误差部分
    F(1:3,1:3) = eye(3) - skew_symmetric(w_clean) * current_dt;
    F(1:3,4:6) = -eye(3) * current_dt;
    
    % 零偏误差部分（随机游走模型）
    F(4:6,4:6) = eye(3);

    % 过程噪声Q
    Q = Q_base;
    Q(1:3,1:3) = Q(1:3,1:3) * current_dt;
    Q(4:6,4:6) = Q(4:6,4:6) * current_dt;

    % 预测误差状态协方差
    P_pred = F * P * F' + Q;

    % 更新步骤
    a_norm = norm(a);
    g_norm = g;
    dynamic_ratio = abs(a_norm - g_norm) / g_norm;
    update_enable = (a_norm > 0.5) && (dynamic_ratio < 0.2);

    if update_enable
        % 计算预测重力向量
        g_hat_body = C_nb' * [0; 0; -g];
        
        % 观测残差（在载体坐标系）
        y = (a / a_norm) - (g_hat_body / norm(g_hat_body));

        % 观测矩阵H（误差状态对观测的雅可比）
        H = zeros(3, 6);
        H(1:3, 1:3) = -skew_symmetric(g_hat_body);
        
        % 动态调整测量噪声R
        R = R_base * (1 + 10*dynamic_ratio);

        % 卡尔曼增益与误差状态更新
        S = H * P_pred * H' + R;
        K = P_pred * H' / S;
        dx = K * y;
        
        % 更新误差状态协方差
        P = (eye(6) - K * H) * P_pred;
        
        % 将误差状态注入到名义状态
        % 姿态误差注入（使用指数映射）
        dtheta = dx(1:3);
        dq = [1; 0.5*dtheta];
        dq = dq / norm(dq);
        q_updated = quat_multiply(q_pred, dq);
        q_updated = q_updated / norm(q_updated);
        
        % 零偏误差注入
        bg_updated = bg_pred + dx(4:6);
        
        % 更新名义状态
        x = [q_updated; bg_updated];
        
        % 重置误差状态
        dx = zeros(6,1);
    else
        % 如果没有更新，只使用预测值
        x = [q_pred; bg_pred];
        P = P_pred;
        dx = zeros(6,1);
    end

    % 结果存储
    q_est(:,i) = x(1:4);
    euler_est(:,i) = quat2euler(x(1:4));
    bg_est(:,i) = x(5:7);
end

% 5. 结果输出与可视化
fprintf('\n===== 误差卡尔曼滤波姿态解算结果 =====\n');
fprintf('平均采样频率: %.2f Hz | 数据点: %d\n', avg_freq, N);
fprintf('最终姿态角(deg): 滚转=%.2f, 俯仰=%.2f, 横摆=%.2f\n', ...
        rad2deg(euler_est(1,end)), rad2deg(euler_est(2,end)), rad2deg(euler_est(3,end)));

% 可视化
figure;
subplot(3,1,1); plot(time, rad2deg(euler_est(1,:))); title('滚转角(deg)'); grid on;
subplot(3,1,2); plot(time, rad2deg(euler_est(2,:))); title('俯仰角(deg)'); grid on;
subplot(3,1,3); plot(time, rad2deg(euler_est(3,:))); title('偏航角(deg)'); grid on;
sgtitle('ESKF姿态估计');

save('euler_est_eskf.mat', 'euler_est');

