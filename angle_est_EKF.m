% 1. 加载数据与参数
load('IMU_1015_1549_Recorder1_20251015165305.mat');  % 加载IMU原始数据
% load('IMU_Recorder1_0919_20250919171410.mat');
% load('IMU_1009_Recorder1_20251009121653.mat');
calib_file = 'IMU_params.mat';
load(calib_file);  % 加载IMU_params结构体



% 2. 数据提取与单位转换
% 时间戳：25us/bit → 秒，起点归零
t = CANFDMessage__AMS330_Timestamp(:,2)*25e-6;
t = t - t(1);
N = length(t);
dt = mean(diff(t));  % 平均采样间隔
accel_raw_x = CANFDMessage__acceleration_mg_0(:,2);
accel_raw_y = CANFDMessage__acceleration_mg_1(:,2);
accel_raw_z = CANFDMessage__acceleration_mg_2(:,2);
gyro_raw_x = CANFDMessage__angular_rate_mdps_0(:,2);
gyro_raw_y = CANFDMessage__angular_rate_mdps_1(:,2);
gyro_raw_z = CANFDMessage__angular_rate_mdps_2(:,2);
% 加速度：mg → m/s² (1mg=9.80665e-3 m/s²)
% accel_raw = [process_outliers(CANFDMessage__acceleration_mg_0(:,2));
%              process_outliers(CANFDMessage__acceleration_mg_1(:,2));
%              process_outliers(CANFDMessage__acceleration_mg_2(:,2))] * 9.80665e-3;
accel_raw = [accel_raw_x, accel_raw_y, accel_raw_z] * 9.80665e-3;
accel_raw = accel_raw';

% 陀螺仪：mdps → rad/s (1mdps=1e-3 deg/s = 1e-3×π/180 rad/s)
% gyro_raw = [process_outliers(CANFDMessage__angular_rate_mdps_0(:,2));
%             process_outliers(CANFDMessage__angular_rate_mdps_1(:,2));
%             process_outliers(CANFDMessage__angular_rate_mdps_2(:,2))] * 1e-3*pi/180;
gyro_raw = [gyro_raw_x, gyro_raw_y, gyro_raw_z] * 1e-3*pi/180;
gyro_raw = gyro_raw';

temp_imu = CANFDMessage__temperature_degC(:,2);  % 温度数据



% 3. 标定参数解包
accel_bias = IMU_params.accel.bias;       % 加速度计零偏 [3×1]
accel_scale_mat = IMU_params.accel.scale_mat;  % 加速度计刻度矩阵 [3×3]
gyro_bias = IMU_params.gyro.bias * 0.001;         % 陀螺仪零偏 [3×1]
gyro_scale_mat = IMU_params.gyro.scale_mat;    % 陀螺仪刻度矩阵 [3×3]


% 4. IMU误差补偿
accel_compensated = zeros(3,N);
gyro_compensated = zeros(3,N);
for i = 1:N
    % 加速度计：去零偏 → 刻度校正
    accel_compensated(:,i) = pinv(accel_scale_mat) * (accel_raw(:,i) - accel_bias);
    % 陀螺仪：去零偏 → 刻度校正
    gyro_compensated(:,i) = pinv(gyro_scale_mat) * (gyro_raw(:,i) - gyro_bias);
end


% 5. EKF参数配置（仅姿态解算：四元数+陀螺零偏，7维状态）
x = [1;0;0;0; 0;0;0];  % 初始状态：单位四元数+零偏为0
P = eye(7);
P(1:4,1:4) = 5e-7*eye(4);        % 四元数协方差
P(5:7,5:7) = (20*pi/180/3600)^2*eye(3);  % 陀螺零偏协方差
% P(1:4, 1:4)取经验值，P(5:7, 5:7)取gyro零偏不稳定性

Q = zeros(7);  % 过程噪声
Q(1:4,1:4) = (1.6568*pi/180)^2*dt^2*eye(4);  % 四元数噪声
Q(5:7,5:7) = (2.4*pi/180/60)^2*dt*eye(3);   % 零偏噪声
% Q(1:4, 1:4)取gyro方差，Q(5:7, 5:7)取gyro随机游走

R = (1.1456 * 9.80665 / 1000)^2*eye(3);  % 测量噪声（加速度计）


% 6. EKF主循环
q_est = zeros(4,N);     % 估计四元数
euler_est = zeros(3,N); % 估计欧拉角(rad)
bg_est = zeros(3,N);    % 估计陀螺零偏
g = 9.81;               % 重力加速度

for i = 1:N
    a = accel_compensated(:,i);  % 补偿后加速度
    w = gyro_compensated(:,i);   % 补偿后角速度
    
    q = x(1:4);   % 当前四元数
    bg = x(5:7);  % 当前零偏
    
    % 预测步骤
    w_clean = w - bg;  % 去零偏角速度
    % 四元数更新矩阵
    Omega = [0, -w_clean(1), -w_clean(2), -w_clean(3);
             w_clean(1), 0, w_clean(3), -w_clean(2);
             w_clean(2), -w_clean(3), 0, w_clean(1);
             w_clean(3), w_clean(2), -w_clean(1), 0];
    F = eye(7);
    F(1:4,1:4) = eye(4) + 0.5*Omega*dt;  % 状态转移
    
    x_pred = F*x;
    x_pred(1:4) = x_pred(1:4)/norm(x_pred(1:4));  % 四元数归一化
    P_pred = F*P*F' + Q;
    
    % 更新步骤（加速度有效时）
    if norm(a) > 0.05
        C_pred = quat2rotmat(x_pred(1:4));  % 旋转矩阵
        g_hat_body = C_pred'*[0;0;g];       % 载体系重力向量
        y = (a/norm(a)) - (g_hat_body/norm(g_hat_body));  % 残差
        
        % 观测矩阵
        H = zeros(3,7);
        H(:,1:4) = jacobian_quat_to_gravity(x_pred(1:4),g);  % 重力雅可比
        H(:,5:7) = jacobian_dy_dbg(x_pred(1:4),H(:,1:4),g,dt);  % 零偏雅可比
        
        K = P_pred*H'/(H*P_pred*H' + R);
        x = x_pred + K*y;
        x(1:4) = x(1:4)/norm(x(1:4));  % 归一化
        P = (eye(7)-K*H)*P_pred;
    else
        x = x_pred;
        P = P_pred;
    end
    
    % 结果存储
    q_est(:,i) = x(1:4);
    euler_est(:,i) = quat2euler(x(1:4));
    bg_est(:,i) = x(5:7);
end


% 7. 结果输出
fprintf('\n===== EKF姿态解算结果 =====\n');
fprintf('采样频率: %.2f Hz | 数据点: %d\n', 1/dt, N);
fprintf('最终姿态角(deg): 滚转=%.2f, 俯仰=%.2f, 偏航=%.2f\n', ...
        rad2deg(euler_est(1,end)), rad2deg(euler_est(2,end)), rad2deg(euler_est(3,end)));

% 保存结果
save('imu_attitude_results.mat', 't', 'q_est', 'euler_est', 'bg_est');
fprintf('✅ 结果保存至: imu_attitude_results.mat\n');

% 可视化
figure;
subplot(3,1,1); plot(t, rad2deg(euler_est(1,:))); title('滚转角(deg)'); grid on;
subplot(3,1,2); plot(t, rad2deg(euler_est(2,:))); title('俯仰角(deg)'); grid on;
subplot(3,1,3); plot(t, rad2deg(euler_est(3,:))); title('偏航角(deg)'); grid on;
sgtitle('EKF姿态估计');