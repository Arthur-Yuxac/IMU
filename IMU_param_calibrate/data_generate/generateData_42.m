% IMU绕中心点连续旋转模拟 - 理想加速度计数据生成
clear; clc;

% 参数设置
g = 9.81;           % 重力加速度 (m/s²)
f = 100;            % 采样频率 (Hz)
dt = 1/f;           % 采样间隔 (s)
n = 1100;           % 总采样点数
w = 2*pi;           % 转动角速率 (rad/s) - 1转/秒
static_points = 20; % 每次改变转动轴后的静态点数
extra_static_points = 10; % 额外增加的静态点数

% 初始方向（四元数表示）
q = [1 0 0 0]';     % 初始四元数 [w x y z]，表示无旋转

% 预分配加速度数据数组
acc_data = zeros(n, 3);
acc_static_data = []; % 静态点数据

% 生成第一个采样点的加速度（初始静止状态）
R0 = quat2rotm(q'); % 初始旋转矩阵
acc_data(1, :) = (R0 * [0; 0; -g])'; % 重力在IMU坐标系中的投影

% 生成均匀分布在球面上的转动轴
num_axes = ceil(n / 50); % 需要的转动轴数量
axes = zeros(3, num_axes);

% 使用球坐标生成均匀分布的转动轴
for i = 1:num_axes
    % 使用球坐标生成均匀分布的点
    theta = 2 * pi * rand();  % 方位角 [0, 2π]
    phi = acos(2 * rand() - 1); % 极角 [0, π]，使用反余弦确保均匀分布
    
    % 转换为笛卡尔坐标
    axes(1, i) = sin(phi) * cos(theta);
    axes(2, i) = sin(phi) * sin(theta);
    axes(3, i) = cos(phi);
end

% 状态标志
is_static_phase = false;
static_counter = 0;
axis_index = 1;
current_axis = axes(:, axis_index);

% 第一个静态点（初始状态）
acc_static_data = [acc_static_data; acc_data(1, :)];

% 生成连续的旋转运动
for i = 2:n
    % 检查是否需要改变转动轴
    if mod(i-1, 50) == 0 % 每50个采样点改变一次转动轴
        % 获取下一个均匀分布的转动轴
        axis_index = axis_index + 1;
        if axis_index > num_axes
            axis_index = 1; % 循环使用转动轴
        end
        current_axis = axes(:, axis_index);
        
        % 进入静态相位
        is_static_phase = true;
        static_counter = 0;
        
        fprintf('在采样点 %d 改变转动轴，开始静态相位\n', i);
    end
    
    % 处理静态相位
    if is_static_phase
        % 静态点：保持当前姿态不变
        acc_data(i, :) = acc_data(i-1, :);
        
        % 将静态点添加到静态数据数组
        acc_static_data = [acc_static_data; acc_data(i, :)];
        
        static_counter = static_counter + 1;
        
        % 检查静态相位是否结束
        if static_counter >= static_points
            is_static_phase = false;
            fprintf('在采样点 %d 结束静态相位，开始旋转\n', i);
        end
    else
        % 正常旋转相位
        % 计算当前时间步的旋转角度
        theta = w * dt;
        
        % 计算当前旋转的四元数
        q_rot = [cos(theta/2); 
                 current_axis * sin(theta/2)];
        
        % 更新IMU方向四元数
        q = quatmultiply(q', q_rot')';
        
        % 归一化四元数防止数值误差积累
        q = q / norm(q);
        
        % 计算当前时刻的旋转矩阵
        R = quat2rotm(q');
        
        % 计算加速度计读数（只有重力分量，因为绕中心旋转无线性加速度）
        acc_data(i, :) = (R * [0; 0; -g])';
    end
end

% 在末尾添加额外的10个静态点
last_acc = acc_data(end, :);
for i = 1:extra_static_points
    acc_data = [acc_data; last_acc];
    acc_static_data = [acc_static_data; last_acc];
end

% 更新总采样点数
n = n + extra_static_points;

% 辅助函数定义
function q = quatmultiply(q1, q2)
    % 四元数乘法
    w1 = q1(1); x1 = q1(2); y1 = q1(3); z1 = q1(4);
    w2 = q2(1); x2 = q2(2); y2 = q2(3); z2 = q2(4);
    
    q = [w1*w2 - x1*x2 - y1*y2 - z1*z2;
         w1*x2 + x1*w2 + y1*z2 - z1*y2;
         w1*y2 - x1*z2 + y1*w2 + z1*x2;
         w1*z2 + x1*y2 - y1*x2 + z1*w2]';
end

function R = quat2rotm(q)
    % 四元数转旋转矩阵
    w = q(1); x = q(2); y = q(3); z = q(4);
    
    R = [1-2*y^2-2*z^2,   2*x*y-2*z*w,   2*x*z+2*y*w;
         2*x*y+2*z*w,   1-2*x^2-2*z^2,   2*y*z-2*x*w;
         2*x*z-2*y*w,     2*y*z+2*x*w, 1-2*x^2-2*y^2];
end
