function [J] = jacobian_euler_quat(q, axis)
% jacobian_euler_quat: 计算欧拉角（滚转/俯仰）对单位四元数的雅可比矩阵
% 输入：
%   q       - 单位四元数 [q0; q1; q2; q3]（列向量，对应姿态：q0实部，q1/q2/q3虚部分量）
%   axis    - 欧拉角类型：1=滚转角φ（绕载体Y轴），2=俯仰角θ（绕载体X轴）
% 输出：
%   J       - 雅可比矩阵（1×4，行向量，J(i) = ∂(欧拉角)/∂q(i)）
% 坐标系适配：载体坐标系（前X，右Y，下Z），导航坐标系（东北天），欧拉角顺序ZYX（偏航→俯仰→滚转）

% 第一步：提取四元数分量（确保输入为列向量）
q0 = q(1, 1);
q1 = q(2, 1);
q2 = q(3, 1);
q3 = q(4, 1);

% 第二步：根据axis选择欧拉角，推导雅可比（基于四元数转欧拉角公式的偏导）
switch axis
    case 1  % 计算滚转角φ对四元数的雅可比（φ = arctan2(2(q0q1+q2q3), 1-2(q1²+q2²))）
        % 定义分子Nφ = 2(q0q1 + q2q3)，分母Dφ = 1 - 2(q1² + q2²)
        N_phi = 2 * (q0*q1 + q2*q3);
        D_phi = 1 - 2 * (q1^2 + q2^2);
        % 雅可比公式：∂φ/∂q = (2/(Nφ²+Dφ²)) * [∂Nφ/∂q * Dφ - Nφ * ∂Dφ/∂q]
        % 因q为单位四元数，Nφ²+Dφ² = 1，简化系数为2
        coeff = 2;
        
        % 逐项计算∂φ/∂q0, ∂φ/∂q1, ∂φ/∂q2, ∂φ/∂q3
        dphi_dq0 = coeff * (q1 * D_phi - N_phi * 0);          % ∂Nφ/∂q0 = q1，∂Dφ/∂q0 = 0
        dphi_dq1 = coeff * (q0 * D_phi - N_phi * (-4*q1));    % ∂Nφ/∂q1 = q0，∂Dφ/∂q1 = -4q1
        dphi_dq2 = coeff * (q3 * D_phi - N_phi * (-4*q2));    % ∂Nφ/∂q2 = q3，∂Dφ/∂q2 = -4q2
        dphi_dq3 = coeff * (q2 * D_phi - N_phi * 0);          % ∂Nφ/∂q3 = q2，∂Dφ/∂q3 = 0
        
        % 组合滚转角雅可比
        J = [dphi_dq0, dphi_dq1, dphi_dq2, dphi_dq3];
        
    case 2  % 计算俯仰角θ对四元数的雅可比（θ = arcsin(2(q0q2 - q3q1))）
        % 定义正弦项Sθ = 2(q0q2 - q3q1)，则cosθ = sqrt(1 - Sθ²)（小角度下cosθ≈1，避免奇点）
        S_theta = 2 * (q0*q2 - q3*q1);
        cos_theta = sqrt(1 - S_theta^2);
        % 避免cosθ=0导致分母为0（俯仰角±90°时，车载场景中极少出现，此处用微小值替代）
        cos_theta = max(cos_theta, 1e-6);
        % 雅可比公式：∂θ/∂q = (2/cosθ) * ∂Sθ/∂q
        coeff = 2 / cos_theta;
        
        % 逐项计算∂θ/∂q0, ∂θ/∂q1, ∂θ/∂q2, ∂θ/∂q3
        dtheta_dq0 = coeff * q2;  % ∂Sθ/∂q0 = 2q2 → 乘以coeff后为2q2/cosθ
        dtheta_dq1 = coeff * (-q3);% ∂Sθ/∂q1 = -2q3 → 乘以coeff后为-2q3/cosθ
        dtheta_dq2 = coeff * q0;  % ∂Sθ/∂q2 = 2q0 → 乘以coeff后为2q0/cosθ
        dtheta_dq3 = coeff * (-q1);% ∂Sθ/∂q3 = -2q1 → 乘以coeff后为-2q1/cosθ
        
        % 组合俯仰角雅可比
        J = [dtheta_dq0, dtheta_dq1, dtheta_dq2, dtheta_dq3];
        
    otherwise
        error('axis参数错误！仅支持axis=1（滚转角）或axis=2（俯仰角）');
end

% 第三步：数值稳定性处理（避免四元数微小偏差导致的异常值）
J = max(min(J, 1e3), -1e3);  % 限制雅可比元素范围，防止数值发散
end