% 四元数更新函数
function q_new = quat_update(q, w, dt)
    Omega = [0, -w(1), -w(2), -w(3);
             w(1), 0, w(3), -w(2);
             w(2), -w(3), 0, w(1);
             w(3), w(2), -w(1), 0];
    q_new = q + 0.5 * Omega * q * dt;
    q_new = q_new / norm(q_new);
end

% % 1. 四元数更新（RK4积分）
% function q_new = quat_update(q, w, dt)
%     q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
%     wx = w(1); wy = w(2); wz = w(3);
%     
%     % RK4积分系数
%     k1 = 0.5 * dt * [ -q1*wx - q2*wy - q3*wz;
%                       q0*wx + q2*wz - q3*wy;
%                       q0*wy + q3*wx - q1*wz;
%                       q0*wz + q1*wy - q2*wx ];
%     
%     q_mid1 = q + 0.5*k1;
%     k2 = 0.5 * dt * [ -q_mid1(2)*wx - q_mid1(3)*wy - q_mid1(4)*wz;
%                       q_mid1(1)*wx + q_mid1(3)*wz - q_mid1(4)*wy;
%                       q_mid1(1)*wy + q_mid1(4)*wx - q_mid1(2)*wz;
%                       q_mid1(1)*wz + q_mid1(2)*wy - q_mid1(3)*wx ];
%     
%     q_mid2 = q + 0.5*k2;
%     k3 = 0.5 * dt * [ -q_mid2(2)*wx - q_mid2(3)*wy - q_mid2(4)*wz;
%                       q_mid2(1)*wx + q_mid2(3)*wz - q_mid2(4)*wy;
%                       q_mid2(1)*wy + q_mid2(4)*wx - q_mid2(2)*wz;
%                       q_mid2(1)*wz + q_mid2(2)*wy - q_mid2(3)*wx ];
%     
%     q_mid3 = q + k3;
%     k4 = 0.5 * dt * [ -q_mid3(2)*wx - q_mid3(3)*wy - q_mid3(4)*wz;
%                       q_mid3(1)*wx + q_mid3(3)*wz - q_mid3(4)*wy;
%                       q_mid3(1)*wy + q_mid3(4)*wx - q_mid3(2)*wz;
%                       q_mid3(1)*wz + q_mid3(2)*wy - q_mid3(3)*wx ];
%     
%     q_new = q + (k1 + 2*k2 + 2*k3 + k4)/6;
%     q_new = q_new / norm(q_new);  % 归一化
% end