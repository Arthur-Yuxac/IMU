% 计算陀螺仪零偏(bg)对观测残差(y)的导数
function dy_dbg = jacobian_dy_dbg(q_pred, dG_dq ,g, dt)
% 输入参数：
% q_pred     -预测步骤后的四元数[4*1]
% dG_dq      -重力向量对四元数的导数[3*4]
% g          -重力加速度(m/s^2)
% dt         -采样时间间隔(s)
% 输出参数：
% dy_dbg     -残差对bg的导数，对应H矩阵的第11-13列[3*3]
C_pred = quat2rotmat(q_pred);
g_hat = C_pred * [0; 0; 1] * g;
g_hat_mag = norm(g_hat);
g_hat_norm = g_hat / g_hat_mag;

dnorm_dg = (eye(3) - g_hat_norm * g_hat_norm') / g_hat_mag;

q0 = q_pred(1); q1 = q_pred(2); q2 = q_pred(3); q3 = q_pred(4);
Omega_q = [0, -q1, -q2, -q3
           q1, 0, -q3, q2
           q2, q3, 0, -q1
           q3, -q2, q1, 0];
dq_dwclean = 0.5 * dt * Omega_q(:, 2:4);

dq_dbg = -dq_dwclean;
dg_hat_dbg = dG_dq * dq_dbg;
dy_dbg = -dnorm_dg * dg_hat_dbg;
end