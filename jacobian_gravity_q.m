function dG_dq = jacobian_gravity_q(q, g)
    % 计算g_hat = C_bn * [0;0;1] * g 对q的偏导（3×4矩阵）
    C_bn = quat2rotmat(q);
    g_hat = C_bn * [0; 0; 1] * g;
    g_norm = norm(g_hat);
    
    % 先求C_bn对q的导数（3×4），再乘以[0;0;1]*g
    dC_dq = jacobian_CbN_q(q);  % 3×4矩阵
    dG_dq_raw = dC_dq * [0; 0; 1] * g;  % 3×4矩阵（d(g_hat)/dq）
    
    % 归一化后的重力向量导数（链式法则）
    dG_dq = (dG_dq_raw * g_norm - g_hat * (g_hat' * dG_dq_raw) / g_norm) / (g_norm^2);
end