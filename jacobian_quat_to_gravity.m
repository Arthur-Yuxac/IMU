function dG_dq = jacobian_quat_to_gravity(q, g)
    % 输入：q-四元数[q0;q1;q2;q3]，g_vec-导航系重力向量[0;0;-g]
    % 输出：3×4雅可比矩阵（载体系重力向量对四元数的偏导）
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

    dG_dq = g * [
        2*q2,  2*q3,  2*q0, -2*q1;
       -2*q1,  2*q0, -2*q3,  2*q2;
        0,    -4*q1, -4*q2,   0
    ];
end