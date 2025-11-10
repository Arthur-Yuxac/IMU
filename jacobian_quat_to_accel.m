% 四元数到加速度的雅可比矩阵
function J = jacobian_quat_to_accel(q, a_meas)
    J = zeros(3, 4);
    J(1, :) = [-a_meas(2)*q(3)-a_meas(3)*q(4), a_meas(2)*q(4)-a_meas(3)*q(3), ...
               -a_meas(1)*q(3)+a_meas(3)*q(2), a_meas(1)*q(4)+a_meas(2)*q(2)];
    J(2, :) = [a_meas(1)*q(3)-a_meas(3)*q(4), -a_meas(1)*q(4)-a_meas(3)*q(3), ...
               -a_meas(2)*q(4)+a_meas(3)*q(1), a_meas(2)*q(3)+a_meas(1)*q(1)];
    J(3, :) = [a_meas(1)*q(4)+a_meas(2)*q(3), -a_meas(1)*q(3)+a_meas(2)*q(4), ...
               -a_meas(1)*q(2)-a_meas(2)*q(1), a_meas(1)*q(1)-a_meas(2)*q(2)];
end