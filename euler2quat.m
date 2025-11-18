function q = euler2quat(roll, pitch, yaw)
    cr = cos(roll/2);
    sr = sin(roll/2);
    cp = cos(pitch/2);
    sp = sin(pitch/2);
    cy = cos(yaw/2);
    sy = sin(yaw/2);
    
    q0 = cr * cp * cy + sr * sp * sy;  % 标量部分
    q1 = sr * cp * cy - cr * sp * sy;  % x分量
    q2 = cr * sp * cy + sr * cp * sy;  % y分量
    q3 = cr * cp * sy - sr * sp * cy;  % z分量
    q = [q0; q1; q2; q3];  % 列向量输出
end