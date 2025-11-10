%% 辅助函数：四元数转欧拉角 (Z-Y-X顺序)
function euler = quat2euler(q)
    % 四元数格式：[q0; q1; q2; q3]
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    % 计算欧拉角（单位：rad）
    roll = atan2(2*(q0*q1 - q2*q3), 1 - 2*(q1^2 + q2^2));  
    pitch = asin(2*(-q0*q2 - q3*q1));                      
    yaw = atan2(2*(q0*q3 - q1*q2), 1 - 2*(q2^2 + q3^2));   
    
    
    euler = [roll; pitch; yaw];
end
