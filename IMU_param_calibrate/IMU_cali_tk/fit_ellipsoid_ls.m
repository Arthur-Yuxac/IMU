function [ellipsoid_params, v] = fit_ellipsoid_ls(stable_means, k)
% 功能：按文献方法拟合椭球（带kJ-I²约束）
% 输入：stable_means-稳定段均值（g），k-约束参数
% 输出：ellipsoid_params-椭球参数[a,b,c,f,g,h,p,q,r,d]，v-文献系数向量

M = size(stable_means, 1);
x = stable_means(:,1); y = stable_means(:,2); z = stable_means(:,3);  % 单位：g

% 1. 构造设计矩阵D（文献3.1节，X_i堆叠）
D = [x.^2, y.^2, z.^2, ...
     2*y.*z, 2*x.*z, 2*x.*y, ...
     2*x, 2*y, 2*z, ...
     ones(M,1)];

% 2. 构造约束矩阵C（文献公式7、8）
C1 = [-1,        k/2-1,    k/2-1, 0, 0, 0;
      k/2-1,    -1,        k/2-1, 0, 0, 0;
      k/2-1,    k/2-1,    -1,    0, 0, 0;
      0,        0,        0,    -k, 0, 0;
      0,        0,        0,     0, -k, 0;
      0,        0,        0,     0, 0, -k];
C = blkdiag(C1, zeros(4,4));  % 10×10约束矩阵

% 3. 计算DD^T与分块矩阵（文献公式9、11）
DDT = D' * D;
S11 = DDT(1:6, 1:6);  % 6×6子矩阵
S12 = DDT(1:6, 7:10); % 6×4子矩阵
S22 = DDT(7:10, 7:10);% 4×4子矩阵

% 4. 处理S22奇异（文献4.1节，用广义逆）
if rank(S22) < 4
    S22_inv = pinv(S22);
else
    S22_inv = inv(S22);
end

% 5. 求解广义特征值问题（文献公式14、15）
S_schur = S11 - S12 * S22_inv * S12';  % Schur补
C1_inv = inv(C1);
[V1, Lambda] = eig(C1_inv * S_schur);

% 6. 选择唯一正特征值对应向量（文献约束解）
eigenvalues = diag(Lambda);
positive_idx = find(eigenvalues > 1e-6);
if length(positive_idx) ~= 1
    [~, max_idx] = max(eigenvalues);
    positive_idx = max_idx;
end
v1 = V1(:, positive_idx);
v2 = -S22_inv * S12' * v1;

% 7. 归一化满足v^T*C*v=1（文献公式10）
v = [v1; v2];
constraint_val = v' * C * v;
v = v / sqrt(constraint_val);

% 输出椭球参数（对应文献方程1）
ellipsoid_params = v';
end
