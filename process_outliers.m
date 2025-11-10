function processed = process_outliers(data)
    [dim, N] = size(data);
    processed = data;
    
    for i = 1:dim
        % 计算当前维度的均值和标准差
        mu = mean(data(i,:));
        sigma = std(data(i,:));
        upper_thresh = mu + 3*sigma;  % 上阈值
        lower_thresh = mu - 3*sigma;  % 下阈值
        
        % 查找异常值位置
        outliers = (data(i,:) < lower_thresh) | (data(i,:) > upper_thresh);
        outlier_idx = find(outliers);
        
        % 处理异常值
        for j = outlier_idx
            if j == 1  % 第一个数据点异常
                processed(i,j) = data(i,j+1);  % 用后一个数据替换
            elseif j == N  % 最后一个数据点异常
                processed(i,j) = data(i,j-1);  % 用前一个数据替换
            else  % 中间数据点异常
                processed(i,j) = (data(i,j-1) + data(i,j+1)) / 2;  % 用前后均值替换
            end
        end
        
        % 显示处理信息
        if ~isempty(outlier_idx)
            fprintf('维度 %d 处理了 %d 个异常值 (3σ范围: [%.4f, %.4f])\n', ...
                i, length(outlier_idx), lower_thresh, upper_thresh);
        end
    end
end