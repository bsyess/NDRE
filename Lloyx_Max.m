function [sigout] = Lloyx_Max(signalin, h, pnob)

%% 输入维度处理（保持与原始代码一致）
if size(signalin,1) > 1
    signalin = signalin.'; % 强制转为行向量
end
sigini = real(signalin);
%% 码表生成
Qbit_num = 2^pnob;
data = sigini(:); 

max_iter = 100;      % 最大迭代次数
tolerance = 1e-6;    % 容差
min_data = min(data);
max_data = max(data);
% 针对ACO信号（瑞利分布）的初始化
sigma = sqrt(mean(data.^2)/(2 - log(2))); % 估计瑞利参数
edges = raylinv(linspace(0,1,Qbit_num+1), sigma); % 瑞利分布分位数
centers = (edges(1:end-1) + edges(2:end)) / 2;

% Lloyd-Max迭代
prev_centers = inf(size(centers));
iter = 0;

while iter < max_iter && max(abs(centers - prev_centers)) > tolerance
    % 计算距离矩阵（隐式扩展）
    distances = abs(data - centers); % N×Qbit_num矩阵
    [~, cluster_idx] = min(distances, [], 2);
    
    % 更新码本
    new_centers = zeros(size(centers));
    for k = 1:Qbit_num
        members = data(cluster_idx == k);
        if ~isempty(members)
            new_centers(k) = mean(members);
        else
            new_centers(k) = centers(k); % 空簇保持原值
        end
    end
    
    % 更新状态
    prev_centers = centers;
    centers = new_centers;
    iter = iter + 1;
  
end
centers = sort(centers(:).'); % 确保码表为行向量

%% 动态步长计算
step_diff = diff(centers); % 各码本间距向量(1×Qbit_num-1)
edge_step = [step_diff(1), step_diff(end)]; % 首尾步长

%% 生成划分区间
ipartition = centers(1:end-1) + step_diff/2; % 相邻码表中点

%% 量化处理（保持原始逻辑）
[Iorder, iroundsig] = quantiz(sigini, ipartition, centers);
sigout=iroundsig;
if size(sigout,2) > 1
    sigout = sigout.';
end
end