function [sigout,sigout2] = VTB_DRE_FCM(signalin, h, pnob)

%% 输入维度处理（保持与原始代码一致）
if size(signalin,1) > 1
    signalin = signalin.'; % 强制转为行向量
end
sigini = real(signalin);

%% 码表生成
Qbit_num = 2^pnob;
data = sigini(:); 
options = [1.1, 300, 1e-8, 0];
[centers, ~] = fcm(data, Qbit_num, options);
centers = sort(centers(:).'); 

%% 动态步长计算
step_diff = diff(centers); % 各码本间距向量(1×Qbit_num-1)
edge_step = [step_diff(1), step_diff(end)]; % 首尾步长

%% 生成划分区间
ipartition = centers(1:end-1) + step_diff/2; % 相邻码表中点

%% 量化处理（保持原始逻辑）
[Iorder, iroundsig] = quantiz(sigini, ipartition, centers);
Iorder = Iorder + 1; % 调整索引偏移

%% 误差补偿（动态步长修正）
m = length(h);
nvtb = length(sigini);
I_error = iroundsig - sigini;

I_error = [I_error, I_error(1:m)];
Iorder = [Iorder, Iorder(1:m)];

Qro1 = zeros(3,m);
Qro2 = zeros(2,m);

for ii = m+1:nvtb+m
    current_order = Iorder(ii);
    
    if current_order > 1 && current_order < Qbit_num
        % 中间区域：使用相邻码本实际步长
        local_step = step_diff(current_order-1); % 索引对应
        Qro1(1,:) = I_error(ii-m+1:ii) - [zeros(1,m-1), local_step];
        Qro1(2,:) = I_error(ii-m+1:ii);
        Qro1(3,:) = I_error(ii-m+1:ii) + [zeros(1,m-1), local_step];
        
        Qsq = convn(Qro1, h, 'same');
        [~, idx] = min(sum(abs(Qsq),2));
        I_error(ii) = Qro1(idx, m);
    else
        % 边界区域：使用首尾步长
        use_step = edge_step(1 + (current_order==Qbit_num));
        adjust_sign = sign(current_order - (Qbit_num+1)/2);
        Qro2(1,:) = I_error(ii-m+1:ii) - [zeros(1,m-1), use_step]*adjust_sign;
        Qro2(2,:) = I_error(ii-m+1:ii);
        
        Qsq = convn(Qro2, h, 'same');
        [~, idx] = min(sum(abs(Qsq),2));
        I_error(ii) = Qro2(idx, m);
    end
end

%% 输出
I_error1 = [I_error(end-m+1:end), I_error(m+1:nvtb)];
sigout = sigini + I_error1;

if size(sigout,2) > 1
    sigout = sigout.';
end
end