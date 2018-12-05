
%% 使用Random Matrix描述目标扩展特性
% 仿真
% 修改方向:
% 1. 强哥用了跟踪前聚类来去除杂点 有必要就留 没必要就去掉 避免跟他重复
% 2. 考虑引入幅度信息筛选量测
% 3. PDA 简单的JPDA即可或者NNCJPDA

clc;
clear;
close all;

%% =================参数初始化================%%
% 总帧数
NFrames = 24;
% 航迹信息结构体
traceInfor = struct('state', [], 'P', [], 'extension', [], 'alpha', [], 'B', []);

T = 1; %帧间采样间隔
Rk = [0.5, 0; ...
    0, 0.5];%与形状无关的那部分的噪声协方差矩阵
P1 = [Rk(1), Rk(1) / T, Rk(1) / (T^2); ...
    Rk(1) / (T^2), 2 * Rk(1) / (T^2), 3 * Rk(1) / (T^3); ...
    Rk(1) / (T^2), 3 * Rk(1) / (T^3), 5 * Rk(1) / (T^4)];% 初始协方差矩阵

% 状态方程参数设置
xita = 8; % 机动相关时间常数 8
F = [1, T, 1 / 2 * T^2; ...
    0, 1, T; ...
    0, 0, exp(-T/xita)];% 状态转移矩阵
xigema = 1.5; % 加速度噪声标准差
Dk = xigema^2 * (1 - exp(-2*T/xita)) ...
    * [0, 0, 0; ...
    0, 0, 0; ...
    0, 0, 1];% 独立高斯过程噪声wk的协方差矩阵

% 观测方程参数设置
H = [1, 0, 0]; % 观测转移矩阵

wid_size = 1; % 椭圆波门半径

% 扩展性动态模型参数
deltak = 2; % 自由度 应该满足：deltak > d-1
Ak = eye(2) / deltak; % 形状更新矩阵 描述大小或方向或形状（取决于其表达式）
lamdak = 0.4; % 描述形状矩阵Xk对噪声的影响
d = 2; % 二维空间

% 假定只有两个目标，航迹起始（后面用检测代替）
traceInfor(1).state = [350; -5.5; 0; 200; 40; 0];
traceInfor(1).P = P1;
traceInfor(1).alpha = 20;
traceInfor(2).state = [0; 10; 0; 0; 10; 0];
traceInfor(2).P = P1;
traceInfor(2).alpha = 20;

% 目标运动范围
xMin = -100;
xMax = 600;
yMin = -100;
yMax = 600;

ifelse = @(a, b, c)(a ~= 0) * b + (a == 0) * c;

%% =================量测产生================ %%
% targetMeas 目标量测信息 四维矩阵 坐标*量测个数*目标数*帧数
[curMeas, targetMeas, numMeas] = CreatMeas(NFrames, [xMin, xMax, yMin, yMax]);

%% =================目标形状矩阵初始化================ %%
% 对于每个目标
for iTrace = 1:length(traceInfor)
    % 形心
    z1 = sum(targetMeas(:, :, iTrace, 1), 2) / numMeas(iTrace, 1);
    % 形状矩阵
    traceInfor(iTrace).extension = zeros(2, 2);
    for indexMeas = 1:numMeas(iTrace, 1)
        traceInfor(iTrace).extension = traceInfor(iTrace).extension ...
            +(targetMeas(:, indexMeas, iTrace, 1) - z1) * (targetMeas(:, indexMeas, iTrace, 1) - z1)';
    end
    % Bk描述的是量测与真实之间的失真情况 由多个量测的协方差矩阵计算得到 (这里Bk计算不知道来源)
    traceInfor(iTrace).B = sqrtm(lamdak*traceInfor(iTrace).extension+Rk) / sqrtm(traceInfor(iTrace).extension);
end


%% 画出量测和对应的椭圆
TrackingFig = figure;
for loop = 1
    plot(curMeas(1, :, loop), curMeas(2, :, loop), '.', 'Color', [0.776, 0.776, 0.776]);
    hold on
    for iTrace = 1:length(traceInfor)
        plot(traceInfor(iTrace).state(1, 1), traceInfor(iTrace).state(4, 1), 'ro');
        % 得到椭圆方程并画出
        syms x y
        z = simplify([x - traceInfor(iTrace).state(1, 1), y - traceInfor(iTrace).state(4, 1)]/traceInfor(iTrace).extension(:, :, loop)*[x - traceInfor(iTrace).state(1, 1); y - traceInfor(iTrace).state(4, 1)]-1);
        h = ezplot(z, [ifelse((xMin - 50) < (yMin - 50), (xMin - 50), (yMin - 50)), xMax + 50]);
        set(h, 'Color', 'k');
    end
    
    axis([xMin, xMax, yMin, yMax]);
    title(['No.', num2str(loop)]);
    grid on
    pause(0.1);
end

%% ================跟踪过程============= %%
% 视频语句1
vidObj = VideoWriter('Tracking by Random Matrix.mp4', 'MPEG-4');
vidObj.FrameRate = 4;
open(vidObj);
for loop = 2:NFrames
    disp(['No. ', num2str(loop)]);
    
    %% 当前帧量测个数
    numberMeas = sum(~isnan(curMeas(1, :, loop)));
% %     % 画量测点
% %     figure(TrackingFig);
% %     hold off
% %     plot(curMeas(1, :, loop), curMeas(2, :, loop), '.', 'Color', [0.776, 0.776, 0.776]);
% %     hold on
    %% 画椭圆波门选取有效点迹
    % 落入波门内的量测
    curComMeas = zeros(2, 1);
    % 新息
    vv = zeros(2, length(traceInfor), numberMeas);
    for iTrace = 1:length(traceInfor)
        pPredict = F * traceInfor(iTrace).P(:, :, loop-1) * F' + Dk;
        SA = H * pPredict * H' * traceInfor(iTrace).extension(:, :, loop-1);
        % 统计在波门内的量测个数
        nk = 0;
        for indexMeas = 1:numberMeas
            vv(:, iTrace, indexMeas) = curMeas(:, indexMeas, loop) - ...
                kron(eye(2), H) * kron(eye(2), F) * traceInfor(iTrace).state(:, loop-1);
            if vv(:, iTrace, indexMeas)' / SA * vv(:, iTrace, indexMeas) < wid_size^2
                nk = nk + 1;
                curComMeas(:, nk, iTrace) = curMeas(:, indexMeas, loop);
            end
        end
%         % 画椭圆波门
%         center = kron(eye(2), H) * kron(eye(2), F) * traceInfor(iTrace).state(:, loop-1);
%         plot(center(1), center(2), 'b*');
%         syms x y
%         ellipse = [x - center(1); y - center(2)]'/ SA * [x - center(1); y - center(2)] - wid_size^2;
%         z = simplify(ellipse);
%         h = ezplot(z, [ifelse((xMin - 50) < (yMin - 50), (xMin - 50), (yMin - 50)), xMax + 50]);
%         set(h, 'Color', 'b');
%         % 画出筛选后的量测点
%         plot(curComMeas(1,1:nk,iTrace), curComMeas(2,1:nk,iTrace), '.');
    end
    
    %% 求量测的形心（即均值）和量测扩展特性
    zk = zeros(2, length(traceInfor));
    Zk = zeros(2, 2, length(traceInfor));
    for iTrace = 1:length(traceInfor)
        % 如果没有量测的话
        if size(curComMeas,2) ~= 1
            nk = length(find(curComMeas(1, :, iTrace)));
        else
            nk = 0;
        end
        if nk == 1
            zk(:, iTrace) = curComMeas(:, 1, iTrace);
        elseif nk >= 2
            zk(:, iTrace) = sum(curComMeas(:, 1:nk, iTrace), 2) / nk;
        elseif nk == 0
            zk(:, iTrace) = kron(eye(2), H) * kron(eye(2), F) * traceInfor(iTrace).state(:, loop-1);
        end
        if nk ~= 0
            for indexMeas = 1:nk
                Zk(:, :, iTrace) = Zk(:, :, iTrace) + ...
                    (curComMeas(:, indexMeas, iTrace) - zk(:, iTrace)) * (curComMeas(:, indexMeas, iTrace) - zk(:, iTrace))';
            end
        end
    end
    
    %% 状态和形状更新
    for iTrace = 1:length(traceInfor)
        % 属于该航迹的量测个数
        if size(curComMeas,2) ~= 1
            nk = length(find(curComMeas(1, :, iTrace)));
        else
            nk = 0;
        end
        % ===== 运动特性 ===== %
        % 一步预测
        statePredict = kron(eye(2), F) * traceInfor(iTrace).state(:, loop-1);
        % 状态误差协方差一步预测
        pPredict = F * traceInfor(iTrace).P(:, :, loop-1) * F' + Dk;
        % 新息
        Gk = zk(:, iTrace) - kron(eye(2), H) * statePredict;
        % 更新因子（标量）
        if nk ~= 0
            sPredict = H * pPredict * H' + det(traceInfor(iTrace).B)^(2 / d) / nk;
        else
            sPredict = H * pPredict * H';
        end
        % 卡尔曼增益
        K = pPredict * H' / sPredict;
        % 目标状态更新
        if nk ~= 0
            traceInfor(iTrace).state(:, loop) = statePredict + kron(eye(2), K) * Gk;
        else
            traceInfor(iTrace).state(:, loop) = statePredict;
        end
        % 协方差矩阵更新
        traceInfor(iTrace).P(:, :, loop) = pPredict - K * sPredict * K';
        
        % ===== 形状特性 ===== %
        % 自由度（标量）
        lamdak = traceInfor(iTrace).alpha(loop-1) - 2 * d - 2;
        % 自由度（标量）
        alphaPredict = 2 * deltak * (lamdak - 1) * (lamdak - 2) / lamdak / (lamdak + deltak) + 2 * d + 4;
        % 形状预测
        if nk ~= 0
            extensionPredict = deltak / lamdak * (alphaPredict - 2 * d - 2) ...
                * Ak * traceInfor(iTrace).extension(:, :, loop-1) * Ak';
        else
            extensionPredict = traceInfor(iTrace).extension(:, :, loop-1);
        end
        % 更新矩阵
        NPredict = sPredict \ Gk * Gk';
        % 形状更新
        traceInfor(iTrace).extension(:, :, loop) = extensionPredict + NPredict ...
            +traceInfor(iTrace).B \ Zk(:, :, iTrace) * inv(traceInfor(iTrace).B)';
        % 自由度更新（标量）
        traceInfor(iTrace).alpha(loop) = alphaPredict + nk;
    end
    
    %% 画图
    figure(TrackingFig);
%     hold off
    % 画量测点
    h1 = plot(curMeas(1, :, loop), curMeas(2, :, loop), '.', 'Color', [0.776, 0.776, 0.776]);
    hold on    
    % 对于每条航迹
    for iTrace = 1:length(traceInfor)      
        % 画出筛选后的量测点
        nk = length(find(curComMeas(1, :, iTrace)));
        plot(curComMeas(1,1:nk,iTrace), curComMeas(2,1:nk,iTrace), '.');
        % 画当前状态
        h2 = plot(traceInfor(iTrace).state(1, loop), traceInfor(iTrace).state(4, loop), 'o', 'Color', 'r');
        % 画椭圆
        syms x y
        ellipse = [x - traceInfor(iTrace).state(1, loop); y - traceInfor(iTrace).state(4, loop)]' ...
            / traceInfor(iTrace).extension(:, :, loop) ...
            * [x - traceInfor(iTrace).state(1, loop); y - traceInfor(iTrace).state(4, loop)] - 1;
        z = simplify(ellipse);
        h = ezplot(z, [ifelse((xMin - 50) < (yMin - 50), (xMin - 50), (yMin - 50)), xMax + 50]);
        set(h, 'Color', 'k');
    end
    axis([xMin, xMax, yMin, yMax]);
    grid on
    title(['No. ', num2str(loop)]);
    pause(0.1);
    
    % 视频语句2
    currFrame = getframe(gcf);
    writeVideo(vidObj, currFrame);
    
end
% 视频语句3
close(vidObj);
