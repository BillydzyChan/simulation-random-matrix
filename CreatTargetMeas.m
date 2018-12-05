function [targetMeas, numMeas] = CreatTargetMeas(frame, minDots, maxDots)

%% 这个程序产生量测

% targetMeas 目标量测信息 四维矩阵 坐标*量测个数*目标数*帧数
% numMeas 量测个数 二维矩阵 目标数*帧数
MaxMeas = 30;
targetMeas = zeros(2, MaxMeas, 2, frame);
numMeas = zeros(2, frame); % 两目标的量测个数

%% 产生目标一量测
% 圆心(200,200) 半径150 圆周运动
angle1 = linspace(0,2*pi,frame);
targetx1 = 150 * cos(angle1) + 200;
targety1 = 150 * sin(angle1) + 200;

numMeas(1, :) = randi([minDots,maxDots],1,frame); % 随机产生每帧量测数
R = [16, 0; 0, 16];

% 生成量测值
for frameIndex = 1:frame
    meas_noise = sample_gaussian(zeros(length(R), 1), R, numMeas(1, frameIndex))'; % 观测噪声
    for meaIndex = 1:numMeas(1, frameIndex)
        targetMeas(:, meaIndex, 1, frameIndex) = [targetx1(frameIndex); targety1(frameIndex)] + meas_noise(:,meaIndex);
    end
end

%% 产生目标二量测
% 参数初始化
% 帧间采样间隔
T = 1; 
% 匀速直线运动的模型的状态转移矩阵
F = [1, T, 0, 0; ...
    0, 1, 0, 0; ...
    0, 0, 1, T; ...
    0, 0, 0, 1];
% 过程噪声协方差矩阵
qs = 0.1;
Q = [T^3 * qs / 3, T^2 * qs / 2, 0, 0; ...
    T^2 * qs / 2, T * qs, 0, 0; ...
    0, 0, T^3 * qs / 3, T^2 * qs / 2; ...
    0, 0, T^2 * qs / 2, T * qs];

%% 生成观测噪声和过程噪声
numMeas(2, :) = randi([minDots,maxDots],1,frame); % 随机产生每帧量测数
process_noise2 = sample_gaussian(zeros(length(Q), 1), Q, frame)';

%% 生成目标观测轨迹
targetState = zeros(4, frame);
targetState(:, 1) = [0; 10; 0; 10];

for i = 2:frame
    targetState(:, i) = F * targetState(:, i-1) + process_noise2(:, i);
end

for frameIndex = 1:frame
    meas_noise = sample_gaussian(zeros(length(R), 1), R, numMeas(2, frameIndex))'; % 观测噪声
    for meaIndex = 1:numMeas(2, frameIndex)
        targetMeas(:, meaIndex, 2, frameIndex) = [targetState(1, frameIndex); targetState(3, frameIndex)] + meas_noise(:,meaIndex);
    end
end
targetMeas = round(targetMeas);

%% 导出视频看量测
% vidObj = VideoWriter('目标量测.mp4', 'MPEG-4');
% vidObj.FrameRate = 4;
% open(vidObj);
% figure;
% for i = 1:size(targetMeas, 4)
%     hold off
%     plot(targetMeas(1, :, 1, i), targetMeas(2, :, 1, i), '.');
%     hold on
%     plot(targetMeas(1, :, 2, i), targetMeas(2, :, 2, i), '.');
%     axis([0, 400, 0, 400]);
%     title(['No.', num2str(i)]);
%     pause(0.1);
%     currFrame = getframe(gcf);
%     writeVideo(vidObj, currFrame);
% end
% close(vidObj);
