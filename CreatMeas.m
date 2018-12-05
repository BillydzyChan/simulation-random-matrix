function [allMeas, targetMeas, numMeas] = CreatMeas(NFrames, range)

% allMeas 全部量测信息 三维矩阵 坐标*量测个数*帧数
% targetMeas 目标量测信息 四维矩阵 坐标*量测个数*目标数*帧数
% range 量测取值范围[xMin, xMax, yMin, yMax]

%% 参数初始化
% 量测取值范围
xMin = range(1);
xMax = range(2);
yMin = range(3);
yMax = range(4);
x_area = xMax - xMin;
y_area = yMax - yMin;
% 根据点迹虚警获取检测门限值
thresh = 4.1;
% 复高斯噪声的标准差
sigma = 1; % 0相当于无噪声
% 虚警概率
Pfa = 0.01; 
% 信噪比
SNR = 15; 

%% 生成目标量测
% 量测点个数的范围
minMeas = 16;
maxMeas = 20;
[targetMeas, numMeas] = CreatTargetMeas(NFrames, minMeas, maxMeas);

% 仿真范围和每一帧回波个数
data_len = ceil(Pfa*(x_area * 1.1)*(y_area * 1.1));
% 幅度矩阵和点迹矩阵初始化
amp_data = zeros(data_len, NFrames);
range_com = zeros(data_len, NFrames);
allMeas = zeros(2, data_len, NFrames) * nan;

%% 生成监测区域的回波幅度
for loop = 1:NFrames
    % 噪声回波（瑞利分布）
    returnData = raylrnd(sqrt(sigma), x_area, y_area);
    % 加入目标回波
    for indexMeas = 1:numMeas(1, loop)
        xIndex = targetMeas(1, indexMeas, 1, loop) - xMin + 1; % 对应回波矩阵的索引
        yIndex = targetMeas(2, indexMeas, 1, loop) - yMin + 1;
        returnData(xIndex, yIndex) = raylrnd(sqrt(1+10^(SNR / 10))); % 幅度服从瑞利分布
    end
    
    for indexMeas = 1:numMeas(2, loop)
        xIndex = targetMeas(1, indexMeas, 2, loop) - xMin + 1; %对应回波矩阵的索引
        yIndex = targetMeas(2, indexMeas, 2, loop) - yMin + 1;
        returnData(xIndex, yIndex) = raylrnd(sqrt(1+10^(SNR / 10)));
    end
    
    % 固定门限检测
    index = find(returnData > thresh); % 取出过门限点迹的全下标
    amp_data(1:length(index), loop) = returnData(index); % 取出过门限点迹的幅度值
    a = ceil(index/x_area); % 取出过门限点迹的列下标，即纵坐标
    b = index - (a - 1) .* x_area; % 取出过门限点迹的行下标，即横坐标
    range_com(1:length(index), loop) = (b + xMin) + 1j * (a + yMin);
    % 全部量测信息
    allMeas(1,1:length(index),loop) = (b + xMin)';
    allMeas(2,1:length(index),loop) = (a + yMin)';
end

%% 导出视频看量测
% vidObj = VideoWriter('全部量测.mp4', 'MPEG-4');
% vidObj.FrameRate = 4;
% open(vidObj);
% figure;
% for loop = 1:NFrames
%     hold off
%     plot(allMeas(1, :, loop), allMeas(2, :, loop), '.');
%     axis(range);
%     title(['No.', num2str(loop)]);
%     pause(0.1);
%     currFrame = getframe(gcf);
%     writeVideo(vidObj, currFrame);
% end
% close(vidObj);

