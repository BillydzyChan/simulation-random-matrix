function [allMeas, targetMeas, numMeas] = CreatMeas(NFrames, range)

% allMeas ȫ��������Ϣ ��ά���� ����*�������*֡��
% targetMeas Ŀ��������Ϣ ��ά���� ����*�������*Ŀ����*֡��
% range ����ȡֵ��Χ[xMin, xMax, yMin, yMax]

%% ������ʼ��
% ����ȡֵ��Χ
xMin = range(1);
xMax = range(2);
yMin = range(3);
yMax = range(4);
x_area = xMax - xMin;
y_area = yMax - yMin;
% ���ݵ㼣�龯��ȡ�������ֵ
thresh = 4.1;
% ����˹�����ı�׼��
sigma = 1; % 0�൱��������
% �龯����
Pfa = 0.01; 
% �����
SNR = 15; 

%% ����Ŀ������
% ���������ķ�Χ
minMeas = 16;
maxMeas = 20;
[targetMeas, numMeas] = CreatTargetMeas(NFrames, minMeas, maxMeas);

% ���淶Χ��ÿһ֡�ز�����
data_len = ceil(Pfa*(x_area * 1.1)*(y_area * 1.1));
% ���Ⱦ���͵㼣�����ʼ��
amp_data = zeros(data_len, NFrames);
range_com = zeros(data_len, NFrames);
allMeas = zeros(2, data_len, NFrames) * nan;

%% ���ɼ������Ļز�����
for loop = 1:NFrames
    % �����ز��������ֲ���
    returnData = raylrnd(sqrt(sigma), x_area, y_area);
    % ����Ŀ��ز�
    for indexMeas = 1:numMeas(1, loop)
        xIndex = targetMeas(1, indexMeas, 1, loop) - xMin + 1; % ��Ӧ�ز����������
        yIndex = targetMeas(2, indexMeas, 1, loop) - yMin + 1;
        returnData(xIndex, yIndex) = raylrnd(sqrt(1+10^(SNR / 10))); % ���ȷ��������ֲ�
    end
    
    for indexMeas = 1:numMeas(2, loop)
        xIndex = targetMeas(1, indexMeas, 2, loop) - xMin + 1; %��Ӧ�ز����������
        yIndex = targetMeas(2, indexMeas, 2, loop) - yMin + 1;
        returnData(xIndex, yIndex) = raylrnd(sqrt(1+10^(SNR / 10)));
    end
    
    % �̶����޼��
    index = find(returnData > thresh); % ȡ�������޵㼣��ȫ�±�
    amp_data(1:length(index), loop) = returnData(index); % ȡ�������޵㼣�ķ���ֵ
    a = ceil(index/x_area); % ȡ�������޵㼣�����±꣬��������
    b = index - (a - 1) .* x_area; % ȡ�������޵㼣�����±꣬��������
    range_com(1:length(index), loop) = (b + xMin) + 1j * (a + yMin);
    % ȫ��������Ϣ
    allMeas(1,1:length(index),loop) = (b + xMin)';
    allMeas(2,1:length(index),loop) = (a + yMin)';
end

%% ������Ƶ������
% vidObj = VideoWriter('ȫ������.mp4', 'MPEG-4');
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

