function [targetMeas, numMeas] = CreatTargetMeas(frame, minDots, maxDots)

%% ��������������

% targetMeas Ŀ��������Ϣ ��ά���� ����*�������*Ŀ����*֡��
% numMeas ������� ��ά���� Ŀ����*֡��
MaxMeas = 30;
targetMeas = zeros(2, MaxMeas, 2, frame);
numMeas = zeros(2, frame); % ��Ŀ����������

%% ����Ŀ��һ����
% Բ��(200,200) �뾶150 Բ���˶�
angle1 = linspace(0,2*pi,frame);
targetx1 = 150 * cos(angle1) + 200;
targety1 = 150 * sin(angle1) + 200;

numMeas(1, :) = randi([minDots,maxDots],1,frame); % �������ÿ֡������
R = [16, 0; 0, 16];

% ��������ֵ
for frameIndex = 1:frame
    meas_noise = sample_gaussian(zeros(length(R), 1), R, numMeas(1, frameIndex))'; % �۲�����
    for meaIndex = 1:numMeas(1, frameIndex)
        targetMeas(:, meaIndex, 1, frameIndex) = [targetx1(frameIndex); targety1(frameIndex)] + meas_noise(:,meaIndex);
    end
end

%% ����Ŀ�������
% ������ʼ��
% ֡��������
T = 1; 
% ����ֱ���˶���ģ�͵�״̬ת�ƾ���
F = [1, T, 0, 0; ...
    0, 1, 0, 0; ...
    0, 0, 1, T; ...
    0, 0, 0, 1];
% ��������Э�������
qs = 0.1;
Q = [T^3 * qs / 3, T^2 * qs / 2, 0, 0; ...
    T^2 * qs / 2, T * qs, 0, 0; ...
    0, 0, T^3 * qs / 3, T^2 * qs / 2; ...
    0, 0, T^2 * qs / 2, T * qs];

%% ���ɹ۲������͹�������
numMeas(2, :) = randi([minDots,maxDots],1,frame); % �������ÿ֡������
process_noise2 = sample_gaussian(zeros(length(Q), 1), Q, frame)';

%% ����Ŀ��۲�켣
targetState = zeros(4, frame);
targetState(:, 1) = [0; 10; 0; 10];

for i = 2:frame
    targetState(:, i) = F * targetState(:, i-1) + process_noise2(:, i);
end

for frameIndex = 1:frame
    meas_noise = sample_gaussian(zeros(length(R), 1), R, numMeas(2, frameIndex))'; % �۲�����
    for meaIndex = 1:numMeas(2, frameIndex)
        targetMeas(:, meaIndex, 2, frameIndex) = [targetState(1, frameIndex); targetState(3, frameIndex)] + meas_noise(:,meaIndex);
    end
end
targetMeas = round(targetMeas);

%% ������Ƶ������
% vidObj = VideoWriter('Ŀ������.mp4', 'MPEG-4');
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
