
%% ʹ��Random Matrix����Ŀ����չ����
% ����
% �޸ķ���:
% 1. ǿ�����˸���ǰ������ȥ���ӵ� �б�Ҫ���� û��Ҫ��ȥ�� ��������ظ�
% 2. �������������Ϣɸѡ����
% 3. PDA �򵥵�JPDA���ɻ���NNCJPDA

clc;
clear;
close all;

%% =================������ʼ��================%%
% ��֡��
NFrames = 24;
% ������Ϣ�ṹ��
traceInfor = struct('state', [], 'P', [], 'extension', [], 'alpha', [], 'B', []);

T = 1; %֡��������
Rk = [0.5, 0; ...
    0, 0.5];%����״�޹ص��ǲ��ֵ�����Э�������
P1 = [Rk(1), Rk(1) / T, Rk(1) / (T^2); ...
    Rk(1) / (T^2), 2 * Rk(1) / (T^2), 3 * Rk(1) / (T^3); ...
    Rk(1) / (T^2), 3 * Rk(1) / (T^3), 5 * Rk(1) / (T^4)];% ��ʼЭ�������

% ״̬���̲�������
xita = 8; % �������ʱ�䳣�� 8
F = [1, T, 1 / 2 * T^2; ...
    0, 1, T; ...
    0, 0, exp(-T/xita)];% ״̬ת�ƾ���
xigema = 1.5; % ���ٶ�������׼��
Dk = xigema^2 * (1 - exp(-2*T/xita)) ...
    * [0, 0, 0; ...
    0, 0, 0; ...
    0, 0, 1];% ������˹��������wk��Э�������

% �۲ⷽ�̲�������
H = [1, 0, 0]; % �۲�ת�ƾ���

wid_size = 1; % ��Բ���Ű뾶

% ��չ�Զ�̬ģ�Ͳ���
deltak = 2; % ���ɶ� Ӧ�����㣺deltak > d-1
Ak = eye(2) / deltak; % ��״���¾��� ������С�������״��ȡ��������ʽ��
lamdak = 0.4; % ������״����Xk��������Ӱ��
d = 2; % ��ά�ռ�

% �ٶ�ֻ������Ŀ�꣬������ʼ�������ü����棩
traceInfor(1).state = [350; -5.5; 0; 200; 40; 0];
traceInfor(1).P = P1;
traceInfor(1).alpha = 20;
traceInfor(2).state = [0; 10; 0; 0; 10; 0];
traceInfor(2).P = P1;
traceInfor(2).alpha = 20;

% Ŀ���˶���Χ
xMin = -100;
xMax = 600;
yMin = -100;
yMax = 600;

ifelse = @(a, b, c)(a ~= 0) * b + (a == 0) * c;

%% =================�������================ %%
% targetMeas Ŀ��������Ϣ ��ά���� ����*�������*Ŀ����*֡��
[curMeas, targetMeas, numMeas] = CreatMeas(NFrames, [xMin, xMax, yMin, yMax]);

%% =================Ŀ����״�����ʼ��================ %%
% ����ÿ��Ŀ��
for iTrace = 1:length(traceInfor)
    % ����
    z1 = sum(targetMeas(:, :, iTrace, 1), 2) / numMeas(iTrace, 1);
    % ��״����
    traceInfor(iTrace).extension = zeros(2, 2);
    for indexMeas = 1:numMeas(iTrace, 1)
        traceInfor(iTrace).extension = traceInfor(iTrace).extension ...
            +(targetMeas(:, indexMeas, iTrace, 1) - z1) * (targetMeas(:, indexMeas, iTrace, 1) - z1)';
    end
    % Bk����������������ʵ֮���ʧ����� �ɶ�������Э����������õ� (����Bk���㲻֪����Դ)
    traceInfor(iTrace).B = sqrtm(lamdak*traceInfor(iTrace).extension+Rk) / sqrtm(traceInfor(iTrace).extension);
end


%% ��������Ͷ�Ӧ����Բ
TrackingFig = figure;
for loop = 1
    plot(curMeas(1, :, loop), curMeas(2, :, loop), '.', 'Color', [0.776, 0.776, 0.776]);
    hold on
    for iTrace = 1:length(traceInfor)
        plot(traceInfor(iTrace).state(1, 1), traceInfor(iTrace).state(4, 1), 'ro');
        % �õ���Բ���̲�����
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

%% ================���ٹ���============= %%
% ��Ƶ���1
vidObj = VideoWriter('Tracking by Random Matrix.mp4', 'MPEG-4');
vidObj.FrameRate = 4;
open(vidObj);
for loop = 2:NFrames
    disp(['No. ', num2str(loop)]);
    
    %% ��ǰ֡�������
    numberMeas = sum(~isnan(curMeas(1, :, loop)));
% %     % �������
% %     figure(TrackingFig);
% %     hold off
% %     plot(curMeas(1, :, loop), curMeas(2, :, loop), '.', 'Color', [0.776, 0.776, 0.776]);
% %     hold on
    %% ����Բ����ѡȡ��Ч�㼣
    % ���벨���ڵ�����
    curComMeas = zeros(2, 1);
    % ��Ϣ
    vv = zeros(2, length(traceInfor), numberMeas);
    for iTrace = 1:length(traceInfor)
        pPredict = F * traceInfor(iTrace).P(:, :, loop-1) * F' + Dk;
        SA = H * pPredict * H' * traceInfor(iTrace).extension(:, :, loop-1);
        % ͳ���ڲ����ڵ��������
        nk = 0;
        for indexMeas = 1:numberMeas
            vv(:, iTrace, indexMeas) = curMeas(:, indexMeas, loop) - ...
                kron(eye(2), H) * kron(eye(2), F) * traceInfor(iTrace).state(:, loop-1);
            if vv(:, iTrace, indexMeas)' / SA * vv(:, iTrace, indexMeas) < wid_size^2
                nk = nk + 1;
                curComMeas(:, nk, iTrace) = curMeas(:, indexMeas, loop);
            end
        end
%         % ����Բ����
%         center = kron(eye(2), H) * kron(eye(2), F) * traceInfor(iTrace).state(:, loop-1);
%         plot(center(1), center(2), 'b*');
%         syms x y
%         ellipse = [x - center(1); y - center(2)]'/ SA * [x - center(1); y - center(2)] - wid_size^2;
%         z = simplify(ellipse);
%         h = ezplot(z, [ifelse((xMin - 50) < (yMin - 50), (xMin - 50), (yMin - 50)), xMax + 50]);
%         set(h, 'Color', 'b');
%         % ����ɸѡ��������
%         plot(curComMeas(1,1:nk,iTrace), curComMeas(2,1:nk,iTrace), '.');
    end
    
    %% ����������ģ�����ֵ����������չ����
    zk = zeros(2, length(traceInfor));
    Zk = zeros(2, 2, length(traceInfor));
    for iTrace = 1:length(traceInfor)
        % ���û������Ļ�
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
    
    %% ״̬����״����
    for iTrace = 1:length(traceInfor)
        % ���ڸú������������
        if size(curComMeas,2) ~= 1
            nk = length(find(curComMeas(1, :, iTrace)));
        else
            nk = 0;
        end
        % ===== �˶����� ===== %
        % һ��Ԥ��
        statePredict = kron(eye(2), F) * traceInfor(iTrace).state(:, loop-1);
        % ״̬���Э����һ��Ԥ��
        pPredict = F * traceInfor(iTrace).P(:, :, loop-1) * F' + Dk;
        % ��Ϣ
        Gk = zk(:, iTrace) - kron(eye(2), H) * statePredict;
        % �������ӣ�������
        if nk ~= 0
            sPredict = H * pPredict * H' + det(traceInfor(iTrace).B)^(2 / d) / nk;
        else
            sPredict = H * pPredict * H';
        end
        % ����������
        K = pPredict * H' / sPredict;
        % Ŀ��״̬����
        if nk ~= 0
            traceInfor(iTrace).state(:, loop) = statePredict + kron(eye(2), K) * Gk;
        else
            traceInfor(iTrace).state(:, loop) = statePredict;
        end
        % Э����������
        traceInfor(iTrace).P(:, :, loop) = pPredict - K * sPredict * K';
        
        % ===== ��״���� ===== %
        % ���ɶȣ�������
        lamdak = traceInfor(iTrace).alpha(loop-1) - 2 * d - 2;
        % ���ɶȣ�������
        alphaPredict = 2 * deltak * (lamdak - 1) * (lamdak - 2) / lamdak / (lamdak + deltak) + 2 * d + 4;
        % ��״Ԥ��
        if nk ~= 0
            extensionPredict = deltak / lamdak * (alphaPredict - 2 * d - 2) ...
                * Ak * traceInfor(iTrace).extension(:, :, loop-1) * Ak';
        else
            extensionPredict = traceInfor(iTrace).extension(:, :, loop-1);
        end
        % ���¾���
        NPredict = sPredict \ Gk * Gk';
        % ��״����
        traceInfor(iTrace).extension(:, :, loop) = extensionPredict + NPredict ...
            +traceInfor(iTrace).B \ Zk(:, :, iTrace) * inv(traceInfor(iTrace).B)';
        % ���ɶȸ��£�������
        traceInfor(iTrace).alpha(loop) = alphaPredict + nk;
    end
    
    %% ��ͼ
    figure(TrackingFig);
%     hold off
    % �������
    h1 = plot(curMeas(1, :, loop), curMeas(2, :, loop), '.', 'Color', [0.776, 0.776, 0.776]);
    hold on    
    % ����ÿ������
    for iTrace = 1:length(traceInfor)      
        % ����ɸѡ��������
        nk = length(find(curComMeas(1, :, iTrace)));
        plot(curComMeas(1,1:nk,iTrace), curComMeas(2,1:nk,iTrace), '.');
        % ����ǰ״̬
        h2 = plot(traceInfor(iTrace).state(1, loop), traceInfor(iTrace).state(4, loop), 'o', 'Color', 'r');
        % ����Բ
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
    
    % ��Ƶ���2
    currFrame = getframe(gcf);
    writeVideo(vidObj, currFrame);
    
end
% ��Ƶ���3
close(vidObj);
