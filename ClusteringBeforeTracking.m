function [halo] = ClusteringBeforeTracking(tarMeas)

%% �汾һ���ڸ���ǰ�Բ�����������з���ȥ���Ӳ���

%% �汾������ֱ�ӶԼ����������з���,֮�����þ������ƻ���

% %% ����������
% for i=1:6
%     tarMeas(:,i)=[12+randn*0.8;12+randn*0.8];
% end
% tarMeas(:,7)=[13;13];
% tarMeas(:,8)=[13.3;12.6];
% figure
% plot(tarMeas(1,:),tarMeas(2,:),'*');
% axis([0 20 0 20]);

%% ����������ĳɾ������ ��һ�к͵ڶ����Ǳ�ţ����������伸�ξ���
[tarMeasMat] = GetDistMat(tarMeas(:, tarMeas(1, :) ~= 0));
tarMeasMat(:, 3) = sqrt(tarMeasMat(:, 3));

%% ����ÿ������ܶȼ�������ܶȴ�Ĵ�֮�����̾���
ND = max(tarMeasMat(:, 2));
NL = max(tarMeasMat(:, 1));
if (NL > ND)
    ND = NL;
end
N = size(tarMeasMat, 1);
dist = zeros(ND, ND);
for i = 1:N
    ii = tarMeasMat(i, 1);
    jj = tarMeasMat(i, 2);
    dist(ii, jj) = tarMeasMat(i, 3);
    dist(jj, ii) = tarMeasMat(i, 3);
end
percent = 2.0;

position = round(N*percent/100);
if position == 0
    position = 1;
end
sda = sort(tarMeasMat(:, 3));
dc = sda(position);

rho = zeros(1, ND);

% Gaussian kernel
%
for i = 1:ND - 1
    for j = i + 1:ND
        rho(i) = rho(i) + exp(-(dist(i, j) / dc)*(dist(i, j) / dc));
        rho(j) = rho(j) + exp(-(dist(i, j) / dc)*(dist(i, j) / dc));
    end
end
%
% "Cut off" kernel
%
%for i=1:ND-1
%  for j=i+1:ND
%    if (dist(i,j)<dc)
%       rho(i)=rho(i)+1.;
%       rho(j)=rho(j)+1.;
%    end
%  end
%end

maxd = max(max(dist));

[~, ordrho] = sort(rho, 'descend');
delta(ordrho(1)) = -1.;
nneigh(ordrho(1)) = 0;

for ii = 2:ND
    delta(ordrho(ii)) = maxd;
    for jj = 1:ii - 1
        if (dist(ordrho(ii), ordrho(jj)) < delta(ordrho(ii)))
            delta(ordrho(ii)) = dist(ordrho(ii), ordrho(jj));
            nneigh(ordrho(ii)) = ordrho(jj);
        end
    end
end
delta(ordrho(1)) = max(delta(:));

% for i = 1:ND
%     ind(i) = i;
%     gamma(i) = rho(i) * delta(i);
% end
% figure
% tt = plot(rho(:), delta(:), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

%% ���߳��ж�����

%% �汾һ  Ĭ�������࣬����һ����Ŀ������ģ�һ�����Ӳ�
[delta_sorted, ~] = sort(delta, 'descend');
deltamin = delta_sorted(2);
NCLUST = 0;
cl = -1 * ones(1, ND);
for i = 1:ND
    if ((delta(i) >= deltamin))
        NCLUST = NCLUST + 1;
        cl(i) = NCLUST;
%         icl(NCLUST) = i;
    end
end

%% �汾�� ֱ�ӶԵõ���������з���
%��δ����

%% ������鵽������

for i = 1:ND
    if (cl(ordrho(i)) == -1)
        cl(ordrho(i)) = cl(nneigh(ordrho(i)));
    end
end
%halo
halo = cl;
if (NCLUST > 1)
    bord_rho = zeros(1, NCLUST);
    for i = 1:ND - 1
        for j = i + 1:ND
            if ((cl(i) ~= cl(j)) && (dist(i, j) <= dc))
                rho_aver = (rho(i) + rho(j)) / 2.;
                if (rho_aver > bord_rho(cl(i)))
                    bord_rho(cl(i)) = rho_aver;
                end
                if (rho_aver > bord_rho(cl(j)))
                    bord_rho(cl(j)) = rho_aver;
                end
            end
        end
    end
    for i = 1:ND
        if (rho(i) < bord_rho(cl(i)))
            halo(i) = 0;
        end
    end
end
for i = 1:NCLUST
    nc = 0;
    nh = 0;
    for j = 1:ND
        if (cl(j) == i)
            nc = nc + 1;
        end
        if (halo(j) == i)
            nh = nh + 1;
        end
    end
end
