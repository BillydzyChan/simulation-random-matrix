function [halo] = ClusteringBeforeTracking(tarMeas)

%% 版本一用于跟踪前对波门内量测进行分类去除杂波点

%% 版本二用于直接对检测后的量测进行分类,之后再用距离外推或者

% %% 测试用数据
% for i=1:6
%     tarMeas(:,i)=[12+randn*0.8;12+randn*0.8];
% end
% tarMeas(:,7)=[13;13];
% tarMeas(:,8)=[13.3;12.6];
% figure
% plot(tarMeas(1,:),tarMeas(2,:),'*');
% axis([0 20 0 20]);

%% 将输入量测改成距离矩阵 第一列和第二列是标号，第三列是其几何距离
[tarMeasMat] = GetDistMat(tarMeas(:, tarMeas(1, :) ~= 0));
tarMeasMat(:, 3) = sqrt(tarMeasMat(:, 3));

%% 计算每个点的密度及与比其密度大的簇之间的最短距离
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

%% 决策出有多少类

%% 版本一  默认有两类，其中一类是目标产生的，一类是杂波
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

%% 版本二 直接对得到的量测进行分类
%尚未开发

%% 将量测归到所在类

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
