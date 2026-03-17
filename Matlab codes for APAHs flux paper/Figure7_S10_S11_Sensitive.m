% Figure7_S10_S11_Sensitivity_with_MediansMeans_new.m
% 读取不同 N 的 Monte Carlo 结果：
% Fig.7  : MeNAP / global(%) 随 N 的敏感性（老样子）
% Fig.S10: 8个非0 APAH（含MeNAP）——通量(kt/year)随N敏感性（2×4）
% Fig.S11: 8个非0 PAH ——通量(kt/year)随N敏感性（2×4）

%% Figure 7 + Fig.S10 + Fig.S11
clearvars; clc;

%  全局参考值（kt/year）
MeNAP_GEMS_Glabal = 69.483;
% S10无全球分母；Fig7只用MeNAP


% 模拟规模列表
Ns   = [1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7];
numN = numel(Ns);

% 用于确定8个非0物种索引与名称
refFile = '';
for j = numN:-1:1
    f = sprintf('Results_%d.mat', Ns(j));
    if exist(f,'file') == 2
        refFile = f;
        break;
    end
end
if isempty(refFile)
    error('No Results_*.mat file found for Ns list. Please check file names.');
end

Sref = load(refFile, 'flux_total_APAH14','flux_total_PAH16','APAH14_names','PAH16_names');
APAH14_names = Sref.APAH14_names;
PAH16_names  = Sref.PAH16_names;

idx_APAH_nz = find(prctile(Sref.flux_total_APAH14,50,1) > 0);   % 期望 8
idx_PAH_nz  = find(prctile(Sref.flux_total_PAH16, 50,1) > 0);   % 期望 8

fprintf('\n[CHECK] Reference file used for non-zero screening: %s\n', refFile);

fprintf('\n[CHECK] Non-zero APAH14 (median>0), col indices:\n');
disp(idx_APAH_nz);
fprintf('[CHECK] Non-zero APAH14 names (order used in S10):\n');
disp(APAH14_names(idx_APAH_nz));

fprintf('\n[CHECK] Non-zero PAH16 (median>0), col indices:\n');
disp(idx_PAH_nz);
fprintf('[CHECK] Non-zero PAH16 names (order used in S11):\n');
disp(PAH16_names(idx_PAH_nz));
% ===  ===

nA = numel(idx_APAH_nz);   % 期望 8
nP = numel(idx_PAH_nz);    % 期望 8

%  预存数据（按 N 存 cell），避免重复 load
ratio_collection_MeNAP = cell(numN,1);   % Fig.7：N×1 (%)
flux_collection_APAH8  = cell(numN,1);   % S10：N×8 (kt/year)
flux_collection_PAH8   = cell(numN,1);   % S11：N×8 (kt/year)

for i = 1:numN
    fname = sprintf('Results_%d.mat', Ns(i));
    if exist(fname,'file') ~= 2
        fprintf('[WARN] Missing file: %s (skip this N)\n', fname);
        ratio_collection_MeNAP{i} = [];
        flux_collection_APAH8{i}  = [];
        flux_collection_PAH8{i}   = [];
        continue;
    end

    S = load(fname, 'flux_total_APAH14','flux_total_PAH16');

    % Fig.7：MeNAP (= APAH14 col1) 比值 (%)
    r = S.flux_total_APAH14(:,1) / MeNAP_GEMS_Glabal * 100;
    ratio_collection_MeNAP{i} = r;

    % Fig.S10：8个非0 APAH 的通量 (kt/year)
    flux_collection_APAH8{i} = S.flux_total_APAH14(:, idx_APAH_nz);

    % Fig.S11：8个非0 PAH 的通量 (kt/year)
    flux_collection_PAH8{i}  = S.flux_total_PAH16(:, idx_PAH_nz);
end
% 

% 颜色定义（保持原样）
ci95Color   = [0, 0.45, 0.74];
medianColor = [1, 0, 0];
meanColor   = [1, 0.6, 0.6];

%% =========================
% Fig.7：MeNAP ratio (%) vs N
%
figure; 
ax = subplot(1,1,1); hold(ax,'on');

% —— 1) Violin 分布 —— 
halfw = diff(log10(Ns(1:2))) / 3;  % 保持原样
for j = 1:numN
    y = ratio_collection_MeNAP{j};
    if isempty(y), continue; end

    % 为 log 绘图安全，仅取正值 
    y = y(y>0);
    if numel(y)<2, continue; end

    xi_log = linspace(log10(min(y)), log10(max(y)), 400);
    [f_log, ~] = ksdensity(log10(y), xi_log);
    xi = 10 .^ xi_log;
    f = f_log / max(f_log) * halfw;
    f([1,end]) = 0;

    logN = log10(Ns(j));
    fill(10.^(logN - f), xi, [0.2 0.2 0.8], ...
         'FaceAlpha',0.2,'EdgeColor','none', 'HandleVisibility','off');
    fill(10.^(logN + f), xi, [0.2 0.2 0.8], ...
         'FaceAlpha',0.2,'EdgeColor','none', 'HandleVisibility','off');
end

% —— 2) 计算各 N 点的 95%/Median/Mean —— 
low95   = nan(numN,1);
high95  = nan(numN,1);
medvals = nan(numN,1);
meanvals= nan(numN,1);
for j = 1:numN
    y = ratio_collection_MeNAP{j};
    if isempty(y), continue; end

    % 为 log 绘图安全，仅取正值
    y = y(y>0);
    if isempty(y), continue; end

    low95(j)   = prctile(y,2.5);
    high95(j)  = prctile(y,97.5);
    medvals(j) = median(y);
    meanvals(j)= mean(y);
end

h25   = plot(Ns, low95,  '-', 'Color',ci95Color,  'LineWidth',1.5);
h975  = plot(Ns, high95, '-', 'Color',ci95Color,  'LineWidth',1.5);
hMed  = plot(Ns, medvals,'-', 'Color',medianColor,'LineWidth',2);
hMean = plot(Ns, meanvals,'-','Color',meanColor,  'LineWidth',2);

set(ax, 'XScale','log','YScale','log');
xlim([3e3, 3e7]);
ylim([0.01, 100]);
xticks([1e4,1e5,1e6,1e7]);
xticklabels({'10^{4}','10^{5}','10^{6}','10^{7}'});
set(ax, 'TickLabelInterpreter','tex');

xlabel('Simulation size N','Interpreter','none');
ylabel('Ratio to Global (%)','Interpreter','none');
title('MeNAP', 'Interpreter','tex');
grid(ax,'on'); box(ax,'on');

legend([h25,hMed,hMean], ...
    {'95% CI','Median','Mean'}, ...
    'Location','southwest','Interpreter','none');


%% 
% Fig.S10：8个非0 APAH（含MeNAP）通量 vs N
% 2行×4列

figure(710);  % 用固定编号避免覆盖
set(gcf,'Position',[100 100 1200 600]);

fprintf('\n===== [S10] Flux sensitivity for APAH8 (kt/year) =====\n');
for s = 1:nA
    k_col = idx_APAH_nz(s);
    fprintf('[S10] subplot %d: APAH14 col %d = %s\n', s, k_col, APAH14_names{k_col});
end

for s = 1:nA
    ax = subplot(2,4,s); hold(ax,'on');

    % —— 1) Violin 分布 —— 
    for j = 1:numN
        M = flux_collection_APAH8{j};
        if isempty(M), continue; end

        y = M(:,s);
        y = y(y>0);
        if numel(y)<2, continue; end

        xi_log = linspace(log10(min(y)), log10(max(y)), 400);
        [f_log, ~] = ksdensity(log10(y), xi_log);
        xi = 10 .^ xi_log;
        f = f_log / max(f_log) * halfw;
        f([1,end]) = 0;

        logN = log10(Ns(j));
        fill(10.^(logN - f), xi, [0.2 0.2 0.8], ...
             'FaceAlpha',0.2,'EdgeColor','none', 'HandleVisibility','off');
        fill(10.^(logN + f), xi, [0.2 0.2 0.8], ...
             'FaceAlpha',0.2,'EdgeColor','none', 'HandleVisibility','off');
    end

    % —— 2) 95%/Median/Mean —— 
    low95   = nan(numN,1);
    high95  = nan(numN,1);
    medvals = nan(numN,1);
    meanvals= nan(numN,1);

    % 收集全体正值用于动态 ylim
    allpos = [];
    for j = 1:numN
        M = flux_collection_APAH8{j};
        if isempty(M), continue; end
        y = M(:,s); y = y(y>0);
        if isempty(y), continue; end

        allpos = [allpos; y]; %#ok<AGROW>
        low95(j)   = prctile(y,2.5);
        high95(j)  = prctile(y,97.5);
        medvals(j) = median(y);
        meanvals(j)= mean(y);
    end

    h25   = plot(Ns, low95,  '-', 'Color',ci95Color,  'LineWidth',1.5);
    h975  = plot(Ns, high95, '-', 'Color',ci95Color,  'LineWidth',1.5);
    hMed  = plot(Ns, medvals,'-', 'Color',medianColor,'LineWidth',2);
    hMean = plot(Ns, meanvals,'-','Color',meanColor,  'LineWidth',2);

    set(ax, 'XScale','log','YScale','log');
    xlim([3e3, 3e7]);
    xticks([1e4,1e5,1e6,1e7]);
    xticklabels({'10^{4}','10^{5}','10^{6}','10^{7}'});
    set(ax, 'TickLabelInterpreter','tex');

    % 动态 ylim（基于该物种所有N的正值） 
    if ~isempty(allpos)
        yl = 10^floor(log10(min(allpos)));
        yh = 10^ceil(log10(max(allpos)));
        ylim([yl, yh]);
    end

    xlabel('Simulation size N','Interpreter','none');
    ylabel('Average annual emission flux (kt/year)','Interpreter','none');
    title(APAH14_names{idx_APAH_nz(s)}, 'Interpreter','tex');
    grid(ax,'on'); box(ax,'on');

    % 只在第一个子图放图例
    if s==1
        legend([h25,hMed,hMean], ...
            {'95% CI','Median','Mean'}, ...
            'Location','southwest','Interpreter','none');
    end
end


%% 
% Fig.S11：8个非0 PAH 通量 vs N（2×4）

figure(711);  % 用固定编号避免覆盖
set(gcf,'Position',[100 100 1200 600]);

fprintf('\n===== [S11] Flux sensitivity for PAH8 (kt/year) =====\n');
for s = 1:nP
    k_col = idx_PAH_nz(s);
    fprintf('[S11] subplot %d: PAH16 col %d = %s\n', s, k_col, PAH16_names{k_col});
end

for s = 1:nP
    ax = subplot(2,4,s); hold(ax,'on');

    % —— 1) Violin 分布 —— 
    for j = 1:numN
        M = flux_collection_PAH8{j};
        if isempty(M), continue; end

        y = M(:,s);
        y = y(y>0);
        if numel(y)<2, continue; end

        xi_log = linspace(log10(min(y)), log10(max(y)), 400);
        [f_log, ~] = ksdensity(log10(y), xi_log);
        xi = 10 .^ xi_log;
        f = f_log / max(f_log) * halfw;
        f([1,end]) = 0;

        logN = log10(Ns(j));
        fill(10.^(logN - f), xi, [0.2 0.2 0.8], ...
             'FaceAlpha',0.2,'EdgeColor','none', 'HandleVisibility','off');
        fill(10.^(logN + f), xi, [0.2 0.2 0.8], ...
             'FaceAlpha',0.2,'EdgeColor','none', 'HandleVisibility','off');
    end

    % —— 2) 95%/Median/Mean —— 
    low95   = nan(numN,1);
    high95  = nan(numN,1);
    medvals = nan(numN,1);
    meanvals= nan(numN,1);

    % 动态 ylim 的 allpos
    allpos = [];
    for j = 1:numN
        M = flux_collection_PAH8{j};
        if isempty(M), continue; end
        y = M(:,s); y = y(y>0);
        if isempty(y), continue; end

        allpos = [allpos; y]; %
        low95(j)   = prctile(y,2.5);
        high95(j)  = prctile(y,97.5);
        medvals(j) = median(y);
        meanvals(j)= mean(y);
    end

    h25   = plot(Ns, low95,  '-', 'Color',ci95Color,  'LineWidth',1.5);
    h975  = plot(Ns, high95, '-', 'Color',ci95Color,  'LineWidth',1.5);
    hMed  = plot(Ns, medvals,'-', 'Color',medianColor,'LineWidth',2);
    hMean = plot(Ns, meanvals,'-','Color',meanColor,  'LineWidth',2);

    set(ax, 'XScale','log','YScale','log');
    xlim([3e3, 3e7]);
    xticks([1e4,1e5,1e6,1e7]);
    xticklabels({'10^{4}','10^{5}','10^{6}','10^{7}'});
    set(ax, 'TickLabelInterpreter','tex');

    %  动态 ylim
    if ~isempty(allpos)
        yl = 10^floor(log10(min(allpos)));
        yh = 10^ceil(log10(max(allpos)));
        ylim([yl, yh]);
    end

    xlabel('Simulation size N','Interpreter','none');
    ylabel('Average annual emission flux (kt/year)','Interpreter','none');
    title(PAH16_names{idx_PAH_nz(s)}, 'Interpreter','tex');
    grid(ax,'on'); box(ax,'on');

    if s==1
        legend([h25,hMed,hMean], ...
            {'95% CI','Median','Mean'}, ...
            'Location','southwest','Interpreter','none');
    end
end

toc;
