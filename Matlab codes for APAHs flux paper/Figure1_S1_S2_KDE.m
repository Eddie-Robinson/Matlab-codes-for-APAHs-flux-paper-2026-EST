% Fig.1 + Fig.S1 + Fig.S2
clearvars; clc;

% === MOD START: 读取新结果文件（变量名为 N×14 / N×16 / 名称列表） ===
load('Results_10000000.mat');   % 包含 flux_total_APAH14, flux_total_PAH16, APAH14_names, PAH16_names

set(groot, ...
    'defaultTextInterpreter',          'tex',  ...
    'defaultAxesTickLabelInterpreter', 'tex',  ...
    'defaultLegendInterpreter',        'tex',  ...
    'defaultAxesFontName',             'Arial',  ...
    'defaultTextFontName',             'Arial',  ...
    'defaultAxesFontWeight',           'bold',   ...
    'defaultTextFontWeight',           'bold'        ...
);

% 绘图前
medianColor = [1 0 0];       % 红色
meanColor   = [1 0.6 0.6];   % 淡红色
ci95Color   = [0, 0.45, 0.74];   % 蓝色
ci100Color  = [0.6 0.8 1];   % 淡蓝色

legendLabels = { ...
    '50% Median' , ...
    'Mean'         , ...
    '95% CI'    , ...
    % '100% Range'   ...
};

% START: 组装 Fig.1 的四类数据（顺序固定为 1)∑APAH8, 2)C1NAP, 3)∑C1–C4 NAPs, 4)∑8 PAHs） ===
APAH14_sum   = sum(flux_total_APAH14, 2);      % N×1
C1NAP_series = flux_total_APAH14(:,1);         % N×1 (APAH14 第1列)
C1_4NAP_sum  = sum(flux_total_APAH14(:,1:4),2);% N×1
PAH16_sum    = sum(flux_total_PAH16, 2);       % N×1

flux4 = [APAH14_sum, C1NAP_series, C1_4NAP_sum, PAH16_sum];  % N×4
PAH_types4 = {'\Sigma_{8} APAHs','MeNAP','\Sigma_{C1-C4} NAPs','\Sigma_{8} PAHs'};

% 图示1：KDE曲线（4类）+ 95%/Mean/Median
figure;
for k = 1:4
    subplot(1,4,k);
    data = flux4(:,k);
  
    % 只取正值部分，排除 ≤0 导致的 log 问题
    pos = data(data>0);
    if isempty(pos)
        title([PAH_types4{k} '（所有值 ≤ 0，无法绘制 KDE）']);
        continue;
    end

    % 1) KDE 曲线（对 log10(pos) 做 KDE）
    [f_log, xi_log] = ksdensity(log10(pos), 'NumPoints', 2000);
    xi = 10 .^ xi_log;
    hKDE = plot(xi, f_log, 'k', 'LineWidth', 1.5);
    hold on;

    % 2) 百分位
    pct = prctile(data, [0, 2.5, 50, 97.5, 100]);

    % 95% 区间
    xline(pct(2), '-', 'Color', ci95Color, 'LineWidth',1.2);
    xline(pct(4), '-', 'Color', ci95Color, 'LineWidth',1.2);

    % Mean / Median
    m = mean(data);
    xline(m, '-', 'Color', meanColor, 'LineWidth',1.2);
    xline(pct(3), '-', 'Color', medianColor, 'LineWidth',1.2);

    if k==1
        hMed   = plot(nan, nan, '-', 'Color', medianColor, 'LineWidth',1.2);
        hMean  = plot(nan, nan, '-', 'Color', meanColor,   'LineWidth',1.2);
        hCI95  = plot(nan, nan, '-', 'Color', ci95Color,   'LineWidth',1.2);
        legend( [hKDE, hMed, hMean, hCI95], ...
                {'KDE Curve', '50% Median', 'Mean', '95% CI'}, ...
                'Location','northeast', 'Interpreter','tex' );
    end

    pos = data(data>0);
    xmin = 10^floor(log10(min(pos)));
    xmax = 10^ceil(log10(max(pos)));
    set(gca, 'XScale', 'log', 'XLim', [xmin, xmax]);

    xlabel('Average annual emission flux (kt/year)');
    ylabel('Probability Density(dP/dlog_{10}F)');
    % 
    title(PAH_types4{k});

    grid on;
end

%% ===START: 图S1 — 14种 APAH 各自 KDE（中位数 = 0 的不画；2列排版） ===
figure;
valid_idx = [];  % 记录要画的物种索引
for k = 1:14
    data_k = flux_total_APAH14(:,k);
    if prctile(data_k,50) > 0   % 中位数>0 才纳入
        valid_idx(end+1) = k; %
    end
end
nvalid = numel(valid_idx);
nrows  = ceil(nvalid/2);
plot_count = 0;

for t = 1:nvalid
    k = valid_idx(t);
    data = flux_total_APAH14(:,k);
    pos  = data(data>0);
    if isempty(pos)
        continue; % 保险：避免全≤0 的情况
    end

    plot_count = plot_count + 1;
    subplot(nrows, 2, plot_count);

    [f_log, xi_log] = ksdensity(log10(pos), 'NumPoints', 2000);
    xi = 10 .^ xi_log;
    hKDE = plot(xi, f_log, 'k', 'LineWidth', 1.5); hold on;

    pct = prctile(data, [0, 2.5, 50, 97.5, 100]);
    xline(pct(2), '-', 'Color', ci95Color, 'LineWidth',1.2);
    xline(pct(4), '-', 'Color', ci95Color, 'LineWidth',1.2);
    m = mean(data);
    xline(m,        '-', 'Color', meanColor,   'LineWidth',1.2);
    xline(pct(3),   '-', 'Color', medianColor, 'LineWidth',1.2);

    if plot_count==1
        hMed   = plot(nan, nan, '-', 'Color', medianColor, 'LineWidth',1.2);
        hMean  = plot(nan, nan, '-', 'Color', meanColor,   'LineWidth',1.2);
        hCI95  = plot(nan, nan, '-', 'Color', ci95Color,   'LineWidth',1.2);
        legend( [hKDE, hMed, hMean, hCI95], ...
                {'KDE Curve', '50% Median', 'Mean', '95% CI'}, ...
                'Location','northeast', 'Interpreter','tex' );
    end

    xmin = 10^floor(log10(min(pos)));
    xmax = 10^ceil(log10(max(pos)));
    set(gca, 'XScale', 'log', 'XLim', [xmin, xmax]);

    xlabel('Average annual emission flux (kt/year)');
    ylabel('Probability Density(dP/dlog_{10}F)');
    title(APAH14_names{k});
    grid on;
end 

%% === START: 图S2 — 16种 PAH 各自 KDE（中位数 = 0 的不画；2列排版） ===
figure;
valid_idx = [];
for k = 1:16
    data_k = flux_total_PAH16(:,k);
    if prctile(data_k,50) > 0
        valid_idx(end+1) = k; %
    end
end
nvalid = numel(valid_idx);
nrows  = ceil(nvalid/2);
plot_count = 0;

for t = 1:nvalid
    k = valid_idx(t);
    data = flux_total_PAH16(:,k);
    pos  = data(data>0);
    if isempty(pos)
        continue;
    end

    plot_count = plot_count + 1;
    subplot(nrows, 2, plot_count);

    [f_log, xi_log] = ksdensity(log10(pos), 'NumPoints', 2000);
    xi = 10 .^ xi_log;
    hKDE = plot(xi, f_log, 'k', 'LineWidth', 1.5); hold on;

    pct = prctile(data, [0, 2.5, 50, 97.5, 100]);
    xline(pct(2), '-', 'Color', ci95Color,  'LineWidth',1.2);
    xline(pct(4), '-', 'Color', ci95Color,  'LineWidth',1.2);
    m = mean(data);
    xline(m,       '-', 'Color', meanColor,  'LineWidth',1.2);
    xline(pct(3),  '-', 'Color', medianColor,'LineWidth',1.2);

    if plot_count==1
        hMed   = plot(nan, nan, '-', 'Color', medianColor, 'LineWidth',1.2);
        hMean  = plot(nan, nan, '-', 'Color', meanColor,   'LineWidth',1.2);
        hCI95  = plot(nan, nan, '-', 'Color', ci95Color,   'LineWidth',1.2);
        legend( [hKDE, hMed, hMean, hCI95], ...
                {'KDE Curve', '50% Median', 'Mean', '95% CI'}, ...
                'Location','northeast', 'Interpreter','tex' );
    end

    xmin = 10^floor(log10(min(pos)));
    xmax = 10^ceil(log10(max(pos)));
    set(gca, 'XScale', 'log', 'XLim', [xmin, xmax]);

    xlabel('Average annual emission flux (kt/year)');
    ylabel('Probability Density(dP/dlog_{10}F)');
    title(PAH16_names{k});
    grid on;
end
%