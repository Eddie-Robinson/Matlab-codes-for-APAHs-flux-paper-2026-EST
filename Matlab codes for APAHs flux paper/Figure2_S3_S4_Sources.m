clearvars; clc;

tic;

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



%% 图示2：小提琴图 + 对数坐标
figure;

% === START: 读取并组装四类的数据矩阵（按 Seepage / Extraction / Transportation / Consumption / Total） ===
% need：flux_nature_APAH14, flux_extraction_APAH14, flux_transport_APAH14, flux_consumption_APAH14, flux_total_APAH14
%          flux_nature_PAH16,  flux_extraction_PAH16,  flux_transport_PAH16,  flux_consumption_PAH16,  flux_total_PAH16


% 1) ΣAPAH14
APAH14_nat   = sum(flux_nature_APAH14,     2);
APAH14_ext   = sum(flux_extraction_APAH14, 2);
APAH14_trans = sum(flux_transport_APAH14,  2);
APAH14_cons  = sum(flux_consumption_APAH14,2);
APAH14_total = sum(flux_total_APAH14,      2);
DM1 = [APAH14_nat, APAH14_ext, APAH14_trans, APAH14_cons, APAH14_total];

% 2) MeNAP （= APAH14 第1列）
MeNAP_nat   = flux_nature_APAH14(:,1);
MeNAP_ext   = flux_extraction_APAH14(:,1);
MeNAP_trans = flux_transport_APAH14(:,1);
MeNAP_cons  = flux_consumption_APAH14(:,1);
MeNAP_total = flux_total_APAH14(:,1);
DM2 = [MeNAP_nat, MeNAP_ext, MeNAP_trans, MeNAP_cons, MeNAP_total];

% 3) ΣC1–C4 NAPs （= APAH14 前4列求和）
C1_4NAP_nat   = sum(flux_nature_APAH14(:,1:4),     2);
C1_4NAP_ext   = sum(flux_extraction_APAH14(:,1:4), 2);
C1_4NAP_trans = sum(flux_transport_APAH14(:,1:4),  2);
C1_4NAP_cons  = sum(flux_consumption_APAH14(:,1:4),2);
C1_4NAP_total = sum(flux_total_APAH14(:,1:4),      2);
DM3 = [C1_4NAP_nat, C1_4NAP_ext, C1_4NAP_trans, C1_4NAP_cons, C1_4NAP_total];

% 4) Σ16 PAHs
PAH16_nat   = sum(flux_nature_PAH16,     2);
PAH16_ext   = sum(flux_extraction_PAH16, 2);
PAH16_trans = sum(flux_transport_PAH16,  2);
PAH16_cons  = sum(flux_consumption_PAH16,2);
PAH16_total = sum(flux_total_PAH16,      2);
DM4 = [PAH16_nat, PAH16_ext, PAH16_trans, PAH16_cons, PAH16_total];

DMs = {DM1, DM2, DM3, DM4};
PAH_types4 = {'\Sigma_{14} APAHs','MeNAP','\Sigma_{C1-C4} NAPs','\Sigma_{16} PAHs'};
% 

for k = 1:4
    subplot(1,4,k);

    data_matrix = DMs{k};  % N×5
    
    % 只保留正值
    data_pos = data_matrix(data_matrix>0);
    if isempty(data_pos)
        title([PAH_types4{k} '（所有值≤0，无法绘制）']);
        continue;
    end
    
    % 动态 log 轴范围
    lo = 10^floor(log10(min(data_pos)));
    hi = 10^ceil(log10(max(data_pos)));
    set(gca,'YScale','log','YLim',[lo,hi]);
    hold on;
    
    labels = {'Seepage','Extraction','Transportation','Consumption','Total'};
    violin_halfwidth = 0.3;   % 整个 Violin 半宽最大值
    rim_width = 0.12;         % 2.5%/97.5% 线的固定半宽度
    
    for s = 1:5
        y = data_matrix(:,s);
        y = y(y>0);
        if numel(y)<2, continue; end
        
        % 在 log 空间等距采样
        xi_log = linspace(log10(min(y)), log10(max(y)), 2000);
        [f_log, ~] = ksdensity(log10(y), xi_log);
        xi = 10 .^ xi_log;
        f = f_log / max(f_log) * violin_halfwidth;
        
        % 画小提琴
        hPatch = fill([s - f, fliplr(s + f)], [xi, fliplr(xi)], ...
                      [0.2 0.2 0.8], 'FaceAlpha',0.2, 'EdgeColor','none');
        if k==1 && s==1
            hViolin = hPatch;  % “小提琴图”
        end

        % 分位点
        q = prctile(y, [0, 2.5, 50, 97.5, 100]);
        
        % 中位数处的半宽
        [~, idx_med] = min(abs(xi - q(3)));
        med_hw = f(idx_med);

        % —— 95% 区间 (蓝色实线) ——
        line([s-med_hw, s+med_hw], [q(2), q(2)], ...
            'Color', ci95Color,  'LineWidth',1.2);
        line([s-med_hw, s+med_hw], [q(4), q(4)], ...
            'Color', ci95Color,  'LineWidth',1.2);

        % —— Mean (淡红色实线) ——
        m = mean(y);
        line([s-med_hw, s+med_hw], [m, m], ...
            'Color', meanColor, 'LineWidth',1.5);

        % —— 50% 中位数 (红色实线) ——
        line([s-med_hw, s+med_hw], [q(3), q(3)], ...
            'Color', medianColor,'LineWidth',1.5);
    
    end

    % —— 只在第一个子图一次性创建图例 —— 
    if k==1
        hMed   = plot(nan, nan, '-', 'Color', medianColor, 'LineWidth',1.2);
        hMean  = plot(nan, nan, '-', 'Color', meanColor,   'LineWidth',1.2);
        hCI95  = plot(nan, nan, '-', 'Color', ci95Color,   'LineWidth',1.2);
        legend( [hViolin, hMed, hMean, hCI95], ...
                {'Density','50% Median','Mean','95% CI'}, ...
                'Location','northeast', 'Interpreter','tex' );
    end

    xticks(1:5); 
    xticklabels(labels);
    xlim([0.5, 5.5]);     
    grid on; box on;               
    
    %
    title([PAH_types4{k}]);
    
    ylabel('Average annual emission flux (kt/year)');
end

%% === START: 图S3 — 14种 APAH 的小提琴图（中位数=0 的不画；2列排版） ===
figure;
valid_idx = [];
for k = 1:14
    if prctile(flux_total_APAH14(:,k),50) > 0
        valid_idx(end+1) = k; %
    end
end
nvalid = numel(valid_idx);
nrows  = ceil(nvalid/4);
plot_count = 0;

for t = 1:nvalid
    k = valid_idx(t);
    data_matrix = [ ...
        flux_nature_APAH14(:,k), ...
        flux_extraction_APAH14(:,k), ...
        flux_transport_APAH14(:,k), ...
        flux_consumption_APAH14(:,k), ...
        flux_total_APAH14(:,k) ];

    data_pos = data_matrix(data_matrix>0);
    if isempty(data_pos)
        continue;
    end

    plot_count = plot_count + 1;
    subplot(nrows,4,plot_count);

    lo = 10^floor(log10(min(data_pos)));
    hi = 10^ceil(log10(max(data_pos)));
    set(gca,'YScale','log','YLim',[lo,hi]);
    hold on;

    labels = {'Seepage','Extraction','Transportation','Consumption','Total'};
    violin_halfwidth = 0.3;
    rim_width = 0.12;

    for s = 1:5
        y = data_matrix(:,s); y = y(y>0);
        if numel(y)<2, continue; end

        xi_log = linspace(log10(min(y)), log10(max(y)), 2000);
        [f_log, ~] = ksdensity(log10(y), xi_log);
        xi = 10 .^ xi_log;
        f = f_log / max(f_log) * violin_halfwidth;

        hPatch = fill([s - f, fliplr(s + f)], [xi, fliplr(xi)], ...
                      [0.2 0.2 0.8], 'FaceAlpha',0.2, 'EdgeColor','none');
        if plot_count==1 && s==1
            hViolin = hPatch;
        end

        q = prctile(y, [0, 2.5, 50, 97.5, 100]);
        [~, idx_med] = min(abs(xi - q(3)));
        med_hw = f(idx_med);

        line([s-med_hw, s+med_hw], [q(2), q(2)], 'Color', ci95Color,  'LineWidth',1.2);
        line([s-med_hw, s+med_hw], [q(4), q(4)], 'Color', ci95Color,  'LineWidth',1.2);

        m = mean(y);
        line([s-med_hw, s+med_hw], [m, m], 'Color', meanColor, 'LineWidth',1.5);

        line([s-med_hw, s+med_hw], [q(3), q(3)], 'Color', medianColor,'LineWidth',1.5);
    end

    if plot_count==1
        hMed   = plot(nan, nan, '-', 'Color', medianColor, 'LineWidth',1.2);
        hMean  = plot(nan, nan, '-', 'Color', meanColor,   'LineWidth',1.2);
        hCI95  = plot(nan, nan, '-', 'Color', ci95Color,   'LineWidth',1.2);
        legend( [hViolin, hMed, hMean, hCI95], ...
                {'Density','50% Median','Mean','95% CI'}, ...
                'Location','northeast', 'Interpreter','tex' );
    end

    xticks(1:5); xticklabels(labels);
    xlim([0.5, 5.5]); grid on; box on;
    title(APAH14_names{k});
    ylabel('Average annual emission flux (kt/year)');
end


%% ===START: 图S4 — 16种 PAH 的小提琴图（中位数=0 的不画；2列排版） ===
figure;
valid_idx = [];
for k = 1:16
    if prctile(flux_total_PAH16(:,k),50) > 0
        valid_idx(end+1) = k; %
    end
end
nvalid = numel(valid_idx);
nrows  = ceil(nvalid/4);
plot_count = 0;

for t = 1:nvalid
    k = valid_idx(t);
    data_matrix = [ ...
        flux_nature_PAH16(:,k), ...
        flux_extraction_PAH16(:,k), ...
        flux_transport_PAH16(:,k), ...
        flux_consumption_PAH16(:,k), ...
        flux_total_PAH16(:,k) ];

    data_pos = data_matrix(data_matrix>0);
    if isempty(data_pos)
        continue;
    end

    plot_count = plot_count + 1;
    subplot(nrows,4,plot_count);

    lo = 10^floor(log10(min(data_pos)));
    hi = 10^ceil(log10(max(data_pos)));
    set(gca,'YScale','log','YLim',[lo,hi]);
    hold on;

    labels = {'Seepage','Extraction','Transportation','Consumption','Total'};
    violin_halfwidth = 0.3;
    rim_width = 0.12;

    for s = 1:5
        y = data_matrix(:,s); y = y(y>0);
        if numel(y)<2, continue; end

        xi_log = linspace(log10(min(y)), log10(max(y)), 2000);
        [f_log, ~] = ksdensity(log10(y), xi_log);
        xi = 10 .^ xi_log;
        f = f_log / max(f_log) * violin_halfwidth;

        hPatch = fill([s - f, fliplr(s + f)], [xi, fliplr(xi)], ...
                      [0.2 0.2 0.8], 'FaceAlpha',0.2, 'EdgeColor','none');
        if plot_count==1 && s==1
            hViolin = hPatch;
        end
 
        q = prctile(y, [0, 2.5, 50, 97.5, 100]);
        [~, idx_med] = min(abs(xi - q(3)));
        med_hw = f(idx_med);

        line([s-med_hw, s+med_hw], [q(2), q(2)], 'Color', ci95Color,  'LineWidth',1.2);
        line([s-med_hw, s+med_hw], [q(4), q(4)], 'Color', ci95Color,  'LineWidth',1.2);

        m = mean(y);
        line([s-med_hw, s+med_hw], [m, m], 'Color', meanColor, 'LineWidth',1.5);

        line([s-med_hw, s+med_hw], [q(3), q(3)], 'Color', medianColor,'LineWidth',1.5);
    end

    if plot_count==1
        hMed   = plot(nan, nan, '-', 'Color', medianColor, 'LineWidth',1.2);
        hMean  = plot(nan, nan, '-', 'Color', meanColor,   'LineWidth',1.2);
        hCI95  = plot(nan, nan, '-', 'Color', ci95Color,   'LineWidth',1.2);
        legend( [hViolin, hMed, hMean, hCI95], ...
                {'Density','50% Median','Mean','95% CI'}, ...
                'Location','northeast', 'Interpreter','tex' );
    end

    xticks(1:5); xticklabels(labels);
    xlim([0.5, 5.5]); grid on; box on;
    title(PAH16_names{k});
    ylabel('Average annual emission flux (kt/year)');
end
%

toc;