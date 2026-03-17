clearvars; clc;

tic;

% 读取新结果文件（变量名为 N×14 / N×16 / 名称列表）
load('Results_10000000.mat');   % 包含 flux_total_APAH14, flux_total_PAH16, APAH14_names, PAH16_names

% 补全全球清单基准
MeNAP_GEMS_Glabal = 69.483;
PAH40_GEMS_Glabal = 678.03;
PAH8_GEMS_Glabal = 455.715;

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


legendLabels = { ...
    '50% Median' , ...
    'Mean'         , ...
    '95% CI'    , ...
    % '100% Range'   ...
};


% === Fig.3：源别比值（MeNAP 与 Σ16PAHs） + Fig.S5（16个PAH逐个） ===

% === START: 若未定义 N，则按结果矩阵行数设定（不影响已有 N） ===
if ~exist('N','var') || isempty(N)
    if exist('flux_total_PAH16','var')
        N = size(flux_total_PAH16,1);
    elseif exist('flux_total_APAH14','var')
        N = size(flux_total_APAH14,1);
    end
end

% 
% 计算各源蒸发占比（仅 MeNAP 与 Σ16PAHs）
% 用新结果变量计算比值矩阵 
% 5个源：1-Seepage, 2-Extraction, 3-Transportation, 4-Consumption, 5-Total
% 2个指标：1-MeNAP（=APAH14第1列），2-Σ16PAHs（16物种求和）
ratios_sources = zeros(5, 2, N);

for i = 1:N
    % MeNAP：与其全球通量 MeNAP_GEMS_Glabal 的比值
    % Seepage
    ratios_sources(1,1,i) = flux_nature_APAH14(i,1)      / MeNAP_GEMS_Glabal;
    % Extraction
    ratios_sources(2,1,i) = flux_extraction_APAH14(i,1)  / MeNAP_GEMS_Glabal;
    % Transportation
    ratios_sources(3,1,i) = flux_transport_APAH14(i,1)   / MeNAP_GEMS_Glabal;
    % Consumption (陆地)
    ratios_sources(4,1,i) = flux_consumption_APAH14(i,1) / MeNAP_GEMS_Glabal;
    % Total
    ratios_sources(5,1,i) = flux_total_APAH14(i,1)       / MeNAP_GEMS_Glabal;

    % Σ16PAHs：对16个物种通量先求和，再与 PAH16_GEMS_Glabal 取比 
    s_nat   = sum(flux_nature_PAH16(i,:));
    s_ext   = sum(flux_extraction_PAH16(i,:));
    s_trans = sum(flux_transport_PAH16(i,:));
    s_cons  = sum(flux_consumption_PAH16(i,:));
    s_total = sum(flux_total_PAH16(i,:));

    ratios_sources(1,2,i) = s_nat   / PAH8_GEMS_Glabal;
    ratios_sources(2,2,i) = s_ext   / PAH8_GEMS_Glabal;
    ratios_sources(3,2,i) = s_trans / PAH8_GEMS_Glabal;
    ratios_sources(4,2,i) = s_cons  / PAH8_GEMS_Glabal;
    ratios_sources(5,2,i) = s_total / PAH8_GEMS_Glabal;
end
%


%% 图示3：MeNAP 与 Σ8PAHs 各源贡献（1×2 组图；）
figure('Position',[100 100 1200 400]);
sources = {'Seepage','Extraction','Transportation','Consumption','Total'};
nSources = numel(sources);
violin_half = 0.3;

% Colors
medianColor = [1 0 0];       % red
meanColor   = [1 0.6 0.6];   % light red
ci95Color   = [0, 0.45, 0.74];   % blue
fillColor   = [0.2 0.2 0.8]; % violin fill color  

set(groot, ...
    'defaultTextInterpreter','tex', ...
    'defaultAxesTickLabelInterpreter','tex', ...
    'defaultLegendInterpreter','tex');

% --- Subplot 1: MeNAP ---
subplot(1,2,1);
hold on;
hViolin = patch(nan, nan, fillColor,'FaceAlpha',0.2,'EdgeColor','none');

for i = 1:nSources
    y = squeeze(ratios_sources(i,1,:));
    y = y(y>0);
    [f_log, xi_log] = ksdensity(log10(y), 'NumPoints', 2000);
    xi = 10.^xi_log;
    f = f_log / max(f_log) * violin_half;
    hPatch = fill([i-f, fliplr(i+f)], [xi, fliplr(xi)], ...
         fillColor, 'FaceAlpha',0.2,'EdgeColor','none');
    
    if i == 1
        hViolin = hPatch;
    end

    pct = prctile(y,[2.5,50,97.5]);
    m   = mean(y);
    hCI95_l = line([i-violin_half,i+violin_half],[pct(1),pct(1)],...
                   'Color',ci95Color,'LineWidth',1.2);
    line([i-violin_half,i+violin_half],[pct(3),pct(3)],...
         'Color',ci95Color,'LineWidth',1.2);
    hMed_l = line([i-violin_half,i+violin_half],[pct(2),pct(2)],...
                  'Color',medianColor,'LineWidth',1.2);
    hMean_l= line([i-violin_half,i+violin_half],[m,m],...
                  'Color',meanColor,'LineWidth',1.2);
end

hold off;
set(gca, 'XTick',1:nSources, 'XTickLabel',sources, 'YScale','log');

% 百分比刻度（将比值乘100等效为指数+2）
yt = get(gca,'YTick');
exponents = round(log10(yt)) + 2;
labels = arrayfun(@(e) sprintf('10^{%d}', e), exponents, 'UniformOutput', false);
set(gca, 'YTickLabel', labels, 'TickLabelInterpreter','tex');
xlabel('Leakage sources','Interpreter','tex');
ylabel('Volatilization-to-global emission ratio (%)','Interpreter','tex');
title('MeNAPs contributions','Interpreter','tex');

legend([hViolin,hMed_l,hMean_l,hCI95_l], ...
       {'Density','50% Median','Mean','95% CI'}, ...
       'Location','northeast','Interpreter','tex');

grid on; box on;

% --- Subplot 2: Σ16 PAHs ---
subplot(1,2,2);
hold on;
hViolin = patch(nan, nan, fillColor,'FaceAlpha',0.2,'EdgeColor','none');
for i = 1:nSources
    y = squeeze(ratios_sources(i,2,:));
    y = y(y>0);
    [f_log, xi_log] = ksdensity(log10(y), 'NumPoints',2000);
    xi = 10.^xi_log;
    f = f_log / max(f_log) * violin_half;
    fill([i-f, fliplr(i+f)], [xi, fliplr(xi)], ...
         fillColor, 'FaceAlpha',0.2,'EdgeColor','none');
    pct = prctile(y,[2.5,50,97.5]);
    m   = mean(y);
    hCI95_l = line([i-violin_half,i+violin_half],[pct(1),pct(1)],...
                   'Color',ci95Color,'LineWidth',1.2);
    line([i-violin_half,i+violin_half],[pct(3),pct(3)],...
         'Color',ci95Color,'LineWidth',1.2);
    hMed_l = line([i-violin_half,i+violin_half],[pct(2),pct(2)],...
                  'Color',medianColor,'LineWidth',1.2);
    hMean_l= line([i-violin_half,i+violin_half],[m,m],...
                  'Color',meanColor,'LineWidth',1.2);
end

hold off;
set(gca, 'XTick',1:nSources, 'XTickLabel',sources, 'YScale','log');

yt = get(gca,'YTick');
exponents = round(log10(yt)) + 2;
labels = arrayfun(@(e) sprintf('10^{%d}', e), exponents, 'UniformOutput', false);
set(gca, 'YTickLabel', labels, 'TickLabelInterpreter','tex');

xlabel('Leakage sources','Interpreter','tex');
ylabel('Volatilization-to-global emission ratio (%)','Interpreter','tex');
title('\Sigma_{8} PAHs contributions','Interpreter','tex');

grid on; box on;


%% ===START: 图S5 — 16种 PAH 的“源别比值”小提琴图（中位数=0 的不画；2列排版）

% 从 Excel 读取 16 个物种全球通量
PAH16_global_2021 = readmatrix('APAHs15 & oil calculated 3.0 20251023.xlsx', ...
    'Sheet','Data from GEMS','Range','AK2:AK17');

% Excel 给定的固定顺序
PAH16_order = {'NAP','PHE','FLO','CHR','ACY','ACE','ANT','FLA','PYR', ...
               'BaA','BbF','BkF','BaP','DahA','IcdP','BghiP'};

% 将 Excel 的全球通量按 PAH16_names 的列顺序映射，得到每一列对应的分母
PAH16_global_denom = nan(1,16);

for kk = 1:16
    nm = PAH16_names{kk};

    id = find(strcmpi(nm, PAH16_order), 1);
    if isempty(id)
        error('Fig.S5 denominator mapping failed: PAH16_names{%d} = %s not found in Excel order list.', kk, PAH16_names{kk});
    end

    PAH16_global_denom(kk) = PAH16_global_2021(id);
end

fprintf('\n[CHECK] Fig.S5 denominators (species-specific global flux):\n');
disp(table(PAH16_names(:), PAH16_global_denom(:), ...
    'VariableNames',{'PAH16_names','GlobalFluxDenom_fromGEMS'}));


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

    % 构造 5×N 的比值（源别）
    den_k = PAH16_global_denom(k);

    r_nat   = flux_nature_PAH16(:,k)      / den_k;
    r_ext   = flux_extraction_PAH16(:,k)  / den_k;
    r_trans = flux_transport_PAH16(:,k)   / den_k;
    r_cons  = flux_consumption_PAH16(:,k) / den_k;
    r_total = flux_total_PAH16(:,k)       / den_k;
    
    R = [r_nat, r_ext, r_trans, r_cons, r_total];  % N×5

    % 若全≤0 则跳过
    data_pos = R(R>0);
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
    violin_half = 0.3;

    for s = 1:5
        y = R(:,s); y = y(y>0);
        if numel(y)<2, continue; end

        [f_log, xi_log] = ksdensity(log10(y), linspace(log10(min(y)), log10(max(y)), 2000));
        xi = 10.^xi_log;
        f  = f_log / max(f_log) * violin_half;

        hPatch = fill([s - f, fliplr(s + f)], [xi, fliplr(xi)], ...
                      fillColor, 'FaceAlpha',0.2, 'EdgeColor','none');
        if plot_count==1 && s==1
            hViolinS5 = hPatch; % 用于图例
        end

        pct = prctile(y,[2.5,50,97.5]);
        m   = mean(y);
        hCI95_lS5 = line([s-violin_half,s+violin_half],[pct(1),pct(1)],'Color',ci95Color,'LineWidth',1.2);
        line([s-violin_half,s+violin_half],[pct(3),pct(3)],'Color',ci95Color,'LineWidth',1.2);
        hMed_lS5  = line([s-violin_half,s+violin_half],[pct(2),pct(2)],'Color',medianColor,'LineWidth',1.2);
        hMean_lS5 = line([s-violin_half,s+violin_half],[m,m],           'Color',meanColor, 'LineWidth',1.2);
    end

    % 百分比刻度标签
    yt = get(gca,'YTick');
    exponents = round(log10(yt)) + 2;
    ylabels = arrayfun(@(e) sprintf('10^{%d}', e), exponents, 'UniformOutput', false);
    set(gca,'YTickLabel',ylabels,'TickLabelInterpreter','tex');

    xticks(1:5); xticklabels(labels);
    xlim([0.5, 5.5]); grid on; box on;
    title(PAH16_names{k}, 'Interpreter','tex');
    ylabel('Volatilization-to-global emission ratio (%)','Interpreter','tex');
end

% 只在第一个子图画一次图例
if exist('hViolinS5','var')
    legend([hViolinS5,hMed_lS5,hMean_lS5,hCI95_lS5], ...
           {'Density','50% Median','Mean','95% CI'}, ...
           'Location','northeast','Interpreter','tex');
end
%

toc;
