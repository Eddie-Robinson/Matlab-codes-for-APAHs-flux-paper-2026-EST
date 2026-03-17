% Figure4_SRC_new.m
% 读取蒙特卡洛模拟结果并绘制 SRC（Fig.4 / Fig.S6 / Fig.S7）及 MeNAP 的 SRC^2 饼图（Fig.5）
clearvars; clc;

tic;
% 1) 加载新结果文件
load('Results_10000000.mat');   


% 2) 补齐全球清单基准
MeNAP_GEMS_Glabal = 69.483;
PAH8_GEMS_Glabal = 455.715;
% === MOD END ===

% 统一的 X 变量名
varNames = { ...
  'Nature spill','Extraction spill','Transportation spill','Consumption spill', ...
  'Ratio_{light oil}','Ratio_{medium oil}', ...
  'Concentration_{light}','Concentration_{medium}','Concentration_{heavy}', ...
  'Evaporation coef_{Sea, light}','Evaporation coef_{Sea, medium}','Evaporation coef_{Sea, heavy}', ...
  'Land-based factor (f_{land})' ...
};

% Bootstrap 次数
B = 1000;

% 初始化 CSV 输出容器与文件名 
SRC_results_all = table();      % 汇总 Fig4 + S6 + S7 的 SRC(中位数+95%CI)
SRC2_share_MeNAP = table();     % Fig5 的 SRC^2 份额（可选，但建议一起导出）

% 用 N / B / 时间戳命名，避免覆盖
if exist('N','var') && ~isempty(N)
    N_tag = sprintf('N%d', N);
else
    N_tag = 'Nunknown';
end
ts_tag = datestr(now,'yyyymmdd_HHMMSS');

csv_SRC_all_name  = sprintf('SRC_All_%s_%s_B%d.csv', N_tag, ts_tag, B);
csv_SRC2_name     = sprintf('SRC2Share_MeNAP_%s_%s_B%d.csv', N_tag, ts_tag, B);


% 筛选“非0物种”（中位数>0）
idx_APAH_nz = find(prctile(flux_total_APAH14,50,1) > 0);   % 期望 8 个
idx_PAH_nz  = find(prctile(flux_total_PAH16, 50,1) > 0);   % 期望 8 个

fprintf('\n[CHECK] Non-zero APAH14 species (median>0):\n');
disp(APAH14_names(idx_APAH_nz));

fprintf('\n[CHECK] Non-zero PAH16 species (median>0):\n');
disp(PAH16_names(idx_PAH_nz));


%% Fig.4：MeNAP（= APAH14 第1列：C1NAP） 的 SRC + 95% CI

% 指定 Fig.4 分析对象 = MeNAP (APAH14 column 1) 
k = 1;  % MeNAP = C1NAP
fprintf('\n[Fig4] SRC target: MeNAP (APAH14 col %d = %s)\n', k, APAH14_names{k});

Y = flux_total_APAH14(:,k);

X = [ ...
    Oil_em_Nature_MC, ...
    Oil_em_Ext_MC, Oil_em_Trans_MC, Oil_em_Cons_MC, ...
    light_ratio, medium_ratio, ...
    oil_APAH14.light(:,k), oil_APAH14.medium(:,k), oil_APAH14.heavy(:,k), ...
    squeeze(sea_env_coef_APAH14(:,1,k)), squeeze(sea_env_coef_APAH14(:,2,k)), squeeze(sea_env_coef_APAH14(:,3,k)), ...
    land_factor ...
];


[N, p] = size(X);

% 标准化（）
X_std = (X - mean(X)) ./ std(X);
Y_std = (Y - mean(Y)) ./ std(Y);

%  防止极端情况 std=0 导致 NaN
if any(std(X)==0) || std(Y)==0
    error('SRC input has zero std (constant variable). Please check inputs for MeNAP.');
end

% Bootstrap
Beta_boot = nan(B,p);
for b = 1:B
    idx = randi(N, N, 1);
    Beta_boot(b,:) = (X_std(idx,:) \ Y_std(idx))';
end

beta_med = median(Beta_boot,1);
ci_lower = prctile(Beta_boot,2.5,1);
ci_upper = prctile(Beta_boot,97.5,1);

% 输出表格核对
disp(table(varNames', beta_med', ci_lower', ci_upper', ...
    'VariableNames',{'Variable','MedianSRC','CI2.5pct','CI97.5pct'}));

% 追加 Fig.4(MeNAP) 的 SRC 结果到总表
figTag     = 'Fig4';
groupTag   = 'APAH14';
speciesCol = k;
speciesTag = 'MeNAP';  % 固定 MeNAP

Ttmp = table();
Ttmp.Figure      = repmat({figTag},   p, 1);
Ttmp.SpeciesGroup= repmat({groupTag}, p, 1);
Ttmp.SpeciesCol  = repmat(speciesCol,p, 1);
Ttmp.SpeciesName = repmat({speciesTag},p,1);
Ttmp.Variable    = varNames';
Ttmp.MedianSRC   = beta_med';
Ttmp.CI2_5pct    = ci_lower';
Ttmp.CI97_5pct   = ci_upper';
Ttmp.BootstrapN  = repmat(B, p, 1);

SRC_results_all = [SRC_results_all; Ttmp];



%% ===  Fig.4 画图 ===
figure; clf; hold on;

barh(1:p, beta_med, 0.6, 'FaceColor',[0.7 0.85 1]);

errorbar(beta_med, 1:p, beta_med-ci_lower, ci_upper-beta_med, ...
    'horizontal','LineStyle','none','Color','k','LineWidth',2,'CapSize',10);

set(gca, ...
    'YTick',       1:p, ...
    'YTickLabel',  varNames, ...
    'YDir','reverse' ...
);
xlabel('Standardized Regression Coefficient');
title(sprintf('SRC for MeNAP'));
xline(0,'--k');
xlim([-0.1, 1.0]);
grid on; box on;


%% Fig.5：MeNAP 的 SRC^2 份额饼图

beta_sq   = beta_med .^ 2;
share_raw = beta_sq ./ sum(beta_sq);

thresh = 0.01;
idxTop = find(share_raw >= thresh);
idxEtc = find(share_raw  < thresh);

share_plot  = [share_raw(idxTop),  sum(share_raw(idxEtc))];
labels_plot = [varNames(idxTop), {'Others (<1% each)'}];

fprintf('\nIndependent variance contributions (SRC^2) for MeNAP:\n');
disp(table(labels_plot', share_plot'*100, ...
     'VariableNames',{'Variable','Share_percent'}));

% Fig.5(MeNAP) 的 SRC^2 份额表
nShare = numel(labels_plot);
Tshare = table();
Tshare.Figure      = repmat({'Fig5'}, nShare, 1);
Tshare.SpeciesName = repmat({'MeNAP'}, nShare, 1);
Tshare.Variable    = labels_plot(:);
Tshare.Share_percent = share_plot(:) * 100;

SRC2_share_MeNAP = Tshare;  % 


figure; clf;
pie(share_plot, labels_plot);
title(sprintf('SRC^2 share for MeNAP'));
%



%% Fig.S6：其余 7 个非0 APAH 的 SRC（2列排版，最后一个居中）


%  除 MeNAP 外的 7 个非0 APAH，按 idx_APAH_nz 顺序
idx_APAH_other = idx_APAH_nz;
idx_APAH_other(idx_APAH_other==1) = [];   % 去掉 MeNAP (col1)
fprintf('\n===== SRC for APAH8 (excluding MeNAP) =====\n');
disp(APAH14_names(idx_APAH_other));
% 

figure;  % Fig.S6
nA = numel(idx_APAH_other);
nrows = ceil(nA/4);
plot_count = 0;

for t = 1:nA
    k = idx_APAH_other(t);
    plot_count = plot_count + 1;

    % 本子图对应的物种
    % S6-APAH #%d: APAH14 col k = %s
    fprintf('\n[S6] SRC for APAH14 species %d/%d (col %d): %s\n', t, nA, k, APAH14_names{k});

    Y = flux_total_APAH14(:,k);

    X = [ ...
        Oil_em_Nature_MC, ...
        Oil_em_Ext_MC, Oil_em_Trans_MC, Oil_em_Cons_MC, ...
        light_ratio, medium_ratio, ...
        oil_APAH14.light(:,k), oil_APAH14.medium(:,k), oil_APAH14.heavy(:,k), ...
        squeeze(sea_env_coef_APAH14(:,1,k)), squeeze(sea_env_coef_APAH14(:,2,k)), squeeze(sea_env_coef_APAH14(:,3,k)), ...
        land_factor ...
    ];

    [N, p] = size(X);

    % X_std = (X - mean(X)) ./ std(X);
    % Y_std = (Y - mean(Y)) ./ std(Y);

    % MEM 原地标准化，避免额外生成 X_std / Y_std 两个大矩阵
    X = (X - mean(X)) ./ std(X);
    Y = (Y - mean(Y)) ./ std(Y);

    % 保险：若该物种实际 std=0，跳过
    if any(std(X)==0) || std(Y)==0
        fprintf('[S6] Skip %s due to zero std.\n', APAH14_names{k});
        continue;
    end

    Beta_boot = nan(B,p);
    for b = 1:B
        idx = randi(N, N, 1);
        % MEM使用原地标准化后的 X/Y，不再引用 X_std/Y_std 
        Beta_boot(b,:) = (X(idx,:) \ Y(idx))';
        
    end

    % Beta_boot = nan(B,p);
    % for b = 1:B
    %     idx = randi(N, N, 1);
    %     Beta_boot(b,:) = (X_std(idx,:) \ Y_std(idx))';
    % end

    beta_med = median(Beta_boot,1);
    ci_lower = prctile(Beta_boot,2.5,1);
    ci_upper = prctile(Beta_boot,97.5,1);

    % 输出表格核对
    disp(table(varNames', beta_med', ci_lower', ci_upper', ...
        'VariableNames',{'Variable','MedianSRC','CI2.5pct','CI97.5pct'}));
    
    % 加 S6 的 SRC 结果到总表
    figTag     = 'FigS6';
    groupTag   = 'APAH14';
    speciesCol = k;
    speciesTag = APAH14_names{k};

    Ttmp = table();
    Ttmp.Figure       = repmat({figTag},    p, 1);
    Ttmp.SpeciesGroup = repmat({groupTag},  p, 1);
    Ttmp.SpeciesCol   = repmat(speciesCol,  p, 1);
    Ttmp.SpeciesName  = repmat({speciesTag},p, 1);
    Ttmp.Variable     = varNames';
    Ttmp.MedianSRC    = beta_med';
    Ttmp.CI2_5pct     = ci_lower';
    Ttmp.CI97_5pct    = ci_upper';
    Ttmp.BootstrapN   = repmat(B, p, 1);

    SRC_results_all = [SRC_results_all; Ttmp];


    subplot(nrows,4,plot_count); hold on;

    barh(1:p, beta_med, 0.6, 'FaceColor',[0.7 0.85 1]);
    errorbar(beta_med, 1:p, beta_med-ci_lower, ci_upper-beta_med, ...
        'horizontal','LineStyle','none','Color','k','LineWidth',2,'CapSize',10);

    set(gca, ...
        'YTick',       1:p, ...
        'YTickLabel',  varNames, ...
        'YDir','reverse' ...
    );
    xlabel('Standardized Regression Coefficient');
    title(sprintf('SRC for %s', APAH14_names{k}));
    xline(0,'--k');
    xlim([-0.1, 1.0]);
    grid on; box on;

    % MEM  清理本物种临时大变量，避免循环累计/内存碎片
    clear X Y Beta_boot idx
    drawnow;
    
end

% 若 S6 子图数量为奇数，让最后一个居中（不改子图内容，只调位置） ===
% if mod(nA,2)==1 && plot_count>=1
%     ax = subplot(nrows,2,plot_count);
%     pos = get(ax,'Position');
%     % 将最后一个子图在横向居中（经验值：稍微加宽/居中）
%     pos(1) = 0.30;
%     pos(3) = 0.40;
%     set(ax,'Position',pos);
% end



%% Fig.S7：8 个非0 PAH 的 SRC（2列排版）

fprintf('\n===== SRC for PAH8 =====\n');
disp(PAH16_names(idx_PAH_nz));

figure;  % Fig.S7
nP = numel(idx_PAH_nz);
nrows = ceil(nP/4);
plot_count = 0;

for t = 1:nP
    k = idx_PAH_nz(t);
    plot_count = plot_count + 1;

    % --- 本子图对应的物种 ---
    % S7-PAH t: PAH16 col k = 名称
    fprintf('\n[S7] SRC for PAH16 species %d/%d (col %d): %s\n', t, nP, k, PAH16_names{k});
    
    Y = flux_total_PAH16(:,k);

    X = [ ...
        Oil_em_Nature_MC, ...
        Oil_em_Ext_MC, Oil_em_Trans_MC, Oil_em_Cons_MC, ...
        light_ratio, medium_ratio, ...
        oil_PAH16.light(:,k), oil_PAH16.medium(:,k), oil_PAH16.heavy(:,k), ...
        squeeze(sea_env_coef_PAH16(:,1,k)), squeeze(sea_env_coef_PAH16(:,2,k)), squeeze(sea_env_coef_PAH16(:,3,k)), ...
        land_factor ...
    ];

    [N, p] = size(X);

    % X_std = (X - mean(X)) ./ std(X);
    % Y_std = (Y - mean(Y)) ./ std(Y);

    % === MOD MEM START: 原地标准化，避免额外生成 X_std / Y_std 两个大矩阵 ===
    X = (X - mean(X)) ./ std(X);
    Y = (Y - mean(Y)) ./ std(Y);

    if any(std(X)==0) || std(Y)==0
        fprintf('[S7] Skip %s due to zero std.\n', PAH16_names{k});
        continue;
    end

    Beta_boot = nan(B,p);
    for b = 1:B
        idx = randi(N, N, 1);
        % === MEM 使用原地标准化后的 X/Y，不再引用 X_std/Y_std ===
        Beta_boot(b,:) = (X(idx,:) \ Y(idx))';
        % Beta_boot(b,:) = (X_std(idx,:) \ Y_std(idx))';
    end

    beta_med = median(Beta_boot,1);
    ci_lower = prctile(Beta_boot,2.5,1);
    ci_upper = prctile(Beta_boot,97.5,1);

    % 输出表格核对
    disp(table(varNames', beta_med', ci_lower', ci_upper', ...
        'VariableNames',{'Variable','MedianSRC','CI2.5pct','CI97.5pct'}));
    
    % === CSV START: 追加 S7 的 SRC 结果到总表 ===
    figTag     = 'FigS7';
    groupTag   = 'PAH16';
    speciesCol = k;
    speciesTag = PAH16_names{k};
    
    Ttmp = table();
    Ttmp.Figure       = repmat({figTag},    p, 1);
    Ttmp.SpeciesGroup = repmat({groupTag},  p, 1);
    Ttmp.SpeciesCol   = repmat(speciesCol,  p, 1);
    Ttmp.SpeciesName  = repmat({speciesTag},p, 1);
    Ttmp.Variable     = varNames';
    Ttmp.MedianSRC    = beta_med';
    Ttmp.CI2_5pct     = ci_lower';
    Ttmp.CI97_5pct    = ci_upper';
    Ttmp.BootstrapN   = repmat(B, p, 1);
    
    SRC_results_all = [SRC_results_all; Ttmp];
    % === CSV END ===

    subplot(nrows,4,plot_count); hold on;

    barh(1:p, beta_med, 0.6, 'FaceColor',[0.7 0.85 1]);
    errorbar(beta_med, 1:p, beta_med-ci_lower, ci_upper-beta_med, ...
        'horizontal','LineStyle','none','Color','k','LineWidth',2,'CapSize',10);

    set(gca, ...
        'YTick',       1:p, ...
        'YTickLabel',  varNames, ...
        'YDir','reverse' ...
    );
    xlabel('Standardized Regression Coefficient');
    title(sprintf('SRC for %s', PAH16_names{k}));
    xline(0,'--k');
    xlim([-0.1, 1.0]);
    grid on; box on;

    % === MEM: 清理本物种临时大变量 ===
    clear X Y Beta_boot idx
    drawnow;

end

% === CSV START: 导出 CSV（含类别字段，便于回溯） ===
writetable(SRC_results_all, csv_SRC_all_name);
fprintf('\n[CSV] SRC results saved: %s\n', csv_SRC_all_name);

% Fig.5 的 SRC^2 份额（若存在）
if ~isempty(SRC2_share_MeNAP)
    writetable(SRC2_share_MeNAP, csv_SRC2_name);
    fprintf('[CSV] SRC^2 share saved: %s\n', csv_SRC2_name);
end
% === CSV END ===


toc;