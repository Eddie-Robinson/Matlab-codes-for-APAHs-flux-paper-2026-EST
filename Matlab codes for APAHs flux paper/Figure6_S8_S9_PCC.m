% Figure6_S8_S9_PCC_new.m
% 读取蒙特卡洛模拟结果并绘制 PCC（Fig.6 / Fig.S8 / Fig.S9）
clearvars; clc;

tic;
% 1) 加载新结果文件
load('Results_10000000.mat');   % <


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

% 初始化 PCC 输出容器与文件名
PCC_results_all = table();   % 汇总 Fig6 + S8 + S9 的 PCC(中位数+95%CI)

% 用 N / B / 时间戳命名，避免覆盖
if exist('N','var') && ~isempty(N)
    N_tag = sprintf('N%d', N);
else
    N_tag = 'Nunknown';
end
ts_tag = datestr(now,'yyyymmdd_HHMMSS');

csv_PCC_all_name = sprintf('PCC_All_%s_%s_B%d.csv', N_tag, ts_tag, B);

% 筛选“非0物种”（中位数>0），避免给全0算 PCC
idx_APAH_nz = find(prctile(flux_total_APAH14,50,1) > 0);   % 期望 8 个
idx_PAH_nz  = find(prctile(flux_total_PAH16, 50,1) > 0);   % 期望 8 个

fprintf('\n[CHECK] Non-zero APAH14 species (median>0):\n');
disp(APAH14_names(idx_APAH_nz));

fprintf('\n[CHECK] Non-zero PAH16 species (median>0):\n');
disp(PAH16_names(idx_PAH_nz));
% 

%% Fig.6：MeNAP（= APAH14 第1列：C1NAP） 的 PCC + 95% CI

% 指定 Fig.6 分析对象 = MeNAP (APAH14 column 1)
k = 1;  % MeNAP = C1NAP
fprintf('\n[Fig6] PCC target: MeNAP (APAH14 col %d = %s)\n', k, APAH14_names{k});

Y = flux_total_APAH14(:,k);

X = [ ...
    Oil_em_Nature_MC, ...
    Oil_em_Ext_MC, ...
    Oil_em_Trans_MC, ...
    Oil_em_Cons_MC, ...
    light_ratio, ...
    medium_ratio, ...
    oil_APAH14.light(:,k), ...
    oil_APAH14.medium(:,k), ...
    oil_APAH14.heavy(:,k), ...
    squeeze(sea_env_coef_APAH14(:,1,k)), ...
    squeeze(sea_env_coef_APAH14(:,2,k)), ...
    squeeze(sea_env_coef_APAH14(:,3,k)), ...
    land_factor ...
];
% 
[N,p] = size(X);

% --- Fig.6：PCC（Partial Correlation Coefficients） ---
% 1) Bootstrap 计算 PCC 分布
P_boot = nan(B,p);
for b = 1:B
    idx = randi(N, N, 1);
    P_boot(b,:) = partialcorri(Y(idx), X(idx,:));
end

% 2) 统计中位数和 95% CI
P_med     = median(P_boot, 1);
ci_lower  = prctile(P_boot, 2.5, 1);
ci_upper  = prctile(P_boot, 97.5,1);

% 3) 绘图
figure; clf; hold on;
hBar = barh(1:p, P_med, 0.6, 'FaceColor',[0.7 0.85 1.0]);

errorbar(P_med, 1:p, P_med-ci_lower, ci_upper-P_med, ...
    'horizontal','LineStyle','none','Color','k','LineWidth',2,'CapSize',10);

disp(table(varNames', P_med', ci_lower', ci_upper', ...
    'VariableNames',{'Var','MedianPCC','CI2.5pct','CI97.5pct'}));

% 追加 Fig.6(MeNAP) 的 PCC 结果到总表 
figTag     = 'Fig6';
groupTag   = 'APAH14';
speciesCol = k;
speciesTag = 'MeNAP';  % 固定写 MeNAP

Ttmp = table();
Ttmp.Figure       = repmat({figTag},    p, 1);
Ttmp.SpeciesGroup = repmat({groupTag},  p, 1);
Ttmp.SpeciesCol   = repmat(speciesCol,  p, 1);
Ttmp.SpeciesName  = repmat({speciesTag},p, 1);
Ttmp.Variable     = varNames';
Ttmp.MedianPCC    = P_med';
Ttmp.CI2_5pct     = ci_lower';
Ttmp.CI97_5pct    = ci_upper';
Ttmp.BootstrapN   = repmat(B, p, 1);

PCC_results_all = [PCC_results_all; Ttmp];
% === CSV END ===

set(gca, ...
    'YDir',        'reverse', ...
    'YTick',       1:p, ...
    'YTickLabel',  varNames);
xlabel('Partial Correlation Coefficient');
title(sprintf('PCC for MeNAP'));   % 
xline(0,'--k');
xlim([-0.1, 1.0]);
grid on; box on;

%  MEM 清理 Fig.6 临时大变量，避免占用内存影响后续 S8/S9
clear X Y P_boot idx
drawnow;

%% Fig.S8：其余 7 个非0 APAH 的 PCC（2列排版，最后一个居中）

% 选择“除 MeNAP 外”的 7 个非0 APAH，按 idx_APAH_nz 顺序
idx_APAH_other = idx_APAH_nz;
idx_APAH_other(idx_APAH_other==1) = [];   % 去掉 MeNAP (col1)
fprintf('\n===== PCC for APAH8 (excluding MeNAP) =====\n');
disp(APAH14_names(idx_APAH_other));
% 

figure;  % Fig.S8
nA = numel(idx_APAH_other);
nrows = ceil(nA/4);
plot_count = 0;

for t = 1:nA
    k = idx_APAH_other(t);
    plot_count = plot_count + 1;

    %  命令行输出用于检查
    fprintf('\n[S8] PCC for APAH14 species %d/%d (col %d): %s\n', t, nA, k, APAH14_names{k});

    Y = flux_total_APAH14(:,k);

    X = [ ...
        Oil_em_Nature_MC, ...
        Oil_em_Ext_MC, ...
        Oil_em_Trans_MC, ...
        Oil_em_Cons_MC, ...
        light_ratio, ...
        medium_ratio, ...
        oil_APAH14.light(:,k), ...
        oil_APAH14.medium(:,k), ...
        oil_APAH14.heavy(:,k), ...
        squeeze(sea_env_coef_APAH14(:,1,k)), ...
        squeeze(sea_env_coef_APAH14(:,2,k)), ...
        squeeze(sea_env_coef_APAH14(:,3,k)), ...
        land_factor ...
    ];

    [N,p] = size(X);

    % Bootstrap
    P_boot = nan(B,p);
    for b = 1:B
        idx = randi(N, N, 1);
        P_boot(b,:) = partialcorri(Y(idx), X(idx,:));
    end

    P_med     = median(P_boot, 1);
    ci_lower  = prctile(P_boot, 2.5, 1);
    ci_upper  = prctile(P_boot, 97.5,1);

    % 输出表格核对
    disp(table(varNames', P_med', ci_lower', ci_upper', ...
        'VariableNames',{'Var','MedianPCC','CI2.5pct','CI97.5pct'}));

    % 追加 Fig.S8 的 PCC 结果到总表
    figTag     = 'FigS8';
    groupTag   = 'APAH14';
    speciesCol = k;
    speciesTag = APAH14_names{k};
    
    Ttmp = table();
    Ttmp.Figure       = repmat({figTag},    p, 1);
    Ttmp.SpeciesGroup = repmat({groupTag},  p, 1);
    Ttmp.SpeciesCol   = repmat(speciesCol,  p, 1);
    Ttmp.SpeciesName  = repmat({speciesTag},p, 1);
    Ttmp.Variable     = varNames';
    Ttmp.MedianPCC    = P_med';
    Ttmp.CI2_5pct     = ci_lower';
    Ttmp.CI97_5pct    = ci_upper';
    Ttmp.BootstrapN   = repmat(B, p, 1);
    
    PCC_results_all = [PCC_results_all; Ttmp];
    % === CSV END ===

    subplot(nrows,4,plot_count); hold on;

    % 绘图风格保持
    hBar = barh(1:p, P_med, 0.6, 'FaceColor',[0.7 0.85 1.0]);
    errorbar(P_med, 1:p, P_med-ci_lower, ci_upper-P_med, ...
        'horizontal','LineStyle','none','Color','k','LineWidth',2,'CapSize',10);

    set(gca, ...
        'YDir',        'reverse', ...
        'YTick',       1:p, ...
        'YTickLabel',  varNames);
    xlabel('Partial Correlation Coefficient');
    title(sprintf('PCC for %s', APAH14_names{k}));   % === MOD: 子图标题 ===
    xline(0,'--k');
    xlim([-0.1, 1.0]);
    grid on; box on;

    % MEM清理本次循环临时大变量，避免循环累计/内存碎片
    clear X Y P_boot idx
    drawnow;

end

% 若 S8 子图数量为奇数，让最后一个居中
% if mod(nA,2)==1 && plot_count>=1
%     ax = subplot(nrows,2,plot_count);
%     pos = get(ax,'Position');
%     pos(1) = 0.30;
%     pos(3) = 0.40;
%     set(ax,'Position',pos);
% end
% 


%% Fig.S9：8 个非0 PAH 的 PCC（2列排版）

fprintf('\n===== PCC for PAH8 =====\n');
disp(PAH16_names(idx_PAH_nz));

figure;  % Fig.S9
nP = numel(idx_PAH_nz);
nrows = ceil(nP/4);
plot_count = 0;

for t = 1:nP
    k = idx_PAH_nz(t);
    plot_count = plot_count + 1;

    % 命令行输出用于检查
    fprintf('\n[S9] PCC for PAH16 species %d/%d (col %d): %s\n', t, nP, k, PAH16_names{k});

    Y = flux_total_PAH16(:,k);

    X = [ ...
        Oil_em_Nature_MC, ...
        Oil_em_Ext_MC, ...
        Oil_em_Trans_MC, ...
        Oil_em_Cons_MC, ...
        light_ratio, ...
        medium_ratio, ...
        oil_PAH16.light(:,k), ...
        oil_PAH16.medium(:,k), ...
        oil_PAH16.heavy(:,k), ...
        squeeze(sea_env_coef_PAH16(:,1,k)), ...
        squeeze(sea_env_coef_PAH16(:,2,k)), ...
        squeeze(sea_env_coef_PAH16(:,3,k)), ...
        land_factor ...
    ];

    [N,p] = size(X);

    % Bootstrap
    P_boot = nan(B,p);
    for b = 1:B
        idx = randi(N, N, 1);
        P_boot(b,:) = partialcorri(Y(idx), X(idx,:));
    end

    P_med     = median(P_boot, 1);
    ci_lower  = prctile(P_boot, 2.5, 1);
    ci_upper  = prctile(P_boot, 97.5,1);

    disp(table(varNames', P_med', ci_lower', ci_upper', ...
        'VariableNames',{'Var','MedianPCC','CI2.5pct','CI97.5pct'}));

    % 追加 Fig.S9 的 PCC 结果到总表
    figTag     = 'FigS9';
    groupTag   = 'PAH16';
    speciesCol = k;
    speciesTag = PAH16_names{k};
    
    Ttmp = table();
    Ttmp.Figure       = repmat({figTag},    p, 1);
    Ttmp.SpeciesGroup = repmat({groupTag},  p, 1);
    Ttmp.SpeciesCol   = repmat(speciesCol,  p, 1);
    Ttmp.SpeciesName  = repmat({speciesTag},p, 1);
    Ttmp.Variable     = varNames';
    Ttmp.MedianPCC    = P_med';
    Ttmp.CI2_5pct     = ci_lower';
    Ttmp.CI97_5pct    = ci_upper';
    Ttmp.BootstrapN   = repmat(B, p, 1);
    
    PCC_results_all = [PCC_results_all; Ttmp];
    

    subplot(nrows,4,plot_count); hold on;

    hBar = barh(1:p, P_med, 0.6, 'FaceColor',[0.7 0.85 1.0]);
    errorbar(P_med, 1:p, P_med-ci_lower, ci_upper-P_med, ...
        'horizontal','LineStyle','none','Color','k','LineWidth',2,'CapSize',10);

    set(gca, ...
        'YDir',        'reverse', ...
        'YTick',       1:p, ...
        'YTickLabel',  varNames);
    xlabel('Partial Correlation Coefficient');
    title(sprintf('PCC for %s', PAH16_names{k}));    %  子图标题 
    xline(0,'--k');
    xlim([-0.1, 1.0]);
    grid on; box on;

    % MEM清理本次循环临时大变量，避免循环累计/内存碎片
    clear X Y P_boot idx
    drawnow;

end

% 导出 PCC 汇总 CSV
writetable(PCC_results_all, csv_PCC_all_name);


fprintf('\n[CSV] PCC results saved: %s\n', csv_PCC_all_name);


toc;
