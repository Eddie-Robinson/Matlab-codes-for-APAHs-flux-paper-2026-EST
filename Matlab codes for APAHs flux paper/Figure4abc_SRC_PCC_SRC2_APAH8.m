% Figure_SRC_PCC_APAH8sum.m
% 功能：对“非零 8 个 APAHs 总和通量（Σ8 APAHs）”计算 SRC & PCC，并绘制 1×3 组图（SRC / PCC / SRC^2 share）
clearvars; clc;
tic;

%  1) 加载结果文件
load('Results_10000000.mat');

% 统一变量名
varNames = { ...
  'Nature spill','Extraction spill','Transportation spill','Consumption spill', ...
  'Ratio_{light oil}','Ratio_{medium oil}', ...
  'Concentration_{light}','Concentration_{medium}','Concentration_{heavy}', ...
  'Evaporation coef_{Sea, light}','Evaporation coef_{Sea, medium}','Evaporation coef_{Sea, heavy}', ...
  'Land-based factor (f_{land})' ...
};

B = 1000;

% 筛选“非0物种”（中位数>0）
idx_APAH_nz = find(prctile(flux_total_APAH14,50,1) > 0);   % 期望 8 个
fprintf('\n[CHECK] Σ8 APAHs (median>0) include:\n');
disp(APAH14_names(idx_APAH_nz));


% 目标：Σ8 APAHs 总通量
% 
fprintf('\n[TARGET] Σ8 APAHs total flux (sum of non-zero species)\n');

%  Y = Σ8 APAHs 总和通量 
Y = zeros(size(flux_total_APAH14,1),1);
for ii = 1:numel(idx_APAH_nz)
    k = idx_APAH_nz(ii);
    Y = Y + flux_total_APAH14(:,k);
end
% 

N = size(Y,1);
p = 13;

% 
% 构造 X（避免 N×8 临时大矩阵，用逐列累加）
%

% 1) 共享参数（原样）
X = zeros(N,p);
X(:,1)  = Oil_em_Nature_MC;
X(:,2)  = Oil_em_Ext_MC;
X(:,3)  = Oil_em_Trans_MC;
X(:,4)  = Oil_em_Cons_MC;
X(:,5)  = light_ratio;
X(:,6)  = medium_ratio;
X(:,13) = land_factor;

% 2) Σ8 浓度（逐列累加，避免 oil_APAH14.light(:,idx) 形成 N×8 临时矩阵）
conc_light_sum  = zeros(N,1);
conc_medium_sum = zeros(N,1);
conc_heavy_sum  = zeros(N,1);

for ii = 1:numel(idx_APAH_nz)
    k = idx_APAH_nz(ii);
    conc_light_sum  = conc_light_sum  + oil_APAH14.light(:,k);
    conc_medium_sum = conc_medium_sum + oil_APAH14.medium(:,k);
    conc_heavy_sum  = conc_heavy_sum  + oil_APAH14.heavy(:,k);
end

X(:,7) = conc_light_sum;
X(:,8) = conc_medium_sum;
X(:,9) = conc_heavy_sum;

% 3) 有效海洋蒸发系数（浓度加权平均）
%    evap_eff = sum(conc_k * coef_k) / sum(conc_k)
%    同样逐列累加，避免 N×8 临时矩阵

num_L = zeros(N,1); den_L = zeros(N,1);
num_M = zeros(N,1); den_M = zeros(N,1);
num_H = zeros(N,1); den_H = zeros(N,1);

for ii = 1:numel(idx_APAH_nz)
    k = idx_APAH_nz(ii);

    wL = oil_APAH14.light(:,k);
    wM = oil_APAH14.medium(:,k);
    wH = oil_APAH14.heavy(:,k);

    cL = squeeze(sea_env_coef_APAH14(:,1,k));
    cM = squeeze(sea_env_coef_APAH14(:,2,k));
    cH = squeeze(sea_env_coef_APAH14(:,3,k));

    num_L = num_L + wL .* cL;  den_L = den_L + wL;
    num_M = num_M + wM .* cM;  den_M = den_M + wM;
    num_H = num_H + wH .* cH;  den_H = den_H + wH;
end

evap_eff_L = num_L ./ max(den_L, eps);
evap_eff_M = num_M ./ max(den_M, eps);
evap_eff_H = num_H ./ max(den_H, eps);

X(:,10) = evap_eff_L;
X(:,11) = evap_eff_M;
X(:,12) = evap_eff_H;

% 简单平均（解释更简单但物理意义弱），用这三行替换：
% X(:,10) = mean(sea_env_coef_APAH14(:,1,idx_APAH_nz), 3);
% X(:,11) = mean(sea_env_coef_APAH14(:,2,idx_APAH_nz), 3);
% X(:,12) = mean(sea_env_coef_APAH14(:,3,idx_APAH_nz), 3);


% ============
% 1×3 组图：左 SRC、右 PCC
figure; clf;

% (A) SRC

% 标准化
muX  = mean(X,1);
sigX = std(X,0,1);
muY  = mean(Y);
sigY = std(Y);

if any(sigX==0) || sigY==0
    error('SRC input has zero std (constant variable). Please check Σ8 APAHs inputs.');
end

for j = 1:p
    X(:,j) = (X(:,j) - muX(j)) ./ sigX(j);
end
Y = (Y - muY) ./ sigY;


Beta_boot = nan(B,p);
for b = 1:B
    idx = randi(N, N, 1);
    Beta_boot(b,:) = (X(idx,:) \ Y(idx))';
end

beta_med = median(Beta_boot,1);
ciL_SRC  = prctile(Beta_boot,2.5,1);
ciU_SRC  = prctile(Beta_boot,97.5,1);

subplot(1,3,1); hold on;
barh(1:p, beta_med, 0.6, 'FaceColor',[0.7 0.85 1]);
errorbar(beta_med, 1:p, beta_med-ciL_SRC, ciU_SRC-beta_med, ...
    'horizontal','LineStyle','none','Color','k','LineWidth',2,'CapSize',10);

set(gca, 'YTick',1:p, 'YTickLabel',varNames, 'YDir','reverse');
xlabel('Standardized Regression Coefficient');
title('SRC for \Sigma_{8} APAHs (non-zero)');
xline(0,'--k'); xlim([-0.1, 1.0]);
grid on; box on;

% ===  第3幅图 — SRC^2 贡献占比（与 Fig.5 同逻辑）===
beta_sq   = beta_med .^ 2;
share_raw = beta_sq ./ sum(beta_sq);

thresh = 0.01;                     % 1%
idxTop = find(share_raw >= thresh);
idxEtc = find(share_raw  < thresh);

share_plot  = [share_raw(idxTop),  sum(share_raw(idxEtc))];
labels_plot = [varNames(idxTop), {'Others (<1% each)'}];

subplot(1,3,3); hold on;
pie(share_plot, labels_plot);
title('SRC^2 share for \Sigma_{8} APAHs');
% 

%  MEM清理 SRC 临时大变量（避免影响 PCC）
clear Beta_boot idx muX sigX muY sigY
drawnow;

% -------
% (B) PCC

P_boot = nan(B,p);
for b = 1:B
    idx = randi(N, N, 1);
    P_boot(b,:) = partialcorri(Y(idx), X(idx,:));
end

P_med   = median(P_boot,1);
ciL_PCC = prctile(P_boot,2.5,1);
ciU_PCC = prctile(P_boot,97.5,1);

subplot(1,3,2); hold on;
barh(1:p, P_med, 0.6, 'FaceColor',[0.7 0.85 1.0]);
errorbar(P_med, 1:p, P_med-ciL_PCC, ciU_PCC-P_med, ...
    'horizontal','LineStyle','none','Color','k','LineWidth',2,'CapSize',10);

set(gca, 'YDir','reverse', 'YTick',1:p, 'YTickLabel',varNames);
xlabel('Partial Correlation Coefficient');
title('PCC for \Sigma_{8} APAHs (non-zero)');
xline(0,'--k'); xlim([-0.1, 1.0]);
grid on; box on;

% 输出检查表（可选）
fprintf('\n[OUTPUT] SRC table for Σ8 APAHs:\n');
disp(table(varNames', beta_med', ciL_SRC', ciU_SRC', ...
    'VariableNames',{'Variable','MedianSRC','CI2.5pct','CI97.5pct'}));

fprintf('\n[OUTPUT] PCC table for Σ8 APAHs:\n');
disp(table(varNames', P_med', ciL_PCC', ciU_PCC', ...
    'VariableNames',{'Variable','MedianPCC','CI2.5pct','CI97.5pct'}));

% MEM 清理 PCC 临时大变量
clear P_boot idx
drawnow;

toc;
