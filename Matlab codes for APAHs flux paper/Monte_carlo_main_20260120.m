%% 蒙特卡洛模拟：估算全球原油泄漏中 30 个 PAH/APAH 物种蒸发通量（kt/year）
% 按泄漏来源区分（自然、开采、运输、消费），再区分海洋与陆地，输出每个物种的总通量

% % —— 并行池设置 —— 
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     parpool('local',8);
% end
tic;
% 设置模拟次数
N = 1e6;

%% 步骤 1：蒙特卡洛采样泄漏原油量（单位：kt/year）
Oil_em_Nature_MC = trunc_lognrnd_mode(Oil_em_Nature,   Oil_em_Nature_min,   Oil_em_Nature_max,   N);
Oil_em_Ext_MC    = trunc_lognrnd_mode(Oil_em_Ext_y_DWH,Oil_em_Ext_y_DWH_min,Oil_em_Ext_y_DWH_max,N);
Oil_em_Trans_MC  = trunc_lognrnd_mode(Oil_em_Trans,    Oil_em_Trans_min,    Oil_em_Trans_max,    N);
Oil_em_Cons_MC   = trunc_lognrnd_mode(Oil_em_Cons,     Oil_em_Cons_min,     Oil_em_Cons_max,     N);

%% 步骤 2：采样轻中重原油比例
light_ratio  = normrnd(Oil_Rate_light,  Oil_Rate_light_err,  [N,1]);
medium_ratio = normrnd(Oil_Rate_medium, Oil_Rate_medium_err, [N,1]);
heavy_ratio  = normrnd(Oil_Rate_heavy,  Oil_Rate_heavy_err,  [N,1]);

light_ratio (light_ratio <0) = 0;
medium_ratio(medium_ratio<0) = 0;
heavy_ratio (heavy_ratio <0) = 0;
sum_ratio = light_ratio + medium_ratio + heavy_ratio;
light_ratio  = light_ratio  ./ sum_ratio;
medium_ratio = medium_ratio ./ sum_ratio;
heavy_ratio  = heavy_ratio  ./ sum_ratio;

%% 步骤 3：采样每类油的 30 物种浓度（单位已在导入时转为 ktPAH/ktOil）
% 逻辑保持：对每次模拟 i，分别在轻/中/重原油中抽取 1 个样本，作为当次三油型代表
nL = size(light_APAH14,  1);  % = 10
nM = size(medium_APAH14, 1);  % = 3
nH = size(heavy_APAH14,  1);  % = 5

% 预分配（存储每次模拟的代表样本浓度向量）
oil_APAH14.light  = zeros(N,14);
oil_APAH14.medium = zeros(N,14);
oil_APAH14.heavy  = zeros(N,14);

oil_PAH16.light   = zeros(N,16);
oil_PAH16.medium  = zeros(N,16);
oil_PAH16.heavy   = zeros(N,16);

for i = 1:N
    idx_L = randi(nL);   % 轻质原油样本索引
    idx_M = randi(nM);   % 中质原油样本索引
    idx_H = randi(nH);   % 重质原油样本索引

    % —— APAH14（14 物种）——
    oil_APAH14.light(i,:)  = light_APAH14 (idx_L,:);
    oil_APAH14.medium(i,:) = medium_APAH14(idx_M,:);
    oil_APAH14.heavy(i,:)  = heavy_APAH14 (idx_H,:);

    % —— PAH16（16 物种）——
    oil_PAH16.light(i,:)   = light_PAH16 (idx_L,:);
    oil_PAH16.medium(i,:)  = medium_PAH16(idx_M,:);
    oil_PAH16.heavy(i,:)   = heavy_PAH16 (idx_H,:);
end

%% 步骤 4：采样蒸发系数（海洋/陆地 × 轻/中/重），对每个物种独立应用
% 海洋环境系数 ~ 截断正态（mu±2σ），陆地系数 = 海洋系数×land_factor
sea_env_coef_APAH14 = zeros(N,3,14);   % [N, (light,med,heavy), species]
sea_env_coef_PAH16  = zeros(N,3,16);

land_factor  = normrnd(land_coef, land_coef_err, [N,1]);
land_factor(land_factor > 1) = 1;
land_factor(land_factor < 0) = 0;

% —— APAH14 —— 
for k = 1:14
    lb = max(light_mu - 2*light_sigma, 0); ub = light_mu  + 2*light_sigma;
    sea_env_coef_APAH14(:,1,k) = truncnormrnd(light_mu, light_sigma, lb, ub, N, 1);
    lb = max(med_mu   - 2*med_sigma,   0); ub = med_mu    + 2*med_sigma;
    sea_env_coef_APAH14(:,2,k) = truncnormrnd(med_mu,   med_sigma,   lb, ub, N, 1);
    lb = max(heavy_mu - 2*heavy_sigma, 0); ub = heavy_mu  + 2*heavy_sigma;
    sea_env_coef_APAH14(:,3,k) = truncnormrnd(heavy_mu, heavy_sigma, lb, ub, N, 1);
end

% —— PAH16 —— 
for k = 1:16
    lb = max(light_mu - 2*light_sigma, 0); ub = light_mu  + 2*light_sigma;
    sea_env_coef_PAH16(:,1,k) = truncnormrnd(light_mu, light_sigma, lb, ub, N, 1);
    lb = max(med_mu   - 2*med_sigma,   0); ub = med_mu    + 2*med_sigma;
    sea_env_coef_PAH16(:,2,k) = truncnormrnd(med_mu,   med_sigma,   lb, ub, N, 1);
    lb = max(heavy_mu - 2*heavy_sigma, 0); ub = heavy_mu  + 2*heavy_sigma;
    sea_env_coef_PAH16(:,3,k) = truncnormrnd(heavy_mu, heavy_sigma, lb, ub, N, 1);
end

%% 步骤 5：并行计算各源通量（单位：kt/year）
% 预分配（APAH14）
flux_nature_APAH14      = zeros(N,14);
flux_extraction_APAH14  = zeros(N,14);
flux_transport_APAH14   = zeros(N,14);
flux_consumption_APAH14 = zeros(N,14);

% 预分配（PAH16）
flux_nature_PAH16      = zeros(N,16);
flux_extraction_PAH16  = zeros(N,16);
flux_transport_PAH16   = zeros(N,16);
flux_consumption_PAH16 = zeros(N,16);

for i = 1:N
    % —— 原油量在三油型上的分配 —— 
    lr = light_ratio(i);  mr = medium_ratio(i);  hr = heavy_ratio(i);

    M_nat_light   = Oil_em_Nature_MC(i) * lr;
    M_nat_medium  = Oil_em_Nature_MC(i) * mr;
    M_nat_heavy   = Oil_em_Nature_MC(i) * hr;

    M_ext_light   = Oil_em_Ext_MC(i)    * lr;
    M_ext_medium  = Oil_em_Ext_MC(i)    * mr;
    M_ext_heavy   = Oil_em_Ext_MC(i)    * hr;

    M_trans_light  = Oil_em_Trans_MC(i) * lr;
    M_trans_medium = Oil_em_Trans_MC(i) * mr;
    M_trans_heavy  = Oil_em_Trans_MC(i) * hr;

    M_cons_light   = Oil_em_Cons_MC(i)  * lr;
    M_cons_medium  = Oil_em_Cons_MC(i)  * mr;
    M_cons_heavy   = Oil_em_Cons_MC(i)  * hr;

    % ===== APAH14：逐物种 =====
    for k = 1:14
        env_coef = squeeze(sea_env_coef_APAH14(i,:,k));  % [light,med,heavy]
        % 海洋最终蒸发比例
        fs = raw_evap_APAH14(k) * env_coef; fs(fs>1) = 1;
        % 陆地最终蒸发比例
        fl = raw_evap_APAH14(k) * (env_coef * land_factor(i));

        cl = oil_APAH14.light(i,k);
        cm = oil_APAH14.medium(i,k);
        ch = oil_APAH14.heavy(i,k);

        flux_nature_APAH14(i,k)      = M_nat_light   * cl * fs(1) + M_nat_medium  * cm * fs(2) + M_nat_heavy  * ch * fs(3);
        flux_extraction_APAH14(i,k)  = M_ext_light   * cl * fs(1) + M_ext_medium  * cm * fs(2) + M_ext_heavy  * ch * fs(3);
        flux_transport_APAH14(i,k)   = M_trans_light * cl * fs(1) + M_trans_medium* cm * fs(2) + M_trans_heavy* ch * fs(3);
        flux_consumption_APAH14(i,k) = M_cons_light  * cl * fl(1) + M_cons_medium * cm * fl(2) + M_cons_heavy * ch * fl(3);
    end

    % ===== PAH16：逐物种 =====
    for k = 1:16
        env_coef = squeeze(sea_env_coef_PAH16(i,:,k));   % [light,med,heavy]
        fs = raw_evap_PAH16(k) * env_coef; fs(fs>1) = 1;
        fl = raw_evap_PAH16(k) * (env_coef * land_factor(i));

        cl = oil_PAH16.light(i,k);
        cm = oil_PAH16.medium(i,k);
        ch = oil_PAH16.heavy(i,k);

        flux_nature_PAH16(i,k)      = M_nat_light   * cl * fs(1) + M_nat_medium  * cm * fs(2) + M_nat_heavy  * ch * fs(3);
        flux_extraction_PAH16(i,k)  = M_ext_light   * cl * fs(1) + M_ext_medium  * cm * fs(2) + M_ext_heavy  * ch * fs(3);
        flux_transport_PAH16(i,k)   = M_trans_light * cl * fs(1) + M_trans_medium* cm * fs(2) + M_trans_heavy* ch * fs(3);
        flux_consumption_PAH16(i,k) = M_cons_light  * cl * fl(1) + M_cons_medium * cm * fl(2) + M_cons_heavy * ch * fl(3);
    end
end

%% 汇总海洋源、陆地源、总通量
% APAH14
flux_marine_APAH14 = flux_nature_APAH14 + flux_extraction_APAH14 + flux_transport_APAH14;
flux_land_APAH14   = flux_consumption_APAH14;
flux_total_APAH14  = flux_marine_APAH14 + flux_land_APAH14;     % N×14

% PAH16
flux_marine_PAH16 = flux_nature_PAH16 + flux_extraction_PAH16 + flux_transport_PAH16;
flux_land_PAH16   = flux_consumption_PAH16;
flux_total_PAH16  = flux_marine_PAH16 + flux_land_PAH16;        % N×16

% 合并（可选）：N×30，物种顺序 = [APAH14, PAH16]
flux_total_all = [flux_total_APAH14, flux_total_PAH16];

% 名称列表（便于输出/绘图）
PAH30_names = [APAH14_names, PAH16_names];

%% —— 输出：每个物种的年蒸发通量统计（kt/year）——
for k = 1:numel(PAH30_names)
    data = flux_total_all(:,k);
    m_mean   = mean(data);
    m_median = prctile(data,50);
    ci95     = prctile(data,[2.5,97.5]);
    ext      = prctile(data,[0,100]);
    fprintf('%-6s 蒸发通量: 平均 = %.2e, 中位 = %.2e, 95%%CI = [%.2e, %.2e], 极值 = [%.2e, %.2e]\n', ...
        PAH30_names{k}, m_mean, m_median, ci95(1), ci95(2), ext(1), ext(2));
end


toc;

% 辅助函数：截断对数正态分布采样
function samples = trunc_lognrnd_mode(mu_mode, minv, maxv, N)
    % True truncated lognormal sampling with mode=mu_mode
    sigma = 0.5;
    m = log(mu_mode) ;  % mode constraint

    % Conditional CDF limits
    pLo = logncdf(minv, m, sigma);
    pHi = logncdf(maxv, m, sigma);

    % Sample uniformly in conditional probability space
    u = pLo + (pHi - pLo) .* rand(N,1);

    % Inverse CDF
    samples = logninv(u, m, sigma);
end
function x = truncnormrnd(mu, sigma, lb, ub, m, n)
    % 生成 m×n 的正态样本，然后截断到 [lb,ub]
    x = normrnd(mu, sigma, m, n);
    x(x<lb) = lb;
    x(x>ub) = ub;
end