%    - 功能：读取并初始化所有输入数据，包括油泄漏量、原油类型比例、APAHs浓度、蒸发率参数等
     clearvars; clc;


% 添加Shen et al., 2021年排放清单数据（单位：kt/year）
MeNAP_GEMS_Glabal = 69.483;
PAH40_GEMS_Glabal = 678.03;
PAH16_GEMS_Glabal = 557.085;
PAH8_GEMS_Glabal = 455.715;

% Step 1: Global Spill Volume Distribution.
% 2010–2019 data from Oil in the sea 4

data_oil_emission = readmatrix('APAHs15 & oil calculated 3.0 20251023.xlsx', ...
    'Sheet', 'Oil emission', 'Range', 'K25:M29');

Oil_em_Nature           = data_oil_emission(1,1);
Oil_em_Nature_min       = data_oil_emission(1,2);
Oil_em_Nature_max       = data_oil_emission(1,3);
Oil_em_Ext_n_DWH        = data_oil_emission(2,1);
Oil_em_Ext_n_DWH_min    = data_oil_emission(2,2);
Oil_em_Ext_n_DWH_max    = data_oil_emission(2,3);
Oil_em_Ext_y_DWH        = data_oil_emission(3,1);
Oil_em_Ext_y_DWH_min    = data_oil_emission(3,2);
Oil_em_Ext_y_DWH_max    = data_oil_emission(3,3);
Oil_em_Trans            = data_oil_emission(4,1);
Oil_em_Trans_min        = data_oil_emission(4,2);
Oil_em_Trans_max        = data_oil_emission(4,3);
Oil_em_Cons             = data_oil_emission(5,1);
Oil_em_Cons_min         = data_oil_emission(5,2);
Oil_em_Cons_max         = data_oil_emission(5,3);

% Step 2: Oil Type Classification and Allocation.
% 依据美国2014-2024年炼油的API分类，重<22.3, 中22.3-31.1，轻>31.1
% rate 是轻中重油相对于开采总量的占比，err后缀表示一个标准差；
data_oil_type = readmatrix('APAHs15 & oil calculated 3.0 20251023.xlsx', ...
    'Sheet', 'light_mid_heavy', 'Range', 'P114:R116');

Oil_Rate_heavy          = data_oil_type(1,1);
Oil_Rate_heavy_err      = data_oil_type(2,1);
Oil_Rate_medium         = data_oil_type(1,2);
Oil_Rate_medium_err     = data_oil_type(2,2);
Oil_Rate_light          = data_oil_type(1,3);
Oil_Rate_light_err      = data_oil_type(2,3);

% Step 3: Chemical Composition – APAH14 & PAH16 (浓度: μg/g)

% === 读入 14 个 APAH 组分 ===
data_APAH14_oil = readmatrix('APAHs15 & oil calculated 3.0 20251023.xlsx', ...
    'Sheet','APAHs15 in Oil','Range','F79:S96');   % 14 列

% === 读入 16 个 PAH16 组分 ===
data_PAH16_oil = readmatrix('APAHs15 & oil calculated 3.0 20251023.xlsx', ...
    'Sheet','APAHs15 in Oil','Range','F130:U147');  % 16 列

% —— 切分三类原油（行）——
light_APAH14  = data_APAH14_oil(1:10,:);     % 轻质
medium_APAH14 = data_APAH14_oil(11:13,:);
heavy_APAH14  = data_APAH14_oil(14:18,:);

light_PAH16   = data_PAH16_oil(1:10,:);
medium_PAH16  = data_PAH16_oil(11:13,:);
heavy_PAH16   = data_PAH16_oil(14:18,:);

% —— 单位统一：μg/g = mg/kg = (kt PAH / kt oil) ——
light_APAH14  = light_APAH14  * 1e-6;
medium_APAH14 = medium_APAH14 * 1e-6;
heavy_APAH14  = heavy_APAH14  * 1e-6;

light_PAH16   = light_PAH16   * 1e-6;
medium_PAH16  = medium_PAH16  * 1e-6;
heavy_PAH16   = heavy_PAH16   * 1e-6;

% —— 列映射（APAH14 共 14 列）——
% 1 C1NAP, 2 C2NAP, 3 C3NAP, 4 C4NAP,
% 5 C1PHE, 6 C2PHE, 7 C3PHE, 8 C4PHE,
% 9 C1FLO,10 C2FLO,11 C3FLO,
% 12 C1CHR,13 C2CHR,14 C3CHR

% 轻质（APAH14）
C1NAP_light = light_APAH14(:,1);  C2NAP_light = light_APAH14(:,2);
C3NAP_light = light_APAH14(:,3);  C4NAP_light = light_APAH14(:,4);
C1PHE_light = light_APAH14(:,5);  C2PHE_light = light_APAH14(:,6);
C3PHE_light = light_APAH14(:,7);  C4PHE_light = light_APAH14(:,8);
C1FLO_light = light_APAH14(:,9);  C2FLO_light = light_APAH14(:,10);
C3FLO_light = light_APAH14(:,11);
C1CHR_light = light_APAH14(:,12); C2CHR_light = light_APAH14(:,13);
C3CHR_light = light_APAH14(:,14);

% 中质（APAH14）
C1NAP_med = medium_APAH14(:,1);  C2NAP_med = medium_APAH14(:,2);
C3NAP_med = medium_APAH14(:,3);  C4NAP_med = medium_APAH14(:,4);
C1PHE_med = medium_APAH14(:,5);  C2PHE_med = medium_APAH14(:,6);
C3PHE_med = medium_APAH14(:,7);  C4PHE_med = medium_APAH14(:,8);
C1FLO_med = medium_APAH14(:,9);  C2FLO_med = medium_APAH14(:,10);
C3FLO_med = medium_APAH14(:,11);
C1CHR_med = medium_APAH14(:,12); C2CHR_med = medium_APAH14(:,13);
C3CHR_med = medium_APAH14(:,14);

% 重质（APAH14）
C1NAP_heavy = heavy_APAH14(:,1);  C2NAP_heavy = heavy_APAH14(:,2);
C3NAP_heavy = heavy_APAH14(:,3);  C4NAP_heavy = heavy_APAH14(:,4);
C1PHE_heavy = heavy_APAH14(:,5);  C2PHE_heavy = heavy_APAH14(:,6);
C3PHE_heavy = heavy_APAH14(:,7);  C4PHE_heavy = heavy_APAH14(:,8);
C1FLO_heavy = heavy_APAH14(:,9);  C2FLO_heavy = heavy_APAH14(:,10);
C3FLO_heavy = heavy_APAH14(:,11);
C1CHR_heavy = heavy_APAH14(:,12); C2CHR_heavy = heavy_APAH14(:,13);
C3CHR_heavy = heavy_APAH14(:,14);

% —— APAH14 汇总（每行求和）——
APAH14_light = sum(light_APAH14, 2);
APAH14_med   = sum(medium_APAH14,2);
APAH14_heavy = sum(heavy_APAH14, 2);

% —— 列映射（PAH16 共 16 列）——
% 1 NAP, 2 PHE, 3 FLO, 4 CHR, 5 ACY, 6 ACE, 7 ANT, 8 FLA,
% 9 PYR, 10 BaA, 11 BbF, 12 BkF, 13 BaP, 14 DahA, 15 IcdP, 16 BghiP

% 轻质（PAH16）
NAP_light  = light_PAH16(:,1);   PHE_light  = light_PAH16(:,2);
FLO_light  = light_PAH16(:,3);   CHR_light  = light_PAH16(:,4);
ACY_light  = light_PAH16(:,5);   ACE_light  = light_PAH16(:,6);
ANT_light  = light_PAH16(:,7);   FLA_light  = light_PAH16(:,8);
PYR_light  = light_PAH16(:,9);   BaA_light  = light_PAH16(:,10);
BbF_light  = light_PAH16(:,11);  BkF_light  = light_PAH16(:,12);
BaP_light  = light_PAH16(:,13);  DahA_light = light_PAH16(:,14);
IcdP_light = light_PAH16(:,15);  BghiP_light= light_PAH16(:,16);

% 中质（PAH16）
NAP_med  = medium_PAH16(:,1);   PHE_med  = medium_PAH16(:,2);
FLO_med  = medium_PAH16(:,3);   CHR_med  = medium_PAH16(:,4);
ACY_med  = medium_PAH16(:,5);   ACE_med  = medium_PAH16(:,6);
ANT_med  = medium_PAH16(:,7);   FLA_med  = medium_PAH16(:,8);
PYR_med  = medium_PAH16(:,9);   BaA_med  = medium_PAH16(:,10);
BbF_med  = medium_PAH16(:,11);  BkF_med  = medium_PAH16(:,12);
BaP_med  = medium_PAH16(:,13);  DahA_med = medium_PAH16(:,14);
IcdP_med = medium_PAH16(:,15);  BghiP_med= medium_PAH16(:,16);

% 重质（PAH16）
NAP_heavy  = heavy_PAH16(:,1);   PHE_heavy  = heavy_PAH16(:,2);
FLO_heavy  = heavy_PAH16(:,3);   CHR_heavy  = heavy_PAH16(:,4);
ACY_heavy  = heavy_PAH16(:,5);   ACE_heavy  = heavy_PAH16(:,6);
ANT_heavy  = heavy_PAH16(:,7);   FLA_heavy  = heavy_PAH16(:,8);
PYR_heavy  = heavy_PAH16(:,9);   BaA_heavy  = heavy_PAH16(:,10);
BbF_heavy  = heavy_PAH16(:,11);  BkF_heavy  = heavy_PAH16(:,12);
BaP_heavy  = heavy_PAH16(:,13);  DahA_heavy = heavy_PAH16(:,14);
IcdP_heavy = heavy_PAH16(:,15);  BghiP_heavy= heavy_PAH16(:,16);


% Step 4: Evaporation Fraction (30 species: 14 APAH14 + 16 PAH16)
% 轻质油 MC252 的纯物理蒸发比例（Yin, 2015）
raw = readmatrix('APAHs15 & oil calculated 3.0 20251023.xlsx', ...
    'Sheet','Evaporation rate','Range','G12:G41');

% 保证比例在 [0,1]
raw(raw < 0) = 0;
raw(raw > 1) = 1;

% === 分组：APAH14(1–14) + PAH16(15–30) ===
raw_evap_APAH14 = raw(1:14).';   % 1×14
raw_evap_PAH16  = raw(15:30).';  % 1×16

% —— APAH14 名称与索引（顺序必须与数据一致）——
APAH14_names = {'C1NAP','C2NAP','C3NAP','C4NAP', ...
                'C1PHE','C2PHE','C3PHE','C4PHE', ...
                'C1FLO','C2FLO','C3FLO','C1CHR','C2CHR','C3CHR'};
idx_APAH14 = cell2struct(num2cell(1:14), APAH14_names, 2);

% —— PAH16 名称与索引（顺序必须与数据一致）——
PAH16_names = {'NAP','PHE','FLO','CHR','ACY','ACE','ANT','FLA', ...
               'PYR','BaA','BbF','BkF','BaP','DahA','IcdP','BghiP'};
idx_PAH16 = cell2struct(num2cell(1:16), PAH16_names, 2);

% —— 单个物种的独立变量（便于单独画图/检查）——
% APAH14:
C1NAP_Evapor = raw_evap_APAH14(1);  C2NAP_Evapor = raw_evap_APAH14(2);
C3NAP_Evapor = raw_evap_APAH14(3);  C4NAP_Evapor = raw_evap_APAH14(4);
C1PHE_Evapor = raw_evap_APAH14(5);  C2PHE_Evapor = raw_evap_APAH14(6);
C3PHE_Evapor = raw_evap_APAH14(7);  C4PHE_Evapor = raw_evap_APAH14(8);
C1FLO_Evapor = raw_evap_APAH14(9);  C2FLO_Evapor = raw_evap_APAH14(10);
C3FLO_Evapor = raw_evap_APAH14(11);
C1CHR_Evapor = raw_evap_APAH14(12); C2CHR_Evapor = raw_evap_APAH14(13);
C3CHR_Evapor = raw_evap_APAH14(14);

% PAH16:
NAP_Evapor   = raw_evap_PAH16(1);   PHE_Evapor   = raw_evap_PAH16(2);
FLO_Evapor   = raw_evap_PAH16(3);   CHR_Evapor   = raw_evap_PAH16(4);
ACY_Evapor   = raw_evap_PAH16(5);   ACE_Evapor   = raw_evap_PAH16(6);
ANT_Evapor   = raw_evap_PAH16(7);   FLA_Evapor   = raw_evap_PAH16(8);
PYR_Evapor   = raw_evap_PAH16(9);   BaA_Evapor   = raw_evap_PAH16(10);
BbF_Evapor   = raw_evap_PAH16(11);  BkF_Evapor   = raw_evap_PAH16(12);
BaP_Evapor   = raw_evap_PAH16(13);  DahA_Evapor  = raw_evap_PAH16(14);
IcdP_Evapor  = raw_evap_PAH16(15);  BghiP_Evapor = raw_evap_PAH16(16);


% —— 读取海洋修正系数的均值和标准差 —— 
% 从 Correction factors 表格 E113:J113 读取 1×6 向量：
% [light_mean, light_std, med_mean, med_std, heavy_mean, heavy_std]
sea_evap_params = readmatrix('APAHs15 & oil calculated 3.0 20251023.xlsx', ...
    'Sheet','Correction factors','Range','E113:J113');  % 1×6

% 方便单独引用
light_mu    = sea_evap_params(1);
light_sigma = sea_evap_params(2);
med_mu      = sea_evap_params(3);
med_sigma   = sea_evap_params(4);
heavy_mu    = sea_evap_params(5);
heavy_sigma = sea_evap_params(6);

% 陆地浮油蒸发系数
land_coef = 0.734;
land_coef_err = 0.026;
