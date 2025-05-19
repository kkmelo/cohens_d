# 混合模型效应量计算
# LMM建模残差，可以直接提取残差标准差，所以可以直接用Group*Time系数的交互作用/残差标准差来计算
estimate_interaction <- summary_model$coefficients["Group:Time", "Estimate"]# 提取模型交互作用未标准化系数
se_interaction <- summary_model$coefficients["Group:Time", "Std. Error"]# 提取模型标准误
residual_sd <- sigma(cleansbqonly_unadj)# 提取模型残差标准差
cohens_d_interaction <- estimate_interaction / residual_sd
se_d_interaction <- se_interaction / residual_sd# 提取cohens d的标准误，目的是将交互作用的标准误（se_interaction）除以残差标准差（residual_sd），目的是将标准误标准化到标准差单位，从而得到标准化效应量（如Cohen’s d）的标准误。
lower_ci_interaction <- cohens_d_interaction - 1.96 * se_d_interaction
upper_ci_interaction <- cohens_d_interaction + 1.96 * se_d_interaction
c(Cohen_s_d = cohens_d_interaction, Lower_95_CI = lower_ci_interaction, Upper_95_CI = upper_ci_interaction)
print(meanSBQR_total_y)
# 但GEE模型无法建模残差，所以无法直接提取残差标准差，参考过往文献只能用Group*Time系数的交互作用/基线时的 pooled SD 来计算
# 举例
GT_StigmaT =geeglm(StigmaT ~Group*Time+t0StigmaT+Age+Gender+Edu+Income+Religious+Workyear, data=GTdata, 
                   id=ID, corstr="ar1")
summary(GT_StigmaT)# 模型拟合
tidy(GT_StigmaT, conf.int = TRUE)# 计算未标准化系数(Estimate)的95CI
standardize_parameters(GT_StigmaT, method = "refit")# 计算标准化系数及95CI
# 例子为4个时间点包括GTdata$Time <- factor(GTdata$Time, levels = c(0,1,2,3), labels = c("T0_Pre", "T1_Post", "T2_6mo", "T3_12mo"))前测、后侧、6个月随访、12个月随访四个时间点
# 计算后测cohens d
interaction_estimate_sti1 <- coef(summary(GT_StigmaT))["GroupIntervention:TimeT1_Post", "Estimate"]# 提取Group*Time系数的交互作用未标准化系数
interaction_se_sti1 <- coef(summary(GT_StigmaT))["GroupIntervention:TimeT1_Post", "Std.err"]# 提取Group*Time系数的标准误
interaction_CI_lower <- interaction_estimate_sti1 - 1.96 * interaction_se_sti1 # 计算未标准化系数的置信区间
interaction_CI_upper <- interaction_estimate_sti1 + 1.96 * interaction_se_sti1 # 计算未标准化系数的置信区间
baseline_SD_sti <- sd(GTdata$StigmaT[GTdata$Time == "T0_Pre"], na.rm = TRUE)# 计算基线合并标准差
cohens_d_sti1 <- interaction_estimate_sti1 / baseline_SD_sti # 交互作用未标准化系数/合并标准差
cohens_d_CI_lower <- interaction_CI_lower / baseline_SD_sti # 把交互作用估计下界（interaction_CI_lower）标准化，得到Cohen’s d尺度的95%置信区间下界
cohens_d_CI_upper <- interaction_CI_upper / baseline_SD_sti # 把交互作用估计下界（interaction_CI_lower）标准化，得到Cohen’s d尺度的95%置信区间上界
cat("Cohen's d =", round(cohens_d_sti1, 3), "\n")# 保留3为小数
cat("95% CI = [", round(cohens_d_CI_lower, 3), ", ", round(cohens_d_CI_upper, 3), "]\n")
# 强制保留位数(可选)，format(..., nsmall = 4)：强制小数点后保留4位，不足补0。
# cat("Cohen's d =", format(round(cohens_d_com1, 4), nsmall = 4), "\n")
# cat("95% CI = [", 
    # format(round(cohens_d_CI_lower_com1, 4), nsmall = 4), ", ", 
    # format(round(cohens_d_CI_upper_com1, 4), nsmall = 4), "]\n")
# 计算6个月随访cohens d
interaction_estimate_sti1 <- coef(summary(GT_StigmaT))["GroupIntervention:TimeT2_6mo", "Estimate"]
interaction_se_sti1 <- coef(summary(GT_StigmaT))["GroupIntervention:TimeT2_6mo", "Std.err"]
interaction_CI_lower <- interaction_estimate_sti1 - 1.96 * interaction_se_sti1
interaction_CI_upper <- interaction_estimate_sti1 + 1.96 * interaction_se_sti1
baseline_SD_sti <- sd(GTdata$StigmaT[GTdata$Time == "T0_Pre"], na.rm = TRUE)
cohens_d_sti1 <- interaction_estimate_sti1 / baseline_SD_sti
cohens_d_CI_lower <- interaction_CI_lower / baseline_SD_sti
cohens_d_CI_upper <- interaction_CI_upper / baseline_SD_sti
cat("Cohen's d =", round(cohens_d_sti1, 3), "\n")
cat("95% CI = [", round(cohens_d_CI_lower, 3), ", ", round(cohens_d_CI_upper, 3), "]\n")
# 计算12个月随访cohens d
interaction_estimate_sti1 <- coef(summary(GT_StigmaT))["GroupIntervention:TimeT3_12mo", "Estimate"]
interaction_se_sti1 <- coef(summary(GT_StigmaT))["GroupIntervention:TimeT3_12mo", "Std.err"]
interaction_CI_lower <- interaction_estimate_sti1 - 1.96 * interaction_se_sti1
interaction_CI_upper <- interaction_estimate_sti1 + 1.96 * interaction_se_sti1
baseline_SD_sti <- sd(GTdata$StigmaT[GTdata$Time == "T0_Pre"], na.rm = TRUE)
cohens_d_sti1 <- interaction_estimate_sti1 / baseline_SD_sti
cohens_d_CI_lower <- interaction_CI_lower / baseline_SD_sti
cohens_d_CI_upper <- interaction_CI_upper / baseline_SD_sti
cat("Cohen's d =", round(cohens_d_sti1, 3), "\n")
cat("95% CI = [", round(cohens_d_CI_lower, 3), ", ", round(cohens_d_CI_upper, 3), "]\n")
# GEE计算cohens d方法参考文章：doi: 10.1037/a0030048
# 值得注意的是，文章中提到关于时间的编码中心化，假设实验有4个测量时点（基线、第1周、第2周、第3周），原始时间编码可能是0、1、2、3，求和除以时间点数量(n=4)做中心化，中心化后为-1.5、-0.5、0.5、1.5。
# 中心化原因：
 # 截距的解释：在回归模型（如公式8）中，截距 a 表示当时间变量取均值（即中心化后的0）时，因变量的预测值。此时，截距对应的是整个研究期间的“平均时间点”的效应，而非基线值（若用原始编码，截距仅代表基线时的效应）。
 # 交互项的解释：中心化后，组别×时间的交互项系数（b3）可直接反映两组斜率差异，避免因时间尺度非对称导致的偏差。
 # 减少共线性：若模型中包含高阶项（如时间平方项），中心化能降低变量间的相关性，提高模型稳定性。
# 但如果将Time因子化后，由于时间点不是相等间隔，可以直接未标准化系数/pooled SD(中心化后需要系数*时间点)。未标准化系数就分别代表，一天后的后测时、6个月的随访时、12个月的随访时分别的变化差异。
# 前提是一定要将Time因子化，如果不因子化会带来时间间隔不均匀问题（1天 vs. 6个月 vs. 6个月），回归系数 b3（交互项）会隐含假设“时间每增加1个单位，效应变化相同”，但实际1天 vs. 6个月 vs. 6个月。
# 若 b3=1，模型会错误地认为“干预组比对照组每天多1分”，而实际上1天后和6个月后的变化完全不同，由于时间尺度混乱，交互项系数无法合理反映“单位时间内的组间差异”。
# 因此，简约来讲，一定要将Time因子化

# LMM自编函数cohensd
extract_cohens_d <- function(model, digits = 3) {
  # 提取模型摘要
  summary_model <- summary(model)
  # 提取固定效应系数和标准误
  coef_table <- summary_model$coefficients
  # 提取残差标准差
  residual_sd <- sigma(model)
  # 提取 term 名称
  terms <- rownames(coef_table)
  # 计算 standardized effect sizes and CIs
  results <- data.frame(
    Term = terms,
    Estimate = coef_table[, "Estimate"],
    Std_Error = coef_table[, "Std. Error"],
    Cohen_d = coef_table[, "Estimate"] / residual_sd,
    SE_d = coef_table[, "Std. Error"] / residual_sd,
    stringsAsFactors = FALSE
  )
  # 添加置信区间
  results$Lower_95_CI <- results$Cohen_d - 1.96 * results$SE_d
  results$Upper_95_CI <- results$Cohen_d + 1.96 * results$SE_d
  # 四舍五入结果
  results <- results %>%
    mutate(across(where(is.numeric), ~ round(.x, digits)))
  return(results)
}








