###README###
### 1. Construct Linear Mix Model(LMM)
### 2.Test Normality for the data:skewness,Kolmogorov-Smirnov test
##-------------If normally distributed:Stick with LMM
##-------------If not, try to transform the data back to normally distributed(Optional step 1)
###Optional step 1:Transform data:log, square root
##-------------If it's normally distributed after transformation:Using LMM on transformed data
##-------------If not, use Gamma Generalized Linear model(GLMM)(Optional step 2)
###Optional step 2: Gamma GLMM(data is right skewed and contains only positive value). Test the fitness of the model.
### 3.Multiple-comparisons correction
### 4. Plotting

library(readxl)
library(lme4)
library(lmerTest)
library(tidyverse)
library(broom.mixed)
library(glmmTMB)
library(moments)
library(performance)


###1.DATA collection:Half width###
df <- read_excel("C:/Users/12868/OneDrive/Desktop/Netsted t test ndnf+ npy+.xlsx")

###2.Construct LMM model:HW###
NT_model_HW <- lmer(Half_width ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_HW)

###3.1 Check normality###
Half_width<-df$Half_width
hist(Half_width)
skewness(df$Half_width, na.rm = TRUE)
check_normality(NT_model_HW)

###3.2 If failed: transformation###
##option1(right_skewed): log##
df$log_HW <- log(df$Half_width)
NT_model_HW_log <- lmer(log_HW ~ Cell_type + (1 | Animal_ID), data = df)
check_normality(NT_model_HW_log)
summary(NT_model_HW_log)

posthoc_HW <-emmeans(NT_model_HW_log, list(pairwise ~ Cell_type), mode = "satterthwaite")
posthoc_HW

edf_value_HW <- df.residual(NT_model_HW_log)
edf_size_HW <- eff_size(posthoc_HW, sigma = sigma(NT_model_HW_log), edf = edf_value_HW)
edf_size_HW

###RT###
NT_model_RT <- lmer(Rise_tau ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_RT)

###3.1 Check normality###
Rise_tau<-df$Rise_tau
hist(Rise_tau)
skewness(df$Rise_tau, na.rm = TRUE)
check_normality(NT_model_RT)

###3.2 If failed: transformation###
##option1(right_skewed): log##
log_RT <- log(df$Rise_tau)
NT_model_RT_log <- lmer(log_RT ~ Cell_type + (1 | Animal_ID), data = df)
check_normality(NT_model_RT_log)
summary(NT_model_RT_log)
skewness(log_RT)
##option2: sqrt##
NT_model_RT_sqrt <- lmer(sqrt(Rise_tau) ~ Cell_type + (1 | Animal_ID), data = df)
check_normality(NT_model_RT_sqrt)

###3.2 If failed: transformation###
model_glmm_RT <- glmmTMB(Rise_tau ~ Cell_type + (1 | Animal_ID), family = Gamma(link = "log"), data = df)
summary(model_glmm_RT)

posthoc_RT <- emmeans(model_glmm_RT, list(pairwise ~ Cell_type), mode = "satterthwaite")
posthoc_RT
##effect_size##
edf_value_RT <- df.residual(model_glmm_RT)
edf_size_RT<-eff_size(posthoc_RT$emmeans, sigma = sigma(model_glmm_RT), edf = edf_value_RT)
edf_size_RT

##DT##

NT_model_DT <- lmer(Decay_tau ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_DT)

###3.1 Check normality###
Decay_tau<-df$Decay_tau
hist(Decay_tau)
skewness(df$Decay_tau, na.rm = TRUE)
check_normality(NT_model_DT)

###3.2 If failed: transformation###
##option1(first option): log##
log_DT <- log(df$Decay_tau)
NT_model_DT_log <- lmer(log_DT ~ Cell_type + (1 | Animal_ID), data = df)
check_normality(NT_model_DT_log)
summary(NT_model_DT_log)
skewness(log_DT)

posthoc_DT <- emmeans(NT_model_DT_log, list(pairwise ~ Cell_type), mode = "satterthwaite")

posthoc_DT
##effect_size##
edf_value_DT <- df.residual(NT_model_DT_log)

edf_size_DT<-eff_size(posthoc_DT$emmeans, sigma = sigma(NT_model_DT_log), edf = edf_value_DT)
edf_size_DT

##A###
NT_model_A <- lmer(Area ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_A)

###3.1 Check normality###
Area <-df$Area
hist(Area)
skewness(df$Area, na.rm = TRUE)
check_normality(NT_model_A)
posthoc_A <- emmeans(NT_model_A, list(pairwise ~ Cell_type), mode = "satterthwaite")

posthoc_A
##effect_size##
edf_value_A <- df.residual(NT_model_A)

edf_size_A<-eff_size(posthoc_A$emmeans, sigma = sigma(NT_model_A), edf = edf_value_A)
edf_size_A

##PA###
NT_model_PA <- lmer(Peak_amplitude ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_PA)

###3.1 Check normality###
Peak_amplitude <-df$Peak_amplitude
hist(Peak_amplitude)
skewness(df$Peak_amplitude, na.rm = TRUE)
check_normality(NT_model_PA)
posthoc_PA <- emmeans(NT_model_PA, list(pairwise ~ Cell_type), mode = "satterthwaite")

posthoc_PA
##effect_size##
edf_value_PA <- df.residual(NT_model_PA)

edf_size_PA<-eff_size(posthoc_PA$emmeans, sigma = sigma(NT_model_PA), edf = edf_value_PA)
edf_size_PA

##MRS###
NT_model_MRS <- lmer(Max_rise_slope ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_MRS)

###3.1 Check normality###
Max_rise_slope <-df$Max_rise_slope
hist(Max_rise_slope)
skewness(df$Max_rise_slope, na.rm = TRUE)
check_normality(NT_model_MRS)
posthoc_MRS <- emmeans(NT_model_MRS, list(pairwise ~ Cell_type), mode = "satterthwaite")

posthoc_MRS

##effect_size##
edf_value_MRS <- df.residual(NT_model_MRS)

edf_size_MRS <-eff_size(posthoc_MRS$emmeans, sigma = sigma(NT_model_MRS), edf = edf_value_MRS)

edf_size_MRS

##MDS###

NT_model_MDS <- lmer(Max_decay_slope ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_MDS)

###3.1 Check normality###
Max_decay_slope <-df$Max_decay_slope
hist(Max_decay_slope)
skewness(df$Max_decay_slope, na.rm = TRUE)
check_normality(NT_model_MDS)

posthoc_MDS <- emmeans(NT_model_MDS, list(pairwise ~ Cell_type), mode = "satterthwaite")

posthoc_MDS


##effect_size##
edf_value_MDS <- df.residual(NT_model_MDS)
edf_size_MDS <-eff_size(posthoc_MDS$emmeans, sigma = sigma(NT_model_MDS), edf = edf_value_MDS)
edf_size_MDS



##Rheobase###
NT_model_RB <- lmer(Rheobase ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_RB)

###3.1 Check normality###
Rheobase <-df$Rheobase
hist(Rheobase)
skewness(df$Rheobase, na.rm = TRUE)
check_normality(NT_model_RB)

###3.2 If failed: transformation###
##option1(first option): log##
log_RB <- log(Rheobase)
NT_model_RB_log <- lmer(log_RB ~ Cell_type + (1 | Animal_ID), data = df)
check_normality(NT_model_RB_log)
summary(NT_model_RB_log)
skewness(log_RB)

posthoc_RB <- emmeans(NT_model_RB_log, list(pairwise ~ Cell_type), mode = "satterthwaite")
posthoc_RB


##effect_size##
edf_value_RB <- df.residual(NT_model_RB_log)
edf_size_RB <-eff_size(posthoc_RB$emmeans, sigma = sigma(NT_model_RB), edf = edf_value_RB)
edf_size_RB




##L###
NT_model_L <- lmer(Latency ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_L)

###3.1 Check normality###
Latency <-df$Latency
hist(Latency)
skewness(df$Latency, na.rm = TRUE)
check_normality(NT_model_L)

###3.2 If failed: transformation###
##option1(first option): log##
log_L <- log(Latency)
NT_model_L_log <- lmer(log_L ~ Cell_type + (1 | Animal_ID), data = df)
check_normality(NT_model_L_log)
summary(NT_model_L_log)
skewness(log_L)

NT_model_L_sqrt <- lmer(sqrt(Latency) ~ Cell_type + (1 | Animal_ID), data = df)
check_normality(NT_model_L_sqrt)

library(glmmTMB)
model_glmm_L <- glmmTMB(Latency ~ Cell_type + (1 | Animal_ID), 
                        family = Gamma(link = "log"), 
                        data = df)
summary(model_glmm_L)

posthoc_L <- emmeans(model_glmm_L, list(pairwise ~ Cell_type), mode = "satterthwaite")
posthoc_L

##effect_size##
edf_value_L <- df.residual(model_glmm_L)
edf_size_L <-eff_size(posthoc_L$emmeans, sigma = sigma(model_glmm_L), edf = edf_value_L)
edf_size_L


##RMP###
NT_model_RMP <- lmer(RMP ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_RMP)


###3.1 Check normality###
RMP <-df$RMP
hist(RMP)
skewness(df$RMP, na.rm = TRUE)
check_normality(NT_model_RMP)

posthoc_RMP <- emmeans(NT_model_RMP, list(pairwise ~ Cell_type), mode = "satterthwaite")
posthoc_RMP

##effect_size##
edf_value_RMP <- df.residual(NT_model_RMP)
edf_size_RMP <-eff_size(posthoc_RMP$emmeans, sigma = sigma(NT_model_RMP), edf = edf_value_RMP)
edf_size_RMP

##R###
NT_model_R <- lmer(Resistance ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_R)

###3.1 Check normality###
R <-df$Resistance
hist(R)
skewness(df$Resistance, na.rm = TRUE)
check_normality(NT_model_R)
posthoc_R <- emmeans(NT_model_R, list(pairwise ~ Cell_type), mode = "satterthwaite")
posthoc_R

##effect_size##
edf_value_R <- df.residual(NT_model_R)
edf_size_R <-eff_size(posthoc_R$emmeans, sigma = sigma(NT_model_R), edf = edf_value_R)
edf_size_R

##Tau###
NT_model_tau <- lmer(tau ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_tau)

###3.1 Check normality###
tau<-df$tau
hist(tau)
skewness(df$tau, na.rm = TRUE)
check_normality(NT_model_tau)
posthoc_tau <- emmeans(NT_model_tau, list(pairwise ~ Cell_type), mode = "satterthwaite")
posthoc_tau

##effect_size##
edf_value_tau <- df.residual(NT_model_tau)
edf_size_tau <-eff_size(posthoc_tau$emmeans, sigma = sigma(NT_model_tau), edf = edf_value_tau)
edf_size_tau

##2.Construct LMM model###

NT_model_C <- lmer(Compacitance ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_C)

###3.1 Check normality###
C<-df$Compacitance
hist(C)
skewness(df$Compacitance, na.rm = TRUE)
check_normality(NT_model_C)

posthoc_C <- emmeans(NT_model_C, list(pairwise ~ Cell_type), mode = "satterthwaite")
posthoc_C

##effect_size##
edf_value_C <- df.residual(NT_model_C)
edf_size_C <-eff_size(posthoc_C$emmeans, sigma = sigma(NT_model_C), edf = edf_value_C)
edf_size_C



##SR###
NT_model_SR <- lmer(Sag_ratio ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_SR)

###3.1 Check normality###
SR<-df$Sag_ratio
hist(SR)
skewness(df$Sag_ratio, na.rm = TRUE)
check_normality(NT_model_SR)

log_SR<-log(SR)
NT_model_SR_log <- lmer(log(Sag_ratio) ~ Cell_type + (1 | Animal_ID), data = df)
check_normality(NT_model_SR_log)
summary(NT_model_SR_log)
skewness(log_SR)

posthoc_SR <- emmeans(NT_model_SR_log, list(pairwise ~ Cell_type), mode = "satterthwaite")
posthoc_SR

##effect_size##
edf_value_SR <- df.residual(NT_model_SR)
edf_size_SR <-eff_size(posthoc_SR$emmeans, sigma = sigma(NT_model_SR), edf = edf_value_SR)
edf_size_SR

###F###

NT_model_F <- lmer(Frequency ~ Cell_type + (1 | Animal_ID), data = df)
summary(NT_model_F)

###3.1 Check normality###
F<-df$Frequency 
hist(F)
skewness(df$Frequency, na.rm = TRUE)
check_normality(NT_model_F)

###3.2 If failed: transformation###
##option1(first option): log##
log_F<-log(F)
NT_model_F_log <- lmer(Frequency ~ Cell_type + (1 | Animal_ID), data = df)
check_normality(NT_model_F_log)
summary(NT_model_F_log)
skewness(log_F)

posthoc_F <- emmeans(NT_model_F_log, list(pairwise ~ Cell_type), mode = "satterthwaite")
posthoc_F
##effect_size##
edf_value_F<- df.residual(NT_model_F_log)
edf_size_F<-eff_size(posthoc_F$emmeans, sigma = sigma(NT_model_F_log), edf = edf_value_F)
edf_size_F

###Multiple-comparison correction###
df <- read_excel("C:/Users/12868/OneDrive/Desktop/Leena/IGL project/Statistical_analysis/Excel/LMM_test_results_1.xlsx")


df$P_adj_BH <- p.adjust(df$P_value, method = "BH")
df$P_adj_Bonferroni <- p.adjust(df$P_value, method = "bonferroni")

View(df)

significant_features_BH <- df[df$P_adj_BH < 0.05, ]
significant_features_BH
significant_features_Bonferroni <- df[df$P_adj_Bonferroni  < 0.05, ]

###Plotting###

library(ggplot2)
library(ggpubr)
library(patchwork)  


ylabels <- c(
  Rise_tau       = "Rise tau (ms)",
  Decay_tau      = "Decay tau (ms)",
  Half_width     = "Half width (ms)",
  Area           = "Area (mV·ms)",
  Peak_amplitude = "Peak amplitude (mV)",
  Max_rise_slope = "Max rise slope (mV/ms)",
  Max_decay_slope= "Max decay slope (mV/ms)",
  Rheobase       = "Rheobase (pA)",
  Latency        = "Latency (ms)",
  Frequency      = "Frequency (Hz)",
  RMP            = "RMP (mV)",
  Resistance     = "Resistance (MΩ)",
  tau            = "Tau (ms)",
  Compacitance   = "Capacitance (pF)",
  Sag_ratio      = "Sag ratio"
)
##Significance##
res <- read_excel("C:/Users/12868/OneDrive/Desktop/Leena/IGL project/Statistical_analysis/Excel/LMM_test_results_1.xlsx")

get_sig_label <- function(p) {
  if (is.na(p)) return("")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  return("ns")
}
res$P_adj_FDR
res$Sig_Label <- sapply(res$P_adj_FDR, get_sig_label)

##Drawing loop##
plot_list <- list()

for (feat in features) {
  
current_sig <- results$Sig_FDR[results$Feature == feat]
y_pos <- max(df[[feat]], na.rm = TRUE) * 1.15 
  
p <- ggplot(df, aes(x = Cell_type, y = .data[[feat]], fill = Cell_type)) +
    
geom_violin(trim = FALSE, color = "grey30", linewidth = 0.2, alpha = 0.6) +
    
geom_jitter(width = 0.1, size = 1, alpha = 0.4, color = "grey10") +
    
annotate("text", x = 1.5, 
  y = max(df[[feat]], na.rm = TRUE) * 1.2, 
  label = current_sig, 
  size = 5, fontface = "bold", color = "grey20")+
  
scale_fill_manual(values = c("NPY+" = "#C8E6C9", "NDNF+" = "#D1C4E9")) +
    
labs(x = NULL, y = ylabels[feat]) +
  
theme_classic(base_size = 10) +
theme(
      legend.position = "none",
      axis.text       = element_text(color = "black",face = "bold"), 
      axis.title.y    = element_text(size = 11, color = "black",face = "bold"),
      axis.line       = element_line(color = "grey40", linewidth = 0.3),
      axis.ticks      = element_line(color = "grey40", linewidth = 0.3),
      plot.margin     = unit(c(5, 5, 5, 5), "pt") 
    )
  
  plot_list[[feat]] <- p
}

#combine figure
combined <- wrap_plots(plot_list, ncol = 4) & 
  theme(plot.margin = margin(10, 10, 10, 10)) 

print(combined)

