###Normal T test###

library(readxl)
library(lme4)
library(lmerTest)
library(tidyverse)
library(broom.mixed)

### 1. DATA collection
df <- read_excel("C:/Users/12868/OneDrive/Desktop/Netsted t test ndnf+ npy+.xlsx")


###Rise_tau 
###2. Normality check
shapiro.test(df$Rise_tau )

###3. Mann-Whitney U test
RT_t_test<- wilcox.test(Rise_tau ~ Cell_type, data = df, 
                        exact = FALSE,      
                        conf.int = TRUE)
print(RT_t_test)

###Decay_tau 
###2. Normality check
shapiro.test(df$Decay_tau)

###3. Mann-Whitney U test
DT_t_test<- wilcox.test(Decay_tau ~ Cell_type, data = df, 
                        exact = FALSE,      
                        conf.int = TRUE)
print(DT_t_test)


###Half width 
###2. Normality check
shapiro.test(df$Half_width)

###3. Mann-Whitney U test
HW_t_test<- wilcox.test(Half_width ~ Cell_type, data = df, 
                        exact = FALSE,      
                        conf.int = TRUE)
print(HW_t_test)

###Area
###2. Normality check
shapiro.test(df$Area)

###3. Welch test
A_t_test <- t.test(Area ~ Cell_type, data = df)
print(A_t_test)

###Peak_amplitude
###2. Normality check
shapiro.test(df$Peak_amplitude)

###3. Mann-Whitney U test
PA_t_test<- wilcox.test(Peak_amplitude ~ Cell_type, data = df, 
                        exact = FALSE,      
                        conf.int = TRUE)
print(PA_t_test)


###Max_rise_slope
###2. Normality check
shapiro.test(df$Max_rise_slope)

###3. Welch test
MRS_t_test <- t.test(Max_rise_slope ~ Cell_type, data = df)
print(MRS_t_test)


###Max_decay_slope
###2. Normality check
shapiro.test(df$Max_decay_slope)

###3. Mann-Whitney U test
MDS_t_test<- wilcox.test(Max_decay_slope ~ Cell_type, data = df, 
                        exact = FALSE,      
                        conf.int = TRUE)
print(MDS_t_test)


###Rheobase
###2. Normality check
shapiro.test(df$Rheobase)

###3. Mann-Whitney U test
R_t_test<- wilcox.test(Rheobase ~ Cell_type, data = df, 
                         exact = FALSE,      
                         conf.int = TRUE)
print(R_t_test)

###Latency
###2. Normality check
shapiro.test(df$Latency)

###3. Mann-Whitney U test
L_t_test<- wilcox.test(Latency ~ Cell_type, data = df, 
                       exact = FALSE,      
                       conf.int = TRUE)
print(L_t_test)


###Frequency
###2. Normality check
shapiro.test(df$Frequency)

###3. Mann-Whitney U test
F_t_test<- wilcox.test(Frequency ~ Cell_type, data = df, 
                       exact = FALSE,      
                       conf.int = TRUE)
print(F_t_test)


###RMP
###2. Normality check
shapiro.test(df$RMP)

###3. Welch t test
RMP_t_test <- t.test(RMP ~ Cell_type, data = df)
print(RMP_t_test)


###Resistance
###2. Normality check
shapiro.test(df$Resistance)

###3. Mann-Whitney U test
R_t_test<- wilcox.test(Resistance ~ Cell_type, data = df, 
                       exact = FALSE,      
                       conf.int = TRUE)
print(R_t_test)


###tau
###2. Normality check
shapiro.test(df$tau)

###3. Welch t test
T_t_test <- t.test(tau ~ Cell_type, data = df)
print(T_t_test)

###Capacitance
###2. Normality check
shapiro.test(df$Compacitance)

###3. Mann-Whitney U test
C_t_test<- wilcox.test(Compacitance ~ Cell_type, data = df, 
                       exact = FALSE,      
                       conf.int = TRUE)
print(C_t_test)

###Sag_ratio
###2. Normality check
shapiro.test(df$Sag_ratio)

###3. Mann-Whitney U test
SR_t_test<- wilcox.test(Sag_ratio ~ Cell_type, data = df, 
                       exact = FALSE,      
                       conf.int = TRUE)
print(SR_t_test)
###Multiple-comparison correction###
df <- read_excel("C:/Users/12868/OneDrive/Desktop/Statistical_Results.xlsx")

df$P_adj_BH <- p.adjust(df$P_value, method = "BH")
df$P_adj_Bonferroni <- p.adjust(df$P_value, method = "bonferroni")

View(df)

significant_features_BH <- df[df$P_adj_BH < 0.05, ]
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
res <- read_excel("C:/Users/12868/OneDrive/Desktop/Statistical_Results.xlsx")

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

