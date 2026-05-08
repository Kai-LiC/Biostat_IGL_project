# IGL Electrophysiological Statistical Analysis

This repository contains R scripts for electrophysiological data analysis.

## Statistical Analysis

### Two-sample comparison

1. Normality was assessed using the Shapiro–Wilk test.
   - p > 0.05: Welch’s t-test was used.
   - p ≤ 0.05: Mann–Whitney U test was used.

2. Multiple-comparisons correction was applied when appropriate.

3. Data visualization and statistical plotting were performed in

### LMM/GLMM
1. Construct Linear Mixed Model (LMM)

2. Test normality
- Skewness < 1
- Shapiro–Wilk normality test (p > 0.05)

If normally distributed:
- Use LMM

Optional step 1: Data transformation
- Log transformation
- Square root transformation

If normally distributed after transformation:
- Use LMM on transformed data

If still not normally distributed:
- Use Gamma Generalized Linear Mixed Model (GLMM)

Optional step 2: Gamma GLMM
- Applied for right-skewed positive data
- Model fitness was evaluated before interpretation

3. Multiple-comparisons correction
   
4. Plotting
- Statistical visualization performed in R
