# ABOUT: Visualize and test distribution of MAF scores in both case and control MDD (Major Depressive Disorder) populations (Wray et al. 2018)
# AUTHOR: Koen Rademaker
# DATE: 25 March 2019


# Install and/or load ggplot2
install.packages("ggplot2")
library(ggplot2)
# Load data
maf_comparison <- read.delim("~/Data/MDD_MAF_comparison/maf_comparison.txt", row.names=NULL)
# Set up DataFrame structure
case_df <- data.frame(group = "MDD cases", value = maf_comparison$MAF_cases)
control_df <- data.frame(group = "MDD controls", value = maf_comparison$MAF_controls)
plot.data <- rbind(case_df, control_df)
# Generate boxplot and violin plot
ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot()
ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_violin()
# Run chi-squared test for independence
t <- table(case_df$value, control_df$value)
chisq.test(t)
