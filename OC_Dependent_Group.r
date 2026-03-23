#########################################
# OC Statistical Analysis 2

# Purpose: Analyze immune checkpoint expression and BMI changes 
#          pre-, during-, and post-treatment
#########################################

# --- Setup and Data Import ---

#Check current working directory
getwd()

#Import dataset
Ass2 <- read.csv("OC_Dataset.csv")

#Inspect dataset structure
colnames(Ass2)
summary(Ass2)
class(Ass2)

#Attach columns for direct referencing
attach(Ass2)


#########################################
# Q1: Difference in immune checkpoint expression pre vs post treatment
#########################################

# Boxplots for each marker pre- and post-treatment
boxplot(PD.1_pre, PD.1_post, names = c("PD1_pre", "PD1_post"),
        ylab = "Expression Level", xlab = "Time", col = c("red", "blue"))

boxplot(PD.L1_pre, PD.L1_post, names = c("PDL1_pre", "PDL1_post"),
        ylab = "Expression Level", xlab = "Time", col = c("red", "blue"))

boxplot(CTLA.4_pre, CTLA.4_post, names = c("CTLA4_pre", "CTLA4_post"),
        ylab = "Expression Level", xlab = "Time", col = c("red", "blue"))

boxplot(SIRP_pre, SIRP_post, names = c("SIRP_pre", "SIRP_post"),
        ylab = "Expression Level", xlab = "Time", col = c("red", "blue"))

#Test for normality of paired differences (Shapiro-Wilk)
PD1_diff   <- PD.1_post   - PD.1_pre
PDL1_diff  <- PD.L1_post  - PD.L1_pre
CTLA4_diff <- CTLA.4_post - CTLA.4_pre
SIRP_diff  <- SIRP_post   - SIRP_pre

shapiro.test(PD1_diff)
shapiro.test(PDL1_diff)
shapiro.test(CTLA4_diff)
shapiro.test(SIRP_diff)

#Non-parametric test (since likely not normal): Wilcoxon signed-rank test
wilcox.test(PD.1_post,  PD.1_pre,  paired = TRUE)
wilcox.test(PD.L1_post, PD.L1_pre, paired = TRUE)
wilcox.test(CTLA.4_post, CTLA.4_pre, paired = TRUE)
wilcox.test(SIRP_post,   SIRP_pre,   paired = TRUE)


#########################################
# Q2: Change in underweight patients (BMI < 20) post-treatment
#########################################

#Determine counts for BMI classification changes
length(BMI_pre_treatment[BMI_pre_treatment > 20 & BMI_post_treatment > 20])  # Pre and post normal
length(BMI_pre_treatment[BMI_pre_treatment > 20 & BMI_post_treatment < 20])  # Became underweight
length(BMI_pre_treatment[BMI_pre_treatment < 20 & BMI_post_treatment > 20])  # Recovered to normal
length(BMI_pre_treatment[BMI_pre_treatment < 20 & BMI_post_treatment < 20])  # Stayed underweight

#Create 2x2 contingency table (McNemar input)
BMI_Matrix <- matrix(c(65, 16, 2, 5), nrow = 2,
                     dimnames = list("BMI_Post" = c("Post Over", "Post Under"),
                                     "BMI_Pre"  = c("Pre Over", "Pre Under")))
print(BMI_Matrix)

#McNemar’s test for paired nominal data (used for pre/post comparison)
mcnemar.test(BMI_Matrix)

#Binomial test on direction of change (16 decreased, 2 increased)
binom.test(x = 16, n = 18, p = 0.5, alternative = "less")

#Interpretation: ~89% chance results not due to random variation.


#########################################
# Q3: Change in protein expression over 3 time points (pre, during, post)
#########################################

#Create data frames for each protein across time points
PD1_Table   <- data.frame(PD.1_pre, PD.1_during, PD.1_post)
PDL1_Table  <- data.frame(PD.L1_pre, PD.L1_during, PD.L1_post)
CTLA4_Table <- data.frame(CTLA.4_pre, CTLA.4_during, CTLA.4_post)
SIRP_Table  <- data.frame(SIRP_pre, SIRP_during, SIRP_post)

#Friedman test (non-parametric equivalent of repeated-measures ANOVA)
friedman.test(as.matrix(PD1_Table[,  c("PD.1_pre",  "PD.1_during",  "PD.1_post")]))
friedman.test(as.matrix(PDL1_Table[, c("PD.L1_pre", "PD.L1_during", "PD.L1_post")]))
friedman.test(as.matrix(CTLA4_Table[, c("CTLA.4_pre", "CTLA.4_during", "CTLA.4_post")]))
friedman.test(as.matrix(SIRP_Table[,  c("SIRP_pre",  "SIRP_during",  "SIRP_post")]))

#Post-hoc Conover comparisons (requires PMCMRplus)
#Omit NA's to prevent errors from occurring
#install.packages("PMCMRplus")  if not installed
library(PMCMRplus)

frdAllPairsConoverTest(as.matrix(na.omit(PD1_Table)))
frdAllPairsConoverTest(as.matrix(na.omit(PDL1_Table)))
frdAllPairsConoverTest(as.matrix(na.omit(CTLA4_Table)))
frdAllPairsConoverTest(as.matrix(na.omit(SIRP_Table)))

#Boxplots for visualization across 3 time points
boxplot(PD.1_pre, PD.1_during, PD.1_post,
        names = c("PD1_pre", "PD1_during", "PD1_post"),
        ylab = "Expression Level", xlab = "Time",
        col = c("red", "green", "blue"))

boxplot(PD.L1_pre, PD.L1_during, PD.L1_post,
        names = c("PDL1_pre", "PDL1_during", "PDL1_post"),
        ylab = "Expression Level", xlab = "Time",
        col = c("red", "green", "blue"))

boxplot(CTLA.4_pre, CTLA.4_during, CTLA.4_post,
        names = c("CTLA4_pre", "CTLA4_during", "CTLA4_post"),
        ylab = "Expression Level", xlab = "Time",
        col = c("red", "green", "blue"))

boxplot(SIRP_pre, SIRP_during, SIRP_post,
        names = c("SIRP_pre", "SIRP_during", "SIRP_post"),
        ylab = "Expression Level", xlab = "Time",
        col = c("red", "green", "blue"))
