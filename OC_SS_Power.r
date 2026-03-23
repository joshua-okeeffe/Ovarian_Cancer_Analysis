###############################################
#Sample Size, Power & OC Analysis
###############################################

# --- Setup and data import ---

getwd()  #Check current working directory

BS3 <- read.csv("OC_Dataset.csv")  #Read in dataset

#Load required packages
library(FSA)          #For effect size utilities / fisheries stats
library(rcompanion)   #For ES.w2 and other effect size functions
library(sjstats)      #For anova_stats
library(effsize)      #For Cohen's d
library(pwr)          #For power and sample size calculations
library(tidyr)        #For data reshaping
library(dplyr)        #For data manipulation
library(dbplyr)       #For database manipulation 

attach(BS3)           #Attach for direct column reference 

colnames(BS3)         #Inspect variable names


###############################################
#Helper vectors for protein labels
###############################################

#For Rucaparib arm (n = 40 per marker)
PD_ruclist   <- rep("PD",   times = 40)
PDL1_ruclist <- rep("PDL1", times = 40)
CTLA_ruclist <- rep("CTLA4", times = 40)
SIRP_ruclist <- rep("SIRP", times = 40)

#For Niraparib arm (n = 48 per marker)
PD_nirlist   <- rep("PD",   times = 48)
PDL1_nirlist <- rep("PDL1", times = 48)
CTLA_nirlist <- rep("CTLA4", times = 48)
SIRP_nirlist <- rep("SIRP", times = 48)


###############################################
# Q1: Between-arm differences at baseline
# Are baseline checkpoint levels different between
# Rucaparib vs Niraparib treatment arms?
###############################################

#Subset baseline values by treatment arm

PD1_ruc   <- PD.1_pre[ Treatment_Arm == "Rucaparib" ]
PD1_nir   <- PD.1_pre[ Treatment_Arm == "niraparib" ]

PDL1_ruc  <- PD.L1_pre[ Treatment_Arm == "Rucaparib" ]
PDL1_nir  <- PD.L1_pre[ Treatment_Arm == "niraparib" ]

CTLA4_ruc <- CTLA.4_pre[ Treatment_Arm == "Rucaparib" ]
CTLA4_nir <- CTLA.4_pre[ Treatment_Arm == "niraparib" ]

SIRP_ruc  <- SIRP_pre[ Treatment_Arm == "Rucaparib" ]
SIRP_nir  <- SIRP_pre[ Treatment_Arm == "niraparib" ]


#Normality checks (Shapiro–Wilk) for each arm/marker

#PD-1
qqnorm(PD.1_pre); qqline(PD.1_pre, col = "red")
pd1_ruc_shap <- shapiro.test(PD1_ruc)
pd1_nir_shap <- shapiro.test(PD1_nir)

#PD-L1
qqnorm(PD.L1_pre); qqline(PD.L1_pre, col = "red")
pdl1_ruc_shap <- shapiro.test(PDL1_ruc)
pdl1_nir_shap <- shapiro.test(PDL1_nir)

#CTLA-4
qqnorm(CTLA.4_pre); qqline(CTLA.4_pre, col = "red")
ctla4_ruc_shap <- shapiro.test(CTLA4_ruc)
ctla4_nir_shap <- shapiro.test(CTLA4_nir)

#SIRP
qqnorm(SIRP_pre); qqline(SIRP_pre, col = "red")
sirp_ruc_shap <- shapiro.test(SIRP_ruc)
sirp_nir_shap <- shapiro.test(SIRP_nir)

#Collect Shapiro P-values per arm
ruc_pval <- c(pd1_ruc_shap$p.value,
              pdl1_ruc_shap$p.value,
              ctla4_ruc_shap$p.value,
              sirp_ruc_shap$p.value)

nir_pval <- c(pd1_nir_shap$p.value,
              pdl1_nir_shap$p.value,
              ctla4_nir_shap$p.value,
              sirp_nir_shap$p.value)

#Bonferroni adjustment for multiple normality tests
adjust_ruc <- p.adjust(ruc_pval, method = "bonferroni")
print(adjust_ruc)

adjust_nir <- p.adjust(nir_pval, method = "bonferroni")
print(adjust_nir)


#Parametric / non-parametric testing between arms at baseline
#PD-1 (baseline): t-test with equal variances

var(PD1_ruc)
var(PD1_nir)
var.test(PD1_ruc, PD1_nir)               #Homogeneity of variances
t.test(PD1_ruc, PD1_nir, var.equal = TRUE)  #Between-arm comparison

#PD-L1 (baseline): t-test, effect size, and power

var(PDL1_ruc)
var(PDL1_nir)
var.test(PDL1_ruc, PDL1_nir)

pdl1_t <- t.test(PDL1_ruc, PDL1_nir, var.equal = TRUE)
pdl1_t$p.value

#Cohen's d for PD-L1 and corresponding power / sample size
pdl1_cohen <- cohen.d(PDL1_ruc, PDL1_nir)
pdl1_cohen

pdl1_power <- pwr.t.test(
  d      = pdl1_cohen$estimate,
  sig    = 0.0125,     #Bonferroni-adjusted alpha for 4 markers
  power  = 0.8,
  type   = "two.sample"
)
pdl1_power

#CTLA-4 (baseline): both Wilcoxon and t-test for robustness
ctla4_w <- wilcox.test(CTLA4_ruc, CTLA4_nir)
ctla4_w

ctla4_t <- t.test(CTLA4_ruc, CTLA4_nir)
ctla4_t

ctla4_cohen <- cohen.d(CTLA4_ruc, CTLA4_nir)
ctla4_cohen

ctla4_power <- pwr.t.test(
  d      = ctla4_cohen$estimate,
  sig    = 0.0125,
  power  = 0.8,
  type   = "two.sample"
)
ctla4_power

#SIRP (baseline): Wilcoxon, t-test, Cohen's d, and power
sirp_w <- wilcox.test(SIRP_ruc, SIRP_nir)
sirp_w

sirp_t <- t.test(SIRP_ruc, SIRP_nir)
sirp_t

sirp_cohen <- cohen.d(SIRP_ruc, SIRP_nir)
sirp_cohen

sirp_power <- pwr.t.test(
  d      = sirp_cohen$estimate,
  sig    = 0.0125,
  power  = 0.8,
  type   = "two.sample"
)
sirp_power

#Vector of required sample sizes per group for non-significant markers
required      <- c(pdl1_power$n, ctla4_power$n, sirp_power$n)
required

#Optionally inflate CTLA4/SIRP requirements by 20% for anticipated drop-out
adj_required  <- c(pdl1_power$n,
                   ctla4_power$n * 1.2,
                   sirp_power$n * 1.2)
adj_required


###############################################
# Q2: Differential marker change by PARP inhibitor
# Are the markers affected differently by Rucaparib vs Niraparib?
###############################################

#Split dataset by treatment arm and prepare change scores

ruc_data <- BS3 %>%
  filter(Treatment_Arm == "Rucaparib")

nir_data <- BS3 %>%
  filter(Treatment_Arm == "niraparib")

#Remove variables not needed for marker-change analysis
ruc_data <- ruc_data %>%
  select(
    -Treatment_Response, -ER_status,
    -BMI_pre_treatment, -BMI_post_treatment,
    -PD.1_during, -PD.L1_during, -CTLA.4_during, -SIRP_during
  )

nir_data <- nir_data %>%
  select(
    -Treatment_Response, -ER_status,
    -BMI_pre_treatment, -BMI_post_treatment,
    -PD.1_during, -PD.L1_during, -CTLA.4_during, -SIRP_during
  )

#Create change-score variables (post - pre) per marker
ruc_data <- ruc_data %>%
  mutate(
    PD1_Diff   = PD.1_post - PD.1_pre,
    PDL1_diff  = PD.L1_post - PD.L1_pre,
    CTLA4_Diff = CTLA.4_post - CTLA.4_pre,
    SIRP_Diff  = SIRP_post - SIRP_pre
  )

nir_data <- nir_data %>%
  mutate(
    PD1_Diff   = PD.1_post - PD.1_pre,
    PDL1_diff  = PD.L1_post - PD.L1_pre,
    CTLA4_Diff = CTLA.4_post - CTLA.4_pre,
    SIRP_Diff  = SIRP_post - SIRP_pre
  )

#variants with only the difference variables retained
ruc_adapted <- ruc_data %>%
  select(-PD.1_pre, -PD.1_post,
         -PD.L1_pre, -PD.L1_post,
         -CTLA.4_pre, -CTLA.4_post,
         -SIRP_pre, -SIRP_post)

nir_adapted <- nir_data %>%
  select(-PD.1_pre, -PD.1_post,
         -PD.L1_pre, -PD.L1_post,
         -CTLA.4_pre, -CTLA.4_post,
         -SIRP_pre, -SIRP_post)


#Rucaparib: one-way ANOVA on change scores across 4 markers

#Combine all change scores into a single vector and label by marker
All_ruc_Val <- c(
  ruc_data$PD1_Diff,
  ruc_data$PDL1_diff,
  ruc_data$CTLA4_Diff,
  ruc_data$SIRP_Diff
)

shapiro.test(All_ruc_Val)  #Normality check of pooled differences

All_ruc_pro <- c(PD_ruclist,
                 PDL1_ruclist,
                 CTLA_ruclist,
                 SIRP_ruclist)

ruc_full <- data.frame(
  All_ruc_pro = factor(All_ruc_pro),
  All_ruc_Val = All_ruc_Val
)

ruc_aov <- aov(All_ruc_Val ~ All_ruc_pro, data = ruc_full)
summary(ruc_aov)

#Post-hoc comparisons between markers within Rucaparib arm
TukeyHSD(ruc_aov)

#Effect size and power calculation for Rucaparib
effectsize::cohens_f(ruc_aov)  #Overall Cohen's f

ruc_a <- anova_stats(ruc_aov)  #From sjstats: ANOVA summary including Cohen's f
ruc_a$cohens.f[1]

#Example power calculation for a relatively large effect (f ≈ 0.8)
pwr.anova.test(f = 0.8, k = 4, power = 0.8, sig.level = 0.0125)


#Niraparib: one-way ANOVA on change scores across 4 markers

All_nir_val <- c(
  nir_data$PD1_Diff,
  nir_data$PDL1_diff,
  nir_data$CTLA4_Diff,
  nir_data$SIRP_Diff
)

shapiro.test(All_nir_val)

All_nir_pro <- c(PD_nirlist,
                 PDL1_nirlist,
                 CTLA_nirlist,
                 SIRP_nirlist)

nir_full <- data.frame(
  All_nir_pro = factor(All_nir_pro),
  All_nir_val = All_nir_val
)

nir_aov <- aov(All_nir_val ~ All_nir_pro, data = nir_full)
summary(nir_aov)

TukeyHSD(nir_aov)

effectsize::cohens_f(nir_aov)
nir_a <- anova_stats(nir_aov)
nir_a$cohens.f[1]

#Example: power calculation for observed effect size (≈ 0.76)
pwr.anova.test(f = 0.76, k = 4, power = 0.8, sig.level = 0.0125)


#Post-treatment expression by marker (within arm)

#Rucaparib: post-treatment levels across markers
All_ruc_post_Val <- c(
  ruc_data$PD.1_post,
  ruc_data$PD.L1_post,
  ruc_data$CTLA.4_post,
  ruc_data$SIRP_post
)

All_ruc_post_pro <- c(PD_ruclist,
                      PDL1_ruclist,
                      CTLA_ruclist,
                      SIRP_ruclist)

ruc_post_full <- data.frame(
  All_ruc_post_pro = factor(All_ruc_post_pro),
  All_ruc_post_Val = All_ruc_post_Val
)

ruc_post_aov <- aov(All_ruc_post_Val ~ All_ruc_post_pro, data = ruc_post_full)
summary(ruc_post_aov)
TukeyHSD(ruc_post_aov)

effectsize::cohens_f(ruc_post_aov)
ruc_post_a <- anova_stats(ruc_post_aov)
ruc_post_a$cohens.f[1]

pwr.anova.test(f = 0.52, k = 4, power = 0.8, sig.level = 0.0125)


#Niraparib: post-treatment levels across markers
All_nir_post_val <- c(
  nir_data$PD.1_post,
  nir_data$PD.L1_post,
  nir_data$CTLA.4_post,
  nir_data$SIRP_post
)

All_nir_post_pro <- c(PD_nirlist,
                      PDL1_nirlist,
                      CTLA_nirlist,
                      SIRP_nirlist)

nir_post_full <- data.frame(
  All_nir_post_pro = factor(All_nir_post_pro),
  All_nir_post_val = All_nir_post_val
)

nir_post_aov <- aov(All_nir_post_val ~ All_nir_post_pro, data = nir_post_full)
summary(nir_post_aov)
TukeyHSD(nir_post_aov)

effectsize::cohens_f(nir_post_aov)
nir_post_a <- anova_stats(nir_post_aov)
nir_post_a$cohens.f[1]

pwr.anova.test(f = 0.50, k = 4, power = 0.8, sig.level = 0.0125)


###############################################
# Q3: Underweight patients and required sample size
# Is there a significant increase in underweight (BMI < 20)
# within/between arms, and what n is needed for 80% power?
###############################################

#Change from >20 to <20 BMI within each arm

#Number moving from BMI > 20 to < 20 (Rucaparib)
ruc_length_change <- length(
  BMI_pre_treatment[
    Treatment_Arm == "Rucaparib" &
      BMI_pre_treatment > 20 &
      BMI_post_treatment < 20
  ]
)
ruc_length_change

#Complement: Rucaparib total - changed
(length(BMI_post_treatment[Treatment_Arm == "Rucaparib"])) - ruc_length_change

#Niraparib
nir_length_change <- length(
  BMI_pre_treatment[
    Treatment_Arm == "niraparib" &
      BMI_pre_treatment > 20 &
      BMI_post_treatment < 20
  ]
)
nir_length_change

(length(BMI_post_treatment[Treatment_Arm == "niraparib"])) - nir_length_change

#2x2 table: treatment vs change/no change (example counts)
New_matrix <- matrix(
  c(10, 6, 30, 42),
  nr = 2,
  dimnames = list(
    "Treat"  = c("Rucaparib", "Niraparib"),
    "Change" = c("Change", "No Change")
  )
)
New_matrix

mcnemar.test(New_matrix)

#Effect size (w2) and power for this table
w2 <- ES.w2(New_matrix / sum(New_matrix))
pwr.chisq.test(w=w2, df=3, sig.level=0.05, power=0.8)


#Pre vs post underweight counts by arm

length(BMI_pre_treatment[BMI_pre_treatment < 20 & Treatment_Arm == "Rucaparib"])
length(BMI_pre_treatment[BMI_pre_treatment < 20 & Treatment_Arm == "niraparib"])

length(BMI_pre_treatment[BMI_post_treatment < 20 & Treatment_Arm == "Rucaparib"])
length(BMI_pre_treatment[BMI_post_treatment < 20 & Treatment_Arm == "niraparib"])

BMI_matrix <- matrix(
  c(4, 3, 12, 9),
  nr = 2,
  dimnames = list(
    "treat" = c("Ruc", "Nir"),
    "Time"  = c("Pre", "Post")
  )
)
BMI_matrix

mcnemar.test(BMI_matrix)

w3 <- ES.w2(BMI_matrix / sum(BMI_matrix))
w3
#w3 too close to zero to calculate a meaningful power size


#Overall pre vs post BMI (collapsed over treatment)

matrixA <- matrix(
  c(36, 28, 4, 12),
  nr = 2,
  dimnames = list(
    "Time" = c("Pre", "Post"),
    "BMI"  = c(">20", "<20")
  )
)
matrixA
mcnemar.test(matrixA)


#Separate McNemar and power for each arm

#Rucaparib: 2x2 pre vs post BMI table
length(BMI_pre_treatment[Treatment_Arm == "Rucaparib" &
                           BMI_pre_treatment > 20 &
                           BMI_post_treatment > 20])
length(BMI_pre_treatment[Treatment_Arm == "Rucaparib" &
                           BMI_pre_treatment > 20 &
                           BMI_post_treatment < 20])
length(BMI_pre_treatment[Treatment_Arm == "Rucaparib" &
                           BMI_pre_treatment < 20 &
                           BMI_post_treatment > 20])
length(BMI_pre_treatment[Treatment_Arm == "Rucaparib" &
                           BMI_pre_treatment < 20 &
                           BMI_post_treatment < 20])

matrix_ruc <- matrix(
  c(26, 10, 2, 2),
  nr = 2,
  dimnames = list(
    "Post" = c(">20", "<20"),
    "Pre"  = c(">20", "<20")
  )
)
matrix_ruc

ruc_mcn <- mcnemar.test(matrix_ruc)
wr      <- ES.w2(matrix_ruc / sum(matrix_ruc))
pwr.chisq.test(w = wr, df = 3, sig.level = 0.05, power = 0.8)


#Niraparib: analogous 2x2 table
length(BMI_pre_treatment[Treatment_Arm == "niraparib" &
                           BMI_pre_treatment > 20 &
                           BMI_post_treatment > 20])
length(BMI_pre_treatment[Treatment_Arm == "niraparib" &
                           BMI_pre_treatment > 20 &
                           BMI_post_treatment < 20])
length(BMI_pre_treatment[Treatment_Arm == "niraparib" &
                           BMI_pre_treatment < 20 &
                           BMI_post_treatment > 20])
length(BMI_pre_treatment[Treatment_Arm == "niraparib" &
                           BMI_pre_treatment < 20 &
                           BMI_post_treatment < 20])

matrix_nir <- matrix(
  c(39, 6, 0, 3),
  nr = 2,
  dimnames = list(
    "Post" = c(">20", "<20"),
    "Pre"  = c(">20", "<20")
  )
)
matrix_nir

nir_mcn <- mcnemar.test(matrix_nir)
wn      <- ES.w2(matrix_nir / sum(matrix_nir))
pwr.chisq.test(w = wn, df = 3, sig.level = 0.05, power = 0.8)

#Bonferroni-adjusted P for within-arm McNemar tests
p_adj       <- c(ruc_mcn$p.value, nir_mcn$p.value)
adjusted_p  <- p.adjust(p_adj, method = "bonferroni")
adjusted_p


###############################################
# Q4: Two-proportion power calculation
###############################################

#Compute Cohen's h for two proportions (e.g. 0.31 vs 0.23)
h <- ES.h(0.31, 0.23)

#Unequal sample-size two-proportion test:
#here n2 fixed at 4,135,000 and power/alpha specified
pwr.2p2n.test(
  h        = h,
  n2       = 4135000,
  sig.level = 0.05,
  power     = 0.8
)
