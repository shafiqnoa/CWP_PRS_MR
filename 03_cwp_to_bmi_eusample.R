

library(TwoSampleMR)
library(MRInstruments)
library(data.table)
library(ggplot2)


# Reverse MR analysis using 5e-08 snps as instruemnt for CWP

# Exposure data
#setwd("~/Downloads/MR_analysis_cwp_paper/CWP_GWAS")
#cwp <- fread("cwp_unification_output_done.csv")
#cwp_ex_ins <- cwp %>% filter(p<1e-05)
#fwrite(cwp_ex_ins, file="cwp_ex_ins.txt", sep = "\t")

setwd("~/Downloads/MR_analysis_cwp_paper/CWP_GWAS/")

cwp_exposure <- fread("cwp_exposure.txt")


cwp_ex_dat <- read_exposure_data("cwp_exposure.txt",
                                 clump = FALSE,
                                 sep = "\t",
                                  snp_col = "SNP",
                                 beta_col = "BETA",
                                 se_col = "SE",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2",
                                 eaf_col = "FREQ",
                                 pval_col = "P",
                                 samplesize_col = "N")

cwp_ex_dat[,"exposure"]=cwp_ex_dat[,"id.exposure"]="CWP" # 47 instuments 


# Outcome data
setwd("~/Downloads/MR_analysis_cwp_paper/bmi")
#bmi <- fread("eu_ancestry_bmi_2015.gz")
bmi_out_dat <- read_outcome_data(snps = NULL,
                                 filename = "eu_ancestry_bmi_2015.gz",
                                 sep = "\t",
                                 snp_col = "SNP",
                                 beta_col = "b",
                                 se_col = "se", 
                                 eaf_col = "Freq1.Hapmap",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2", 
                                 pval_col = "p",
                                 samplesize_col = "N")

bmi_out_dat[,"outcome"]=bmi_out_dat[,"id.outcome"]="BMI"


# Harmonising again: 
dat11 <- harmonise_data(exposure_dat = cwp_ex_dat, outcome_dat = bmi_out_dat) 
dat11 <- dat11[(dat11$'eaf.exposure' < 0.95) & (dat11$'eaf.exposure' > 0.05) & (dat11$'eaf.outcome' < 0.95) & (dat11$'eaf.outcome' > 0.05), ]
dat11=dat11[dat11$mr_keep==TRUE,] # 82 instruments identified.
dat11 <- dat11[complete.cases(dat11),] # 81 instruments identified.

#dat11$f_stat <- ((dat11$beta.exposure)^2)/((dat11$se.exposure)^2)
#a <- dat11[,31]
#mean(a) #18.60

summary(cwp_exposure$FREQ)
mr_report(dat11, output_type = "html",output_path = "~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/")

cwp_to_bmi <-  mr(dat11, method_list =c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                       "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))

fwrite(cwp_to_bmi, file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/cwp_to_bmi.csv")

het<-mr_heterogeneity(dat11,method_list=c("mr_egger_regression", "mr_ivw"))
plt<-mr_pleiotropy_test(dat11)

set.seed(10000)
mr_presso <- run_mr_presso(dat11, NbDistribution = 3000, SignifThreshold = 0.05) 

mr_presso 
[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
[1] 0.02266667

[[1]]$`MR-PRESSO results`$`Distortion Test`
[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
[1]  2 13
#MR.RAPS (Robust Adjusted Profile Score) is a recently proposed method that considers the 
#measurement error in SNP-exposure effects, is unbiased when there are many (e.g. hundreds of) weak instruments, 
#and is robust to systematic and idiosyncratic pleiotropy.
mr_raps <- mr_raps(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#$b
#[1] 1.33965
#$se
#[1] 0.2487647
#$pval
#[1] 7.236075e-08
#$nsnp
#[1] 81

# mr ivw radial: for quantifying heterogeneity and outlier detection
# Raidal MR is an simple alternative to MR Presso
mr_ivw_radial <- mr_ivw_radial(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#F-statistic: 25.21 on 1 and 80 DF, p-value: 3.05e-06
#Q-Statistic for heterogeneity: 106.7509 on 80 DF , p-value: 0.02453638


## High Resulation PLotting 

# leave one out plot
leave_one_out <- mr_leaveoneout(dat11)
p3 <- mr_leaveoneout_plot(leave_one_out)
ggsave(p3[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/leave_one_out_cwp_to_bmi.pdf", width=7, height=7, dpi = 300)

# funnel plot
funnel_plot <- mr_singlesnp(dat11,all_method=c("mr_egger_regression","mr_ivw_fe",
                                              "mr_weighted_median","mr_weighted_mode"))
p4 <- mr_funnel_plot(funnel_plot)
ggsave(p4[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/funnel_plot_cwp_to_bmi.pdf", width=7, height=7, dpi = 300)


# scatter plot
p1 <- mr_scatter_plot(cwp_to_bmi, dat11)
length(p1)
ggsave(p1[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/scatter_plot_cwp_to_bmi.pdf", width=7, height=7, dpi = 300)


# forrest plot 
forr_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_fe",
                                             "mr_weighted_median","mr_weighted_mode"))
p2 <- mr_forest_plot(forr_plot)
p2[[1]]
ggsave(p2[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/forest_plot_cwp_to_bmi.pdf", width=7, height=7, dpi = 300)


a <- directionality_test(dat11)

## plots created from this points are added in the manuscript
# Excluding three oulier snps (rs1015055, rs12480881, rs583514)
which(dat11$SNP=="rs583514")
which(dat11$SNP=="rs1015055")
which(dat11$SNP=="rs12480881")
which(dat11$SNP=="rs920986")
which(dat11$SNP=="rs9492472")

dat12 <- dat11[-c(2,13,48,67,76),]

cwp_to_bmi_sensitivity <-  mr(dat12, method_list =c("mr_egger_regression","mr_ivw_fe",
                                        "mr_weighted_median","mr_weighted_mode"))

fwrite(cwp_to_bmi_sensitivity, file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/cwp_to_bmi_sensitivity.csv")

het<-mr_heterogeneity(dat12,method_list=c("mr_egger_regression", "mr_ivw"))
plt<-mr_pleiotropy_test(dat12)

set.seed(10000)
mr_presso <- run_mr_presso(dat12, NbDistribution = 3000, SignifThreshold = 0.05) 
mr_presso 
mr_ivw_radial <- mr_ivw_radial(dat12$beta.exposure, dat12$beta.outcome, dat12$se.exposure, dat12$se.outcome, parameters = default_parameters())

#mr_raps <- mr_raps(dat12$beta.exposure, dat12$beta.outcome, dat12$se.exposure, dat12$se.outcome, parameters = default_parameters())

# High Resulation PLotting 

# leave one out plot
leave_one_out <- mr_leaveoneout(dat12)
p3 <- mr_leaveoneout_plot(leave_one_out)
ggsave(p3[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/leave_one_out_cwp_to_bmi_outlierexcluded.pdf", width=7, height=7, dpi = 300)

# funnel plot
funnel_plot <- mr_singlesnp(dat12,all_method=c("mr_egger_regression","mr_ivw_fe",
                                               "mr_weighted_median","mr_weighted_mode"))
p4 <- mr_funnel_plot(funnel_plot)
ggsave(p4[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/funnel_plot_cwp_to_bmi_outlierexcluded.pdf", width=7, height=7, dpi = 300)


# scatter plot
p1 <- mr_scatter_plot(cwp_to_bmi_sensitivity, dat12)
length(p1)
ggsave(p1[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/scatter_plot_cwp_to_bmi_outlierexcluded.pdf", width=7, height=7, dpi = 300)


# forrest plot 
forr_plot <- mr_singlesnp(dat12, all_method=c("mr_egger_regression","mr_ivw_fe",
                                              "mr_weighted_median","mr_weighted_mode"))
p2 <- mr_forest_plot(forr_plot)
p2[[1]]
ggsave(p2[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/forest_plot_cwp_to_bmi_outlierexcluded.pdf", width=7, height=7, dpi = 300)



exp(0.58)

#write.csv(dat11, file = "~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_BMI/cwp_to_bmi_81ins.csv")
getwd()

colnames(dat11)
dat11$r2_cwp <- 0
## varriance explained by each 81 instruments in CWP:
for (i in 1:nrow(dat11)){
  dat11$r2_cwp[i] <- (2*(dat11[i,6]^2)*(dat11[i,8])*(1-dat11[i,8]))/((2*(dat11[i,6]^2)*(dat11[i,8])) + (dat11[i,21]^2)*(2*dat11[i,23])*dat11[i,8]*(1-dat11[i,8]))
}
sum(dat11$r2_cwp)*100


## varriance explained by each 81 instruments in CWP:

# Formulla: F <- (n-k-1/k) * (R^2/(1-R^2))

F <- ((249843-81-1)/81)*((sum(dat11$r2_cwp))^2)/(1-(sum(dat11$r2_cwp))^2)

F*100
11.21198

#F <- ((233588-81-1)/81)*((sum(dat11$r2_bmi))^2)/(1-(sum(dat11$r2_bmi))^2)
#F*100
## varriance explained by each 81 instruments in BMI:
dat11$r2_bmi <- 0
## varriance explained by each 81 instruments in CWP:
for (i in 1:nrow(dat11)){
  dat11$r2_bmi[i] <- (2*(dat11[i,7]^2)*(dat11[i,9])*(1-dat11[i,9]))/((2*(dat11[i,7]^2)*(dat11[i,9])) + (dat11[i,14]^2)*(2*dat11[i,16])*dat11[i,9]*(1-dat11[i,9]))
}
sum(dat11$r2_bmi)*100

F <- ((249843-78-1)/78)*((sum(dat12$r2_cwp))^2)/(1-(sum(dat12$r2_cwp))^2)
F*100







dat12$r2_cwp <- 0
## varriance explained by each 81 instruments in CWP:
for (i in 1:nrow(dat12)){
  dat12$r2_cwp[i] <- (2*(dat12[i,6]^2)*(dat12[i,8])*(1-dat12[i,8]))/((2*(dat12[i,6]^2)*(dat12[i,8])) + (dat12[i,21]^2)*(2*dat12[i,23])*dat12[i,8]*(1-dat12[i,8]))
}
sum(dat12$r2_cwp)*100 #0.57%


## varriance explained by each 81 instruments in CWP:

# Formulla: F <- (n-k-1/k) * (R^2/(1-R^2))

F <- ((249843-76-1)/76)*((sum(dat12$r2_cwp))^2)/(1-(sum(dat12$r2_cwp))^2) #10.64

