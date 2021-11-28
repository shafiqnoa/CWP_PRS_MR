

library(TwoSampleMR)
library(MRInstruments)
library(data.table)


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

cwp_ex_dat[,"exposure"]=cwp_ex_dat[,"id.exposure"]="CWP" # 39 instuments 


# Outcome data
setwd("~/Downloads/MR_analysis_cwp_paper/AgeFirstBirth/")
#afb <- fread("AFB_ctg_vl.txt") # beta and se calculated from Z and P value.

afb_out_dat <- read_outcome_data(snps = NULL,
                                 filename = "AFB_ctg_vl.txt",
                                 sep = "\t",
                                 snp_col = "SNP",
                                 beta_col = "BETA",
                                 se_col = "SE", 
                                 eaf_col = "FREQ",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2", 
                                 pval_col = "P",
                                 samplesize_col = "N")

afb_out_dat[,"outcome"]=afb_out_dat[,"id.outcome"]="afb"


# Harmonising again: 
dat11 <- harmonise_data(exposure_dat = cwp_ex_dat, outcome_dat = afb_out_dat) # 39 SNPs remained
dat11 <- dat11[(dat11$'eaf.exposure' < 0.95) & (dat11$'eaf.exposure' > 0.05) & (dat11$'eaf.outcome' < 0.95) & (dat11$'eaf.outcome' > 0.05), ]
dat11=dat11[dat11$mr_keep==TRUE,] # 28 instruments identified.


mr_report(dat11, output_type = "html",output_path = "~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_AFB/")

cwp_to_afb <- mr(dat11, method_list =c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                       "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))

fwrite(cwp_to_afb, file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_AFB/cwp_to_afb.txt")

mr_presso <- run_mr_presso(dat11, NbDistribution = 3000, SignifThreshold = 0.05) 

mr_presso 
#[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
#[1] 0.072

#MR.RAPS (Robust Adjusted Profile Score) is a recently proposed method that considers the 
#measurement error in SNP-exposure effects, is unbiased when there are many (e.g. hundreds of) weak instruments, 
#and is robust to systematic and idiosyncratic pleiotropy.
mr_raps <- mr_raps(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#b -1.512536
#se  0.3081601
#pval 9.18779e-07
#nsnp 28

# mr ivw radial: for quantifying heterogeneity and outlier detection
# Raidal MR is an simple alternative to MR Presso
mr_ivw_radial <- mr_ivw_radial(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#F-statistic: 29.13 on 1 and 27 DF, p-value: 1.05e-05
#Q-Statistic for heterogeneity: 38.012 on 27 DF , p-value: 0.07767237


## High Resulation PLotting 

# leave one out plot
leave_one_out <- mr_leaveoneout(dat11,method_list=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                                    "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p3 <- mr_leaveoneout_plot(leave_one_out)
ggsave(p3[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_AFB/leave_one_out_cwp_to_afb.pdf", width=7, height=7, dpi = 300)

# funnel plot
funnel_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                                "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p4 <- mr_funnel_plot(funnel_plot)
ggsave(p4[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_AFB/funnel_plot_cwp_to_afb.pdf", width=7, height=7, dpi = 300)


# scatter plot
p1 <- mr_scatter_plot(cwp_to_afb, dat11)
length(p1)
ggsave(p1[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_AFB/scatter_plot_cwp_to_afb.pdf", width=7, height=7, dpi = 300)


# forrest plot 
forr_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                              "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p2 <- mr_forest_plot(forr_plot)
p2[[1]]
ggsave(p2[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_AFB/forest_plot_cwp_to_afb.pdf", width=7, height=7, dpi = 300)


a <- directionality_test(dat11)
#id.exposure id.outcome exposure outcome snp_r2.exposure snp_r2.outcome correct_causal_direction
#1         CWP        afb      CWP     afb     0.002544524    0.000336976                     TRUE
#steiger_pval
#1 5.909498e-30
