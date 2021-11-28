

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
setwd("~/Downloads/MR_analysis_cwp_paper/major_depression/")
#mdd <- fread("daner_pgc_mdd_meta_w2_no23andMe_rmUKBB_ PMID29700475.gz")
#mdd <- mdd[,c(1:7,9:11)]
#mdd$FREQ <- (mdd$FRQ_A_45396+mdd$FRQ_U_97250)/2
#mdd$N <- 142646
#mdd$BETA <- log(mdd$OR) # se and p-value not converted.
#mdd <- mdd[,c(1:5,11,13,9,10,12)]
#fwrite(mdd, file = "mdd_noukb_23me_readyforMR.txt", quote = FALSE, sep = "\t")

mdd_out_dat <- read_outcome_data(snps = NULL,
                                 filename = "mdd_noukb_23me_readyforMR.txt",
                                 sep = "\t",
                                 snp_col = "SNP",
                                 beta_col = "BETA",
                                 se_col = "SE", 
                                 eaf_col = "FREQ",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2", 
                                 pval_col = "P",
                                 samplesize_col = "N")

mdd_out_dat[,"outcome"]=mdd_out_dat[,"id.outcome"]="MDD"


# Harmonising again: 
dat11 <- harmonise_data(exposure_dat = cwp_ex_dat, outcome_dat = mdd_out_dat) # 39 SNPs remained
dat11 <- dat11[(dat11$'eaf.exposure' < 0.95) & (dat11$'eaf.exposure' > 0.05) & (dat11$'eaf.outcome' < 0.95) & (dat11$'eaf.outcome' > 0.05), ]
dat11=dat11[dat11$mr_keep==TRUE,] # 35 instruments identified.


mr_report(dat11, output_type = "html",output_path = "~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_MDD/")

cwp_to_mdd <- mr(dat11, method_list =c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                       "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))

fwrite(cwp_to_mdd, file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_MDD/cwp_to_mdd.txt")

mr_presso <- run_mr_presso(dat11, NbDistribution = 3000, SignifThreshold = 0.05) 

mr_presso 

#MR.RAPS (Robust Adjusted Profile Score) is a recently proposed method that considers the 
#measurement error in SNP-exposure effects, is unbiased when there are many (e.g. hundreds of) weak instruments, 
#and is robust to systematic and idiosyncratic pleiotropy.
mr_raps <- mr_raps(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#b 4.970905
#se 0.9233407
#pval 7.300694e-08
#nsnp 35

# mr ivw radial: for quantifying heterogeneity and outlier detection
# Raidal MR is an simple alternative to MR Presso
mr_ivw_radial <- mr_ivw_radial(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#F-statistic: 31.06 on 1 and 34 DF, p-value: 3.09e-06
#Q-Statistic for heterogeneity: 52.55606 on 34 DF , p-value: 0.02201132










## High Resulation PLotting 

# leave one out plot
leave_one_out <- mr_leaveoneout(dat11)
p3 <- mr_leaveoneout_plot(leave_one_out)
ggsave(p3[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_MDD/leave_one_out_cwp_to_mdd.pdf", width=7, height=7, dpi = 300)

# funnel plot
funnel_plot <- mr_singlesnp(dat11,method_list=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                                "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p4 <- mr_funnel_plot(funnel_plot)
ggsave(p4[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_MDD/funnel_plot_cwp_to_mdd.pdf", width=7, height=7, dpi = 300)


# scatter plot
p1 <- mr_scatter_plot(cwp_to_mdd, dat11)
length(p1)
ggsave(p1[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_MDD/scatter_plot_cwp_to_mdd.pdf", width=7, height=7, dpi = 300)


# forrest plot 
forr_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                              "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p2 <- mr_forest_plot(forr_plot)
p2[[1]]
ggsave(p2[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/CWP_to_BMI_MR/CWP_to_MDD/forest_plot_cwp_to_mdd.pdf", width=7, height=7, dpi = 300)


a <- directionality_test(dat11)
#id.exposure id.outcome exposure outcome snp_r2.exposure snp_r2.outcome correct_causal_direction
#        CWP        MDD      CWP     MDD     0.003192998   0.0007507444                     TRUE
#steiger_pval
# 1.536851e-18
