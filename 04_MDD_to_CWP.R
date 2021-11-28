# Direction of effect MDD to CWP


# load MDD data (exposure)
setwd("~/Downloads/MR_analysis_cwp_paper/major_depression/")

#mdd <- fread('mdd_exposure.txt') #67 instruments

eu_mdd_dat <- read_exposure_data("mdd_exposure.txt",
                                 clump = FALSE,
                                 sep = "\t",
                                 snp_col = "SNP",
                                 beta_col = "BETA",
                                 se_col = "SE", 
                                 eaf_col = "FREQ",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2", 
                                 pval_col = "P",
                                 samplesize_col = "N")

eu_mdd_dat [,"exposure"]=eu_mdd_dat[,"id.exposure"]="MDD" # 77 SNPs



# Load CWP data (outcome)
setwd("~/Downloads/MR_analysis_cwp_paper/CWP_GWAS")
cwp <- fread("all_samples_final.csv")
cwp <- cwp[,1:10]
colnames(cwp)[6] <- "FREQ"
cwp$N <- 249843
fwrite(cwp, file = "all_samples_final.csv",  quote = FALSE, row.names = FALSE)

cwp_out_dat <- read_outcome_data(snps = NULL,
                                 filename = "all_samples_final.csv",
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "BETA",
                                 se_col = "SE",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2",
                                 eaf_col = "FREQ",
                                 pval_col = "P",
                                 samplesize_col = "N")

cwp_out_dat[,"outcome"]=cwp_out_dat[,"id.outcome"]="CWP"


# Harmonization: 
dat11 <- harmonise_data(exposure_dat = eu_mdd_dat, outcome_dat = cwp_out_dat) 
dat11 <- dat11[(dat11$'eaf.exposure' < 0.95) & (dat11$'eaf.exposure' > 0.05) & (dat11$'eaf.outcome' < 0.95) & (dat11$'eaf.outcome' > 0.05), ]
dat11=dat11[dat11$mr_keep==TRUE,] # 51 SNPs left as instrument.


### Run MR-analysis ###

het<-mr_heterogeneity(dat11,method_list=c("mr_egger_regression", "mr_ivw"))
#id.exposure id.outcome outcome exposure                    method        Q Q_df    Q_pval
#1         MDD     2mCFMo outcome      MDD                  MR Egger 58.96238   49 0.1558680
#2         MDD     2mCFMo outcome      MDD Inverse variance weighted 62.18220   50 0.1156912
plt<-mr_pleiotropy_test(dat11)
#id.exposure id.outcome outcome exposure egger_intercept           se      pval
#1         MDD     2mCFMo outcome      MDD    0.0008129336 0.0004969684 0.1082928

mr_report(dat11,output_type = "html",
          output_path = "~/Downloads/MR_analysis_cwp_paper/major_depression/MDD_to_CWP/")

mdd_to_cwp <- mr(dat11, method_list =c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                       "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))

fwrite(mdd_to_cwp, file="~/Downloads/MR_analysis_cwp_paper/major_depression/MDD_to_CWP/mdd_to_cwp.txt")


mr_presso <- run_mr_presso(dat11, NbDistribution = 3000, SignifThreshold = 0.05) 
#[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
#[1] 0.167


#MR.RAPS (Robust Adjusted Profile Score) is a recently proposed method that considers the 
#measurement error in SNP-exposure effects, is unbiased when there are many (e.g. hundreds of) weak instruments, 
#and is robust to systematic and idiosyncratic pleiotropy.
mr_raps <- mr_raps(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#$b 0.01124917
#$se 0.001775911
#$pval 2.384046e-10
#$nsnp 51

# mr ivw radial: for quantifying heterogeneity and outlier detection
# Raidal MR is an simple alternative to MR Presso
mr_ivw_radial <- mr_ivw_radial(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())

#F-statistic: 42.75 on 1 and 50 DF, p-value: 3.17e-08
#Q-Statistic for heterogeneity: 59.51213 on 50 DF , p-value: 0.1678876

## High Resulation PLotting 

# leave one out plot
leave_one_out <- mr_leaveoneout(dat11)
p3 <- mr_leaveoneout_plot(leave_one_out)
ggsave(p3[[1]], file="~/Downloads/MR_analysis_cwp_paper/major_depression/MDD_to_CWP/leave_one_out_mdd.pdf", width=7, height=7, dpi = 300)

# funnel plot
funnel_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                                 "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p4 <- mr_funnel_plot(funnel_plot)
ggsave(p4[[1]], file="~/Downloads/MR_analysis_cwp_paper/major_depression/MDD_to_CWP/funnel_plot_mdd.pdf", width=7, height=7, dpi = 300)


# scatter plot
p1 <- mr_scatter_plot(mdd_to_cwp, dat11)
length(p1)
ggsave(p1[[1]], file="~/Downloads/MR_analysis_cwp_paper/major_depression/MDD_to_CWP/scatter_plot_mdd.pdf", width=7, height=7, dpi = 300)


# forrest plot 
forr_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                              "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p2 <- mr_forest_plot(forr_plot)
p2[[1]]
ggsave(p2[[1]], file="~/Downloads/MR_analysis_cwp_paper/major_depression/MDD_to_CWP/forest_plot_mdd.pdf", width=7, height=7, dpi = 300)

a <- directionality_test(dat11)
#id.exposure id.outcome exposure outcome snp_r2.exposure snp_r2.outcome correct_causal_direction
#1         MDD     2mCFMo      MDD outcome     0.008236166   0.0004623377                     TRUE
#steiger_pval
#1 2.214457e-97














