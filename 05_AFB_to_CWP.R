# Direction of effect MDD to CWP


# load MDD data (exposure)
setwd("~/Downloads/MR_analysis_cwp_paper/AgeFirstBirth/")

#afb<- fread("afb_exposure.txt")
#afb$N <- 251151
#fwrite(afb, file="afb_exposure.txt", quote = FALSE, sep = "\t")

eu_afb_dat <- read_exposure_data("afb_exposure.txt",
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

eu_afb_dat [,"exposure"]=eu_afb_dat[,"id.exposure"]="Age at First Birth" # 77 SNPs



# Load CWP data (outcome)
setwd("~/Downloads/MR_analysis_cwp_paper/CWP_GWAS")

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
dat11 <- harmonise_data(exposure_dat = eu_afb_dat, outcome_dat = cwp_out_dat) 
dat11 <- dat11[(dat11$'eaf.exposure' < 0.95) & (dat11$'eaf.exposure' > 0.05) & (dat11$'eaf.outcome' < 0.95) & (dat11$'eaf.outcome' > 0.05), ]
dat11=dat11[dat11$mr_keep==TRUE,] # 41 SNPs left as instrument.


### Run MR-analysis ###

het<-mr_heterogeneity(dat11,method_list=c("mr_egger_regression", "mr_ivw"))
#id.exposure id.outcome outcome           exposure                    method        Q Q_df
#1 Age at First Birth     2mCFMo outcome Age at First Birth                  MR Egger 50.93027   39
#2 Age at First Birth     2mCFMo outcome Age at First Birth Inverse variance weighted 55.40118   40
#Q_pval
#1 0.09556459
#2 0.05340253
plt<-mr_pleiotropy_test(dat11)
#id.exposure id.outcome outcome           exposure egger_intercept           se       pval
#1 Age at First Birth     2mCFMo outcome Age at First Birth    0.0008989321 0.0004858305 0.07185676

mr_report(dat11,output_type = "html",
          output_path = "~/Downloads/MR_analysis_cwp_paper/AgeFirstBirth/AFB_to_CWP/")

afb_to_cwp <- mr(dat11, method_list =c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                       "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))

fwrite(afb_to_cwp, file="~/Downloads/MR_analysis_cwp_paper/AgeFirstBirth/AFB_to_CWP/afb_to_cwp.txt")


mr_presso <- run_mr_presso(dat11, NbDistribution = 3000, SignifThreshold = 0.05) 
#[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
#[1] 0.075


#MR.RAPS (Robust Adjusted Profile Score) is a recently proposed method that considers the 
#measurement error in SNP-exposure effects, is unbiased when there are many (e.g. hundreds of) weak instruments, 
#and is robust to systematic and idiosyncratic pleiotropy.
mr_raps <- mr_raps(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#$b -0.03175809
#$se 0.006327279
#$pval 5.187782e-07
#$nsnp 41

# mr ivw radial: for quantifying heterogeneity and outlier detection
# Raidal MR is an simple alternative to MR Presso
mr_ivw_radial <- mr_ivw_radial(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())

#F-statistic: 29.45 on 1 and 40 DF, p-value: 3.01e-06
#Q-Statistic for heterogeneity: 53.3568 on 40 DF , p-value: 0.07690085


## High Resulation PLotting 

# leave one out plot
leave_one_out <- mr_leaveoneout(dat11)
p3 <- mr_leaveoneout_plot(leave_one_out)
ggsave(p3[[1]], file="~/Downloads/MR_analysis_cwp_paper/AgeFirstBirth/AFB_to_CWP/leave_one_out_afb.pdf", width=7, height=7, dpi = 300)

# funnel plot
funnel_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                                "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p4 <- mr_funnel_plot(funnel_plot)
ggsave(p4[[1]], file="~/Downloads/MR_analysis_cwp_paper/AgeFirstBirth/AFB_to_CWP/funnel_plot_afb.pdf", width=7, height=7, dpi = 300)


# scatter plot
p1 <- mr_scatter_plot(afb_to_cwp, dat11)
length(p1)
ggsave(p1[[1]], file="~/Downloads/MR_analysis_cwp_paper/AgeFirstBirth/AFB_to_CWP/scatter_plot_mdd.pdf", width=7, height=7, dpi = 300)


# forrest plot 
forr_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                              "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p2 <- mr_forest_plot(forr_plot)
p2[[1]]
ggsave(p2[[1]], file="~/Downloads/MR_analysis_cwp_paper/AgeFirstBirth/AFB_to_CWP/forest_plot_afb.pdf", width=7, height=7, dpi = 300)

a <- directionality_test(dat11)
#id.exposure id.outcome           exposure outcome snp_r2.exposure snp_r2.outcome
#1 Age at First Birth        CWP Age at First Birth     CWP     0.004203479   0.0003847578
#correct_causal_direction steiger_pval
#1                     TRUE  7.34845e-58














