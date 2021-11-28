# Direction of effect MDD to CWP


# load exposure
setwd("~/Downloads/MR_analysis_cwp_paper/school_yrs_pmid_23722424/")


eu_school_dat <- read_exposure_data("school_exposure.txt",
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

eu_school_dat [,"exposure"]=eu_school_dat[,"id.exposure"]="Years of Schooling" # 77 SNPs



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
dat11 <- harmonise_data(exposure_dat = eu_school_dat, outcome_dat = cwp_out_dat) 
dat11 <- dat11[(dat11$'eaf.exposure' < 0.95) & (dat11$'eaf.exposure' > 0.05) & (dat11$'eaf.outcome' < 0.95) & (dat11$'eaf.outcome' > 0.05), ]
dat11=dat11[dat11$mr_keep==TRUE,] # 41 SNPs left as instrument.


### Run MR-analysis ###

het<-mr_heterogeneity(dat11,method_list=c("mr_egger_regression", "mr_ivw"))
#id.exposure id.outcome outcome           exposure                    method        Q Q_df
#1 Years of Schooling        CWP     CWP Years of Schooling                  MR Egger 54.09539   39
#2 Years of Schooling        CWP     CWP Years of Schooling Inverse variance weighted 57.50971   40
#Q_pval
#1 0.05463301
#2 0.03589527
plt<-mr_pleiotropy_test(dat11)
#id.exposure id.outcome outcome           exposure egger_intercept           se      pval
#1 Years of Schooling        CWP     CWP Years of Schooling     0.001390549 0.0008863022 0.1247424

mr_report(dat11,output_type = "html",
          output_path = "~/Downloads/MR_analysis_cwp_paper/school_yrs_pmid_23722424/Schooling_to_CWP/")

school_to_cwp <- mr(dat11, method_list =c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                       "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))

fwrite(school_to_cwp, file="~/Downloads/MR_analysis_cwp_paper/school_yrs_pmid_23722424/Schooling_to_CWP/school_to_cwp.txt")


mr_presso <- run_mr_presso(dat11, NbDistribution = 3000, SignifThreshold = 0.05) 
#[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
#[1] 0.04333333
#[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
#[1] "No significant outliers"


#MR.RAPS (Robust Adjusted Profile Score) is a recently proposed method that considers the 
#measurement error in SNP-exposure effects, is unbiased when there are many (e.g. hundreds of) weak instruments, 
#and is robust to systematic and idiosyncratic pleiotropy.
mr_raps <- mr_raps(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#$b -0.02099797
#$se 0.005014083
#$pval 2.816724e-05
#$nsnp 41

# mr ivw radial: for quantifying heterogeneity and outlier detection
# Raidal MR is an simple alternative to MR Presso
mr_ivw_radial <- mr_ivw_radial(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())

#F-statistic: 16.69 on 1 and 40 DF, p-value: 0.000205
#Q-Statistic for heterogeneity: 56.15118 on 40 DF , p-value: 0.04647709


## High Resulation PLotting 

# leave one out plot
leave_one_out <- mr_leaveoneout(dat11)
p3 <- mr_leaveoneout_plot(leave_one_out)
ggsave(p3[[1]], file="~/Downloads/MR_analysis_cwp_paper/school_yrs_pmid_23722424/Schooling_to_CWP/leave_one_out_school.pdf", width=7, height=7, dpi = 300)

# funnel plot
funnel_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                                "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p4 <- mr_funnel_plot(funnel_plot)
ggsave(p4[[1]], file="~/Downloads/MR_analysis_cwp_paper/school_yrs_pmid_23722424/Schooling_to_CWP/funnel_plot_school.pdf", width=7, height=7, dpi = 300)


# scatter plot
p1 <- mr_scatter_plot(school_to_cwp, dat11)
length(p1)
ggsave(p1[[1]], file="~/Downloads/MR_analysis_cwp_paper/school_yrs_pmid_23722424/Schooling_to_CWP/scatter_plot_school.pdf", width=7, height=7, dpi = 300)


# forrest plot 
forr_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                              "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p2 <- mr_forest_plot(forr_plot)
p2[[1]]
ggsave(p2[[1]], file="~/Downloads/MR_analysis_cwp_paper/school_yrs_pmid_23722424/Schooling_to_CWP/forest_plot_school.pdf", width=7, height=7, dpi = 300)

a <- directionality_test(dat11)
#id.exposure id.outcome           exposure outcome snp_r2.exposure snp_r2.outcome
#1 Years of Schooling        CWP Years of Schooling     CWP     0.009552572   0.0003263334
#correct_causal_direction  steiger_pval
#1                     TRUE 4.055314e-102





