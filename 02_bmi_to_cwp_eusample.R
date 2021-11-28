
#install.packages("devtools")


#devtools::install_github("MRCIEU/TwoSampleMR")

#devtools::install_github("MRCIEU/MRInstruments")

#install.packages("Cairo")

library(TwoSampleMR)
library(MRInstruments)
library(data.table)
library(ggplot2)


setwd("~/Downloads/MR_analysis_cwp_paper/bmi")
# MR analysis with all 77 BMI instrument selected from EU ancestry data set:

eu_bmi_dat <- read_exposure_data("eu_snp_ins.txt",
                                  clump = FALSE,
                                  sep = "\t",
                                  snp_col = "SNP",
                                  beta_col = "b",
                                  se_col = "se", 
                                  eaf_col = "Freq1.Hapmap",
                                  effect_allele_col = "A1",
                                  other_allele_col = "A2", 
                                  pval_col = "p",
                                  samplesize_col = "N")

eu_bmi_dat [,"exposure"]=eu_bmi_dat[,"id.exposure"]="BMI" # 77 SNPs

## instrument stength approximation: 

#eu_bmi_dat$f_stat <- ((eu_bmi_dat$beta.exposure)^2)/((eu_bmi_dat$se.exposure)^2)
#a <- round(eu_bmi_dat[,14], digits = 2) #64.24
#mean(a)
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


# Harmonization: Exposure 1: BMI, outcome 1: Back pain
dat11 <- harmonise_data(exposure_dat = eu_bmi_dat, outcome_dat = cwp_out_dat) 
dat11 <- dat11[(dat11$'eaf.exposure' < 0.95) & (dat11$'eaf.exposure' > 0.05) & (dat11$'eaf.outcome' < 0.95) & (dat11$'eaf.outcome' > 0.05), ]
dat11=dat11[dat11$mr_keep==TRUE,] # 67 SNPs left as instrument.

dim(dat11)



#dat11$f_stat <- ((dat11$beta.exposure)^2)/((dat11$se.exposure)^2)
#a <- dat11[,31]
#mean(a) #58.03


### Run MR-analysis ###

het<-mr_heterogeneity(dat11,method_list=c("mr_egger_regression", "mr_ivw"))
#id.exposure id.outcome outcome exposure                    method        Q Q_df       Q_pval
#1         BMI        CWP     CWP      BMI                  MR Egger 115.2426   65 0.0001242564
#2         BMI        CWP     CWP      BMI Inverse variance weighted 116.3887   66 0.0001293945
plt<-mr_pleiotropy_test(dat11)
#id.exposure id.outcome outcome exposure egger_intercept           se      pval
#1         BMI        CWP     CWP      BMI    0.0002036849 0.0002533358 0.4243212
#mr_method_list() 

mr_report(dat11,output_type = "html",
          output_path = "~/Downloads/MR_analysis_cwp_paper/bmi/77_snp_no_clumped/")

bmi_to_cwp <- mr(dat11, method_list =c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe","mr_weighted_median","mr_weighted_mode"))
fwrite(bmi_to_cwp, file="~/Downloads/MR_analysis_cwp_paper/bmi/77_snp_no_clumped/bmi_to_cwp.csv")

set.seed(10000)
mr_presso <- run_mr_presso(dat11, NbDistribution = 3000, SignifThreshold = 0.05) # no significant outliers detected.
#[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
#[1] 0.001

#[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
#[1] 38

#[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
#[1] 0.7603333



# After elimination of rs205262 as pleiotrophic outlier
#dat13 <- dat11[-38,] # excluding the outliers 
#set.seed(10000)
#mr_presso_corrected <- run_mr_presso(dat13, NbDistribution = 3000, SignifThreshold = 0.05)

#bmi_to_cwp_outlier_corrected <- mr(dat13, method_list =c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
#                                       "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))





#MR.RAPS (Robust Adjusted Profile Score) is a recently proposed method that considers the 
#measurement error in SNP-exposure effects, is unbiased when there are many (e.g. hundreds of) weak instruments, 
#and is robust to systematic and idiosyncratic pleiotropy.

mr_raps <- mr_raps(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#b 0.0143823
#se 0.003211219
#p 7.507532e-06
#snps 67

# mr ivw radial: for quantifying heterogeneity and outlier detection
# Raidal MR is an simple alternative to MR Presso
mr_ivw_radial <- mr_ivw_radial(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#F-statistic: 20.2 on 1 and 66 DF, p-value: 2.88e-05
#Q-Statistic for heterogeneity: 115.3469 on 66 DF , p-value: 0.0001649589


# https://mrcieu.github.io/TwoSampleMR/


## High Resulation PLotting 

# leave one out plot
leave_one_out <- mr_leaveoneout(dat11)
p3 <- mr_leaveoneout_plot(leave_one_out)
ggsave(p3[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/77_snp_no_clumped/leave_one_out_77noclump.pdf", width=7, height=7, dpi = 300)

# funnel plot
funnel_plot <- mr_singlesnp(dat11, all_method =c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                                  "mr_weighted_median","mr_weighted_mode"))
p4 <- mr_funnel_plot(funnel_plot)
ggsave(p4[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/77_snp_no_clumped/funnel_plot_77noclump.pdf", width=7, height=7, dpi = 300)


# scatter plot
p1 <- mr_scatter_plot(bmi_to_cwp, dat11)
length(p1)[1]
ggsave(p1[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/77_snp_no_clumped/scatter_plot_77noclump.pdf", width=7, height=7, dpi = 300)


# forrest plot 
forr_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                              "mr_weighted_median","mr_weighted_mode"))
p2 <- mr_forest_plot(forr_plot)
p2[[1]]
ggsave(p2[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/77_snp_no_clumped/forest_plot_77noclump.pdf", width=7, height=7, dpi = 300)



colnames(dat11)
dat11$r2_bmi <- 0
## varriance explained by each 81 instruments in CWP:
for (i in 1:nrow(dat11)){
  dat11$r2_bmi[i] <- (2*(dat11[i,6]^2)*(dat11[i,8])*(1-dat11[i,8]))/((2*(dat11[i,6]^2)*(dat11[i,8])) + (dat11[i,21]^2)*(2*dat11[i,23])*dat11[i,8]*(1-dat11[i,8]))
}
sum(dat11$r2_bmi)*100


## varriance explained by each 81 instruments in CWP:

# Formulla: F <- (n-k-1/k) * (R^2/(1-R^2))

F <- ((249843-67-1)/67)*((sum(dat11$r2_bmi))^2)/(1-(sum(dat11$r2_bmi))^2)

F*100


0.0006083*100












#####################################################################################################################
            
                  # Replicated BMI SNPs were used withput clumpting #

###################################################################################################################
setwd("~/Downloads/MR_analysis_cwp_paper/bmi")

eu_bmi_dat <- read_exposure_data("eu_repli_snps_ins.txt",    #41 formerly associated SNPs with obesity relarted traits in locke et al 2015
                                 clump = FALSE,
                                 sep = "\t",
                                 snp_col = "SNP",
                                 beta_col = "b",
                                 se_col = "se", 
                                 eaf_col = "Freq1.Hapmap",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2", 
                                 pval_col = "p",
                                 samplesize_col = "N")

eu_bmi_dat [,"exposure"]=eu_bmi_dat[,"id.exposure"]="BMI"


# Harmonization: Exposure 1: BMI, outcome 1: Back pain
dat11 <- harmonise_data(exposure_dat = eu_bmi_dat, outcome_dat = cwp_out_dat) # 39 snps remained, 2 SNP removed becoz paliandromic

# filtering based on effect allele frequency: 
dat11 <- dat11[(dat11$'eaf.exposure' < 0.95) & (dat11$'eaf.exposure' > 0.05) & (dat11$'eaf.outcome' < 0.95) & (dat11$'eaf.outcome' > 0.05), ]

dat11 = dat11[dat11$mr_keep==TRUE,] # 36 SNPs left as instrument.


a <- run_mr_presso(dat11) # no significant outliers detected.


mr_report(dat11,output_type = "html",
          output_path = "~/Downloads/MR_analysis_cwp_paper/bmi/41_replicated_snp_no_clump/")

bmi_to_cwp_replicationSNPS <- mr(dat11, method_list =c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                                       "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
fwrite(bmi_to_cwp_replicationSNPS, file="~/Downloads/MR_analysis_cwp_paper/bmi/41_replicated_snp_no_clump/bmi_to_cwp_replicationSNPS.txt")


mr_ivw_radial <- mr_ivw_radial(dat11$beta.exposure, dat11$beta.outcome, dat11$se.exposure, dat11$se.outcome, parameters = default_parameters())
#F-statistic: 16.33 on 1 and 35 DF, p-value: 0.000277
#Q-Statistic for heterogeneity: 60.63158 on 35 DF , p-value: 0.004585004

## High Resulation PLotting 

# leave one out plot
leave_one_out <- mr_leaveoneout(dat11)
p3 <- mr_leaveoneout_plot(leave_one_out)
ggsave(p3[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/41_replicated_snp_no_clump/leave_one_out_41_noclumped.pdf", width=7, height=7, dpi = 300)

# funnel plot
funnel_plot <- mr_singlesnp(dat11, all_method =c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                                  "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p4 <- mr_funnel_plot(funnel_plot)
ggsave(p4[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/41_replicated_snp_no_clump/funnel_plot_41_noclumped.pdf", width=7, height=7, dpi = 300)


# scatter plot
p1 <- mr_scatter_plot(bmi_to_cwp_replicationSNPS, dat11)
length(p1)
ggsave(p1[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/41_replicated_snp_no_clump/scatter_plot_41_noclumped.pdf", width=7, height=7, dpi = 300)


# forrest plot 
forr_plot <- mr_singlesnp(dat11, all_method=c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe",
                                              "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
p2 <- mr_forest_plot(forr_plot)
p2[[1]]
ggsave(p2[[1]], file="~/Downloads/MR_analysis_cwp_paper/bmi/41_replicated_snp_no_clump/forest_plot_41_noclumped.pdf", width=7, height=7, dpi = 300)

###################################################################################################################



## to calculate F-statistics, use formula available in this article 
#http://www.phpc.cam.ac.uk/ceu/files/2012/11/paper3four131110.pdf





