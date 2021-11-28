
library(MRInstruments)
library(TwoSampleMR)
library(MRInstruments)
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)


# MDD GWAS without "UKB and 23&ME data"
setwd("~/Downloads/MR_analysis_cwp_paper/major_depression")
mdd <- fread("daner_pgc_mdd_meta_w2_no23andMe_rmUKBB_ PMID29700475.gz")

# No of loci in the paper: 44
# PMID: 29700475 ( used in multisite chronic pain)
mdd_ins <- mdd %>% filter(P<1e-05)
colnames(mdd_ins)
mdd_ins <- mdd_ins[,c(1:7,9:11,17,18)]
mdd_ins$FREQ <- (mdd_ins$FRQ_A_45396+mdd_ins$FRQ_U_97250)/2
mdd_ins$N <- 142646
mdd_ins$BETA <- log(mdd_ins$OR) # se and p-value not converted.
mdd_ins <- mdd_ins[,c(1:5,13,15,9,10,14)] 
mdd_ins<- mdd_ins %>% filter_all(any_vars(str_detect(., pattern = "rs"))) # 2002 obs
summary(mdd_ins)
fwrite(mdd_ins, file = "mdd_ins_e05.txt",  quote = FALSE, sep = "\t")       

## All GWAS sample was used including UKB samples: whole samole for MDD GWAS
ao <- available_outcomes()
mdddat <- extract_instruments(outcomes='ieu-a-1187', p1=5e-08, r2=0.1)
fwrite(mdddat, file = "mdd_instrument_44loci_wholesample.txt", sep = "\t")

# CWP instrument selection
setwd("~/Downloads/MR_analysis_cwp_paper/CWP_GWAS/")
cwp <- fread("all_samples_final.csv")
cwp <- cwp[,1:10]
colnames(cwp)[6] <- "FREQ" #A1FREQ to FREQ for clumping.
cwp_ins_e05 <- cwp %>% filter(P<1e-05)
cwp_ins_e05$N <- 249843
fwrite(cwp_ins_e05, file = "cwp_ins_e05.txt",  quote = FALSE, sep = "\t")



# Age at first Birth GWAS
setwd("~/Downloads/MR_analysis_cwp_paper/AgeFirstBirth")
AFB <-  fread("AgeFirstBirth_Pooled.txt")
head(AFB)
AFB$N <- 251151
#PMID: 27798627
# 10 loci (but over all with other reproducticve trats the no of loci is 12)

#AFB$BETA <- 0

#for (i in 1:nrow(AFB)) {
 # AFB$BETA[i]<- AFB[i,7]/sqrt((2*AFB[i,6]*(1-AFB[i,6]))*(AFB[i,9] + (AFB[i,7])^2))
 #}

#AFB$SE <- 0
#for (i in 1:nrow(AFB)) {
  #AFB$SE[i] <- 1/sqrt((2*AFB[,6]*(1-AFB$Freq_HapMap))*(AFB[,9]+(AFB[,7])^2)))
#}

#https://www.nature.com/articles/ng.3538
beta <- AFB[,7]/sqrt((2*AFB[,6]*(1-AFB$Freq_HapMap))*(AFB[,9] + (AFB[,7])^2)) #Beta = z/sqrt(2p(1− p)(n + z^2))
se <- 1/sqrt((2*AFB[,6]*(1-AFB$Freq_HapMap))*(AFB[,9] + (AFB[,7])^2)) # SE =1/sqrt(2p(1− p)(n + z^2))
beta_se <- cbind(beta,se)
afb <- cbind(AFB,beta_se)
colnames(afb)[c(10:11)] <- c("BETA","SE")

AFB_ctg_vl <- afb[,c(2,3,1,4,5,6,10,11,8,9)]
colnames(AFB_ctg_vl) <- c("CHR","BP","SNP","A1","A2","FREQ","BETA","SE","P","N")
fwrite(AFB_ctg_vl, file = "AFB_ctg_vl.txt", quote = FALSE, row.names=FALSE,sep = "\t")

AFB_1e05 <- AFB_ctg_vl %>% filter(P<1E-05)
fwrite(AFB_1e05, file = "AFB_1e05.txt", quote = FALSE, sep = "\t")







# Years of schooling
setwd("~/Downloads/MR_analysis_cwp_paper/school_yrs_pmid_23722424")
schooling <- fread("SSGAC_EduYears_Rietveld2013_publicrelease.txt")
# Yrs of schooling GWAS: month of schooling per allele 
# PMID: 23722424 
# 3 SNP at p < 5E-08

colnames(schooling)
yr_school_1e05 <- schooling %>% filter(Pvalue<1e-05)
yr_school_1e05$N <- 101069
colnames(yr_school_1e05) <- c("SNP","A1","A2","FREQ","BETA","SE","P","N")
fwrite(yr_school_1e05, file = "yr_school_1e05.txt", quote = FALSE, sep = "\t")




## Clumping for all of these were done using plink separately.

setwd("~/Downloads/MR_analysis_cwp_paper")

mdd_clump <- fread("mdd_clumped_snp_1e05.txt") #67 snps
mdd_clump <- data.frame(mdd_clump)

cwp_clump <- fread("cwp_clumped_snp_1e05.txt", h=F) #39 snps
str(cwp_clump)
cwp_clump <- data.frame(cwp_clump)
c

afb_clump <- fread("afb_clumped_snp_1e05.txt") #52 snps
afb_clump <- data.frame(afb_clump)

school_clump <- fread("schooling_clumped_snp_1e05.txt") #67 snps
school_clump <- data.frame(school_clump)

#Creating ready to use data for MR base: 
setwd("~/Downloads/MR_analysis_cwp_paper/major_depression")
mdd_ins <- fread("mdd_ins_e05.txt")
mdd_exposure <- mdd_ins[which(mdd_ins$SNP%in%mdd_clump$SNP),]
fwrite(mdd_exposure, file = "mdd_exposure.txt", sep = "\t", quote = FALSE)


setwd("~/Downloads/MR_analysis_cwp_paper/CWP_GWAS")


cwp<- fread("all_samples_final.csv")

cwp<- cwp %>% filter_all(any_vars(str_detect(., pattern = "rs")))
str(cwp)
cwp_clump <- fread("cwp_ins_clump.clumped")
cwp_clump  <- cwp_clump [,3]
cwp_exposure <- cwp[which(cwp$SNP %in% cwp_clump$SNP),]
cwp_exposure <- cwp_exposure[!duplicated(cwp_exposure$SNP),]
cwp_exposure <- cwp_exposure %>% filter(P<0.0001)

fwrite(cwp_exposure, file = "cwp_exposure.txt", sep = "\t")





setwd("~/Downloads/MR_analysis_cwp_paper/AgeFirstBirth")
AFB_1e05 <- fread("AFB_1e05.txt")
afb_exposure <- AFB_1e05[which(AFB_1e05$SNP%in%afb_clump$SNP),]
fwrite(afb_exposure, file = "afb_exposure.txt", sep = "\t", quote = FALSE)


setwd("~/Downloads/MR_analysis_cwp_paper/school_yrs_pmid_23722424")
yr_school_1e05 <- fread("yr_school_1e05.txt")
school_exposure <- yr_school_1e05[which(yr_school_1e05$SNP%in%school_clump$SNP),]
fwrite(school_exposure, file = "school_exposure.txt", sep = "\t", quote = FALSE)






## BMI

setwd("~/Downloads/MR_analysis_cwp_paper/bmi")
library(readxl)

# 97 snps extractiin from all ancestry data
all_bmi <- fread("all_ancestry_bmi_2015.gz") # 2555510
a <- read_excel("locke_etal_2015.xls", sheet = 5)
all_bmi_ins <- all_bmi[all_bmi$SNP %in% a$SNP,]
fwrite(all_bmi_ins, file = "all_bmi_ins.txt", sep = "\t")


# 77 snps extraction from european ancestry data
eu_bmi <- fread("exposure_BMI/eu_ancestry_bmi_2015.gz") # 2554637
eu_snp_ins <- read_excel("locke_etal_2015.xls", sheet=4)
eu_snp_ins <- eu_bmi[eu_bmi$SNP %in% eu_snp_ins$SNP, ]
eu_snp_ins <- eu_snp_ins [eu_snp_ins$p < 5e-08,]
fwrite(eu_snp_ins, file="eu_snp_ins.txt", sep = "\t")

# replicated snps extraction from european ancestry data
nobel_snps_locke <- read_excel("locke_etal_2015.xls", sheet = 3) #56
replicated_snps_locke <- all_bmi_ins[!all_bmi_ins$SNP %in% nobel_snps_locke$SNP,]
mean(replicated_snps_locke$p)
eu_repli_snps_ins <- eu_bmi[eu_bmi$SNP %in% replicated_snps_locke$SNP,]
fwrite(eu_repli_snps_ins, file = "eu_repli_snps_ins.txt", sep = "\t") #41











