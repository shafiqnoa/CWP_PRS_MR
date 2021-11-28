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






