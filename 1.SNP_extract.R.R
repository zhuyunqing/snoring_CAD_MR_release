########
#By Yunqing Zhu, Canqing Yu
#Email: zhuyun_qing@126.com yucanqing@pku.edu.cn
#########

library(data.table)
library(TwoSampleMR)
library(tidyverse)
library(dbplyr)
library(dplyr)
library(LDlinkR)
library(MRPRESSO)
library(data.table)
library(ggplot2)
library(mr.raps)
library('gsmr')

setwd("C:/Users/86132/Desktop/复核")

#snoring-GWAS
#CKB-snoring  
ckb_snoring_summary <- fread("C:/Users/86132/Desktop/复核/snoring2CAD_check/data/snoring.summarydata.tsmr.txt")
ckb_snoring_summary <- ckb_snoring_summary[,c(1,2,3,4,5,6,7,8)]

snoring <- read_exposure_data(filename = "C:/Users/86132/Desktop/复核/snoring2CAD_check/data/togsmr_snoring3.txt",
                              snp_col = "SNP",
                              sep = "\t",
                              beta_col = "b",
                              se_col = "se",
                              effect_allele_col = "A1",
                              other_allele_col = "A2",
                              eaf_col = "freq",
                              pval_col = "p",
                              samplesize_col = "N"
)
head(snoring)
dim(snoring)


exp_dat_snoring <- clump_data(snoring, clump_r2 = 0.001, 
                              clump_kb = 10000, pop = "EAS")

exp_dat_snoring <- exp_dat_snoring[,-c(8:13)]
colnames(exp_dat_snoring) <- c("SNP", "EA", "OA", "EAF", "beta.snoring", "se.snoring", "P.snoring")




#habitual-GWAS
ckb_habitual_summary <- fread("C:/Users/86132/Desktop/复核/snoring2CAD_check/data/habitual.summarydata.tsmr.txt")
ckb_habitual_summary <- ckb_habitual_summary[,c(1,2,3,4,5,6,7,8)]

habitual <- read_exposure_data(filename = "C:/Users/86132/Desktop/复核/snoring2CAD_check/data/togsmr_habitual.txt",
                               snp_col = "SNP",
                               sep = "\t",
                               beta_col = "b",
                               se_col = "se",
                               effect_allele_col = "A1",
                               other_allele_col = "A2",
                               eaf_col = "freq",
                               pval_col = "p",
                               samplesize_col = "N"
)
head(habitual)
dim(habitual)

exp_dat_habitual <- clump_data(habitual, clump_r2 = 0.001, 
                               clump_kb = 10000, pop = "EAS")
exp_dat_habitual <- exp_dat_habitual[,-c(8:13)]
colnames(exp_dat_habitual) <- c("SNP", "P.habitual", "EA", "OA", "beta.habitual", "se.habitual", "EAF")


######
#BMI-GWAS
#CKB-BMI
ckb_bmi_summary <- fread("C:/Users/86132/Desktop/复核/snoring2CAD_check/data/bmi.summarydata.tsmr.txt")
ckb_bmi_summary <- ckb_bmi_summary[,c(1,2,3,4,5,6,7,8)]

BMI <- read_exposure_data(filename = "C:/Users/86132/Desktop/复核/snoring2CAD_check/data/bmi.summarydata.tsmr1.txt",
                          snp_col = "SNP",
                          sep = "\t",
                          beta_col = "beta",
                          se_col = "se",
                          effect_allele_col = "EA",
                          other_allele_col = "OA",
                          eaf_col = "EAF",
                          pval_col = "pval.exposure",
                          samplesize_col = "N",
                          phenotype_col = "exposure"
)
head(BMI)

#因为CKB的SNP有5个位点在BBJ的summary data中找不到，所以排除之，找代理位点

#"rs12885071" "rs6504568"  "rs35560038" "rs6044722"  "rs2744475"
#"rs2103785"  "rs28855509"

BMI <- BMI[which(BMI$SNP!="rs12885071" & 
                   BMI$SNP!="rs6504568" &
                   BMI$SNP!="rs35560038" &
                   BMI$SNP!="rs6044722" &
                   BMI$SNP!="rs2744475" &
                   BMI$SNP!="rs2103785" &
                   BMI$SNP!="rs28855509"),]

exp_dat_BMI <- clump_data(BMI, clump_r2 = 0.001, 
                          clump_kb = 10000, pop = "EAS")

exp_dat_BMI <- exp_dat_BMI[,-c(8:13)]
colnames(exp_dat_BMI) <- c("SNP", "P.BMI", "EA", "OA", "beta.BMI", "se.BMI", "EAF")

nrow(exp_dat_BMI)
exp_dat_BMI <- subset(exp_dat_BMI, SNP!="rs476828" & SNP!="rs6265")





#寻找另一个暴露的beta se p
#snoringSNP-BMI
exp_dat_snoring_merge <- merge(exp_dat_snoring, ckb_bmi_summary, by.x = "SNP", 
                               by.y = "SNP", all = F)
exp_dat_snoring_merge <- subset(exp_dat_snoring_merge, EA!="NA")
exp_dat_snoring_merge$BETA <- ifelse(exp_dat_snoring_merge$EA==exp_dat_snoring_merge$A1, 
                                     exp_dat_snoring_merge$BETA, -1*exp_dat_snoring_merge$BETA)
exp_dat_snoring_merge <- exp_dat_snoring_merge[,c(1,8,9,2,3,4,5,6,7,13,14,10)]
colnames(exp_dat_snoring_merge) <- c("SNP", "Chromosome", "Position", "EA", "OA", "EAF", 
                                     "BETA.snoring", "SE.snoring", "P.snoring",
                                     "BETA.BMI", "SE.BMI", "P.BMI")
exp_dat_snoring_merge$F <- round((exp_dat_snoring_merge$BETA.snoring/exp_dat_snoring_merge$SE.snoring)^2,2)
write.csv(exp_dat_snoring_merge,
          "C:/Users/86132/Desktop/复核/snoring2CAD_check/result/Supplementary Table1.3snps_snoring.csv", row.names=FALSE )


#habitualSNP-BMI
exp_dat_habitual_merge <- merge(exp_dat_habitual, ckb_bmi_summary, by.x = "SNP", 
                                by.y = "SNP", all = F)
exp_dat_habitual_merge <- subset(exp_dat_habitual_merge, EA!="NA")
exp_dat_habitual_merge$BETA <- ifelse(exp_dat_habitual_merge$EA==exp_dat_habitual_merge$A1, 
                                      exp_dat_habitual_merge$BETA, -1*exp_dat_habitual_merge$BETA)
exp_dat_habitual_merge <- exp_dat_habitual_merge[,c(1,8,9,3,4,7,5,6,2,13,14,10)]
colnames(exp_dat_habitual_merge) <- c("SNP", "Chromosome", "Position", "EA", "OA", "EAF", 
                                      "BETA.habitual", "SE.habitual", "P.habitual",
                                      "BETA.BMI", "SE.BMI", "P.BMI")
exp_dat_habitual_merge$F <- round((exp_dat_habitual_merge$BETA.habitual/exp_dat_habitual_merge$SE.habitual)^2,2)
write.csv(exp_dat_habitual_merge,
          "C:/Users/86132/Desktop/复核/snoring2CAD_check/result/Supplementary Table2.3snps_habitual.csv", row.names=FALSE )


#BMISNP-snoring-habitual
exp_dat_BMI_merge0 <- merge(exp_dat_BMI, ckb_snoring_summary, by="SNP", all=F)
exp_dat_BMI_merge0$BETA <- ifelse(exp_dat_BMI_merge0$EA==exp_dat_BMI_merge0$A1, 
                                  exp_dat_BMI_merge0$BETA, -1*exp_dat_BMI_merge0$BETA)
exp_dat_BMI_merge0 <- exp_dat_BMI_merge0[,c(1,8,9,3,4,7,5,6,2,13,14,10)]
colnames(exp_dat_BMI_merge0) <- c("SNP", "Chromosome", "Position", "EA", "OA", "EAF", 
                                  "beta.BMI", "se.BMI", "P.BMI",
                                  "beta.snoring", "se.snoring", "P.snoring")

exp_dat_BMI_merge <- merge(exp_dat_BMI_merge0, ckb_habitual_summary, by="SNP", all=F)
exp_dat_BMI_merge$BETA <- ifelse(exp_dat_BMI_merge$EA==exp_dat_BMI_merge$A1, 
                                 exp_dat_BMI_merge$BETA, -1*exp_dat_BMI_merge$BETA)
exp_dat_BMI_merge <- exp_dat_BMI_merge[,c(1:12,18,19,15)]
colnames(exp_dat_BMI_merge) <- c("SNP", "Chromosome", "Position", "EA", "OA", "EAF", 
                                 "BETA.BMI", "SE.BMI", "P.BMI",
                                 "BETA.snoring", "SE.snoring", "P.snoring",
                                 "BETA.habitual", "SE.habitual", "P.habitual")
exp_dat_BMI_merge$F <- round((exp_dat_BMI_merge$BETA.BMI/exp_dat_BMI_merge$SE.BMI)^2,2)
write.csv(exp_dat_BMI_merge,
          "C:/Users/86132/Desktop/复核/snoring2CAD_check/result/Supplementary Table3.57snps_BMI.csv", row.names=FALSE )