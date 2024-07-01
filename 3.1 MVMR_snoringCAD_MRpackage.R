########
#By Yunqing Zhu, Canqing Yu
#Email: zhuyun_qing@126.com yucanqing@pku.edu.cn
#########

library(data.table)
library(TwoSampleMR)
library(MendelianRandomization)
library(tidyverse)
library(dbplyr)
library(dplyr)
library(LDlinkR)
library(MRPRESSO)
library(data.table)
library(xtable)
library(flextable)
library(officer)
library('gsmr')

#CKB-BMI
ckb_bmi_summary <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/bmi.summarydata.tsmr.txt")
ckb_bmi_summary1 <- ckb_bmi_summary[,c(1,5,4,7,8)]
setnames(ckb_bmi_summary1, old = c("SNP","A1","P","BETA","SE"), 
         new = c("SNP","ea","pval.exposure","beta","se"))
ckb_bmi_summary1$exposure <- 'bmi'

#BMI-IVs
ckb_bmi_iv0 <- subset(ckb_bmi_summary1, pval.exposure<5e-08 )
ckb_bmi_iv1 <- clump_data(ckb_bmi_iv0, clump_r2 = 0.001, 
                          clump_kb = 10000, pop = "EAS")
ckb_bmi_iv1 <- ckb_bmi_iv1[,-7]



#CKB-snoring  
ckb_snoring_summary <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/snoring.summarydata.tsmr.txt")
ckb_snoring_summary1 <- ckb_snoring_summary[,c(1,5,4,7,8)]
setnames(ckb_snoring_summary1, old = c("SNP","A1","P","BETA","SE"), 
         new = c("SNP","ea","pval.exposure","beta","se"))
ckb_snoring_summary1$exposure <- 'snoring'

#CKB-snoring-IVs
ckb_snoring_iv0 <- subset(ckb_snoring_summary1, pval.exposure<5e-08 )
ckb_snoring_iv1 <- clump_data(ckb_snoring_iv0, clump_r2 = 0.001, 
                              clump_kb = 10000, pop = "EAS")
ckb_snoring_iv1 <- ckb_snoring_iv1[,-7]


#CKB-habitual
ckb_habitual_summary <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/habitual.summarydata.tsmr.txt")
ckb_habitual_summary1 <- ckb_habitual_summary[,c(1,5,4,7,8)]
setnames(ckb_habitual_summary1, old = c("SNP","A1","P","BETA","SE"), 
         new = c("SNP","ea","pval.exposure","beta","se"))
ckb_habitual_summary1$exposure <- 'habitual'

#CKB-habitual-IVs
ckb_habitual_iv0 <- subset(ckb_habitual_summary1, pval.exposure<5e-08 )
ckb_habitual_iv1 <- clump_data(ckb_habitual_iv0, clump_r2 = 0.001, 
                               clump_kb = 10000, pop = "EAS")
ckb_habitual_iv1 <- ckb_habitual_iv1[,-7]


#×¢?â£ºBBJ??A1 A2Êµ??????Ref ?? ALT??????Êµ????Òª????A2??????A1
#????CAD???Ñ¾??????ÃµÄ£??????Ù»???

#load outcome summary data
bbj_CAD_summary <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/CAD_BBJ_tsmr.txt")
bbj_CAD_summary1 <- bbj_CAD_summary[,c(1,5,9,7,8)]
setnames(bbj_CAD_summary1, old = c("SNP","A1","P","b","se"),
         new = c("SNP","ea","pval.exposure","beta","se"))

bbj_MI_summary <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/MI_BBJ_tsmr.txt")
bbj_MI_summary1 <- bbj_MI_summary[,c(1,6,10,8,9)]
setnames(bbj_MI_summary1, old = c("SNP","A2","P","b","se"),
         new = c("SNP","ea","pval.exposure","beta","se"))

bbj_CHF_summary <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/CHF_BBJ_tsmr.txt")
bbj_CHF_summary1 <- bbj_CHF_summary[,c(1,6,10,8,9)]
setnames(bbj_CHF_summary1, old = c("SNP","A2","P","b","se"),
         new = c("SNP","ea","pval.exposure","beta","se"))


bbj_Angina_summary <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/Angina_BBJ_tsmr.txt")
bbj_Angina_summary1 <- bbj_Angina_summary[,c(1,6,10,8,9)]
setnames(bbj_Angina_summary1, old = c("SNP","ALT","P","b","se"),
         new = c("SNP","ea","pval.exposure","beta","se"))


bbj_UAP_summary <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/UAP_BBJ_tsmr.txt")
bbj_UAP_summary1 <- bbj_UAP_summary[,c(1,6,10,8,9)]
setnames(bbj_UAP_summary1, old = c("SNP","ALT","P","b","se"),
         new = c("SNP","ea","pval.exposure","beta","se"))


bbj_SAP_summary <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/SAP_BBJ_tsmr.txt")
bbj_SAP_summary1 <- bbj_SAP_summary[,c(1,6,10,8,9)]
setnames(bbj_SAP_summary1, old = c("SNP","ALT","P","b","se"),
         new = c("SNP","ea","pval.exposure","beta","se"))







#check summary data
#SNP, ea (in capitals), beta, se, pval.expsoure
head(ckb_snoring_summary1)
head(ckb_habitual_summary1)
head(ckb_bmi_summary1)
head(bbj_CAD_summary1)


#load instruments
#change_heading
ckb_bmi_iv <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/ckb_bmi_iv.txt")
#BMI??Î»????Òª?Å³?2????????Ç¿??Áª?Äµ?
ckb_bmi_iv <- subset(ckb_bmi_iv, SNP!="rs476828" & SNP!="rs6265")

ckb_snoring_iv <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/ckb_snoring_iv.txt")
ckb_habitual_iv <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/ckb_habitual_iv.txt")



#############################
#ALL IVs-snoring
#############################

#joint snps from each exposure
joint_snoring <- as.data.frame(c(ckb_snoring_iv$SNP, ckb_bmi_iv$SNP))
colnames(joint_snoring)[1]<-'SNP'
joint_snoring <- clump_data(joint_snoring, clump_r2 = 0.001, 
                            clump_kb = 10000, pop = "EAS")


joint_habitual <- as.data.frame(c(ckb_habitual_iv$SNP, ckb_bmi_iv$SNP))
colnames(joint_habitual)[1]<-'SNP'
joint_habitual <- clump_data(joint_habitual, clump_r2 = 0.001, 
                             clump_kb = 10000, pop = "EAS")

#CKB-BMI??4??Î»????BBJ???Ò²?????
#"rs1532127" "rs6044722" "rs4973506" "rs2744475"
# snps59 <- ckb_snoring_summary1$SNP[which(paste(ckb_snoring_summary1$SNP) %in% paste(snoring_bmi_T2D_joint_gx$SNP))]
# snps55 <- bbj_T2D_summary1$SNP[which(paste(bbj_T2D_summary1$SNP) %in% paste(snoring_bmi_T2D_joint_gx$SNP))]
# snps59[which(!snps59 %in% snps55)]

#Ñ°?Ò´???Î»??
#rs1532127-rs55731973??r2=0.96??
#rs6044722-rs1361511??r2=0.98??
#rs4973506-rs61146033??r2=0.99??
#rs2744475-rs2206271 ??r2=1??

#snoring??Î»??rs8043757&rs11642015Á¬??-????Ç°??

proxy_bmi <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/bmi_proxy.txt",h=T)
snoring_1snp <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/snoring_1snp.txt",h=T)




#?Ãµ????Õµ?joint-snps
joint_snoring1 <- joint_snoring[which(joint_snoring$SNP!="rs1532127" & 
                                        joint_snoring$SNP!="rs6044722" &
                                        joint_snoring$SNP!="rs4973506" &
                                        joint_snoring$SNP!="rs2744475" &
                                        joint_snoring$SNP!="rs11642015"),]
joint_snoring2 <- rbind(joint_snoring1, proxy_bmi, snoring_1snp)
joint_snoring = joint_snoring2








#MI
#extract SNP infomation from GWAS
snoring_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                   paste(joint_snoring$SNP))]
snoring_bmi_MI_joint_gx <- ckb_snoring_summary1[ which( paste(ckb_snoring_summary1$SNP) %in% 
                                                          paste(snoring_bmi_joint_ga$SNP)),]
snoring_bmi_MI_joint_gx <- snoring_bmi_MI_joint_gx[!duplicated(snoring_bmi_MI_joint_gx$SNP),]
snoring_bmi_MI_joint_gy <- bbj_MI_summary1[ which( paste(bbj_MI_summary1$SNP) %in% 
                                                     paste(snoring_bmi_MI_joint_gx$SNP))]
snoring_bmi_MI_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                    %in% paste(snoring_bmi_MI_joint_gx$SNP)),]
#arrange
snoring_bmi_MI_joint_gx <- arrange(snoring_bmi_MI_joint_gx, SNP)
snoring_bmi_MI_joint_gy <- arrange(snoring_bmi_MI_joint_gy, SNP)
snoring_bmi_MI_joint_ga <- arrange(snoring_bmi_MI_joint_ga, SNP)

#????Á´beta
snoring_bmi_MI_joint_gy$gy <- ifelse(snoring_bmi_MI_joint_gy$ea == snoring_bmi_MI_joint_gx$ea, 
                                     snoring_bmi_MI_joint_gy$beta,
                                     -1*snoring_bmi_MI_joint_gy$beta)
snoring_bmi_MI_joint_ga$ga <- ifelse(snoring_bmi_MI_joint_ga$ea == snoring_bmi_MI_joint_gx$ea,
                                     snoring_bmi_MI_joint_ga$beta,
                                     -1*snoring_bmi_MI_joint_ga$beta)

#MVMR???Ý¿?
snoring_bmi_MI_joint <- data.frame("V1"=1:57)
snoring_bmi_MI_joint_gx$SNP <- snoring_bmi_MI_joint_gx$SNP
snoring_bmi_MI_joint$ea <- snoring_bmi_MI_joint_gx$ea
snoring_bmi_MI_joint$gx <- snoring_bmi_MI_joint_gx$beta
snoring_bmi_MI_joint$gx_se <- snoring_bmi_MI_joint_gx$se
snoring_bmi_MI_joint$gy <- snoring_bmi_MI_joint_gy$gy
snoring_bmi_MI_joint$gy_se <- snoring_bmi_MI_joint_gy$se
snoring_bmi_MI_joint$ga <- snoring_bmi_MI_joint_ga$ga
snoring_bmi_MI_joint$ga_se <- snoring_bmi_MI_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(snoring_bmi_MI_joint$gx, snoring_bmi_MI_joint$ga), 
                           bxse = cbind(snoring_bmi_MI_joint$gx_se, snoring_bmi_MI_joint$ga_se),
                           by = snoring_bmi_MI_joint$gy,
                           byse = snoring_bmi_MI_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(snoring_bmi_MI_joint$gx, snoring_bmi_MI_joint$ga), 
                                   bxse = cbind(snoring_bmi_MI_joint$gx_se, snoring_bmi_MI_joint$ga_se),
                                   by = snoring_bmi_MI_joint$gy,
                                   byse = snoring_bmi_MI_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(snoring_bmi_MI_joint$gx, snoring_bmi_MI_joint$ga), 
                               bxse = cbind(snoring_bmi_MI_joint$gx_se, snoring_bmi_MI_joint$ga_se),
                               by = snoring_bmi_MI_joint$gy,
                               byse = snoring_bmi_MI_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(snoring_bmi_MI_joint$gx, snoring_bmi_MI_joint$ga), 
                               bxse = cbind(snoring_bmi_MI_joint$gx_se, snoring_bmi_MI_joint$ga_se),
                               by = snoring_bmi_MI_joint$gy,
                               byse = snoring_bmi_MI_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Snoring", "BMI", "Snoring", "BMI")
table$Outcome <- c("MI","MI","MI","MI")
table0 = table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-snoring-MI")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test0 = test






##CHF
#extract SNP infomation from GWAS
snoring_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                   paste(joint_snoring$SNP))]
snoring_bmi_CHF_joint_gx <- ckb_snoring_summary1[ which( paste(ckb_snoring_summary1$SNP) %in% 
                                                           paste(snoring_bmi_joint_ga$SNP)),]
snoring_bmi_CHF_joint_gx <- snoring_bmi_CHF_joint_gx[!duplicated(snoring_bmi_CHF_joint_gx$SNP),]
snoring_bmi_CHF_joint_gy <- bbj_CHF_summary1[ which( paste(bbj_CHF_summary1$SNP) %in% 
                                                       paste(snoring_bmi_CHF_joint_gx$SNP)),]
snoring_bmi_CHF_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                     %in% paste(snoring_bmi_CHF_joint_gx$SNP)),]
#arrange
snoring_bmi_CHF_joint_gx <- arrange(snoring_bmi_CHF_joint_gx, SNP)
snoring_bmi_CHF_joint_gy <- arrange(snoring_bmi_CHF_joint_gy, SNP)
snoring_bmi_CHF_joint_ga <- arrange(snoring_bmi_CHF_joint_ga, SNP)

#????Á´beta
snoring_bmi_CHF_joint_gy$gy <- ifelse(snoring_bmi_CHF_joint_gy$ea == snoring_bmi_CHF_joint_gx$ea, 
                                      snoring_bmi_CHF_joint_gy$beta,
                                      -1*snoring_bmi_CHF_joint_gy$beta)
snoring_bmi_CHF_joint_ga$ga <- ifelse(snoring_bmi_CHF_joint_ga$ea == snoring_bmi_CHF_joint_gx$ea,
                                      snoring_bmi_CHF_joint_ga$beta,
                                      -1*snoring_bmi_CHF_joint_ga$beta)

#MVMR???Ý¿?
snoring_bmi_CHF_joint <- data.frame("V1"=1:57)
snoring_bmi_CHF_joint_gx$SNP <- snoring_bmi_CHF_joint_gx$SNP
snoring_bmi_CHF_joint$ea <- snoring_bmi_CHF_joint_gx$ea
snoring_bmi_CHF_joint$gx <- snoring_bmi_CHF_joint_gx$beta
snoring_bmi_CHF_joint$gx_se <- snoring_bmi_CHF_joint_gx$se
snoring_bmi_CHF_joint$gy <- snoring_bmi_CHF_joint_gy$gy
snoring_bmi_CHF_joint$gy_se <- snoring_bmi_CHF_joint_gy$se
snoring_bmi_CHF_joint$ga <- snoring_bmi_CHF_joint_ga$ga
snoring_bmi_CHF_joint$ga_se <- snoring_bmi_CHF_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(snoring_bmi_CHF_joint$gx, snoring_bmi_CHF_joint$ga), 
                           bxse = cbind(snoring_bmi_CHF_joint$gx_se, snoring_bmi_CHF_joint$ga_se),
                           by = snoring_bmi_CHF_joint$gy,
                           byse = snoring_bmi_CHF_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(snoring_bmi_CHF_joint$gx, snoring_bmi_CHF_joint$ga), 
                                   bxse = cbind(snoring_bmi_CHF_joint$gx_se, snoring_bmi_CHF_joint$ga_se),
                                   by = snoring_bmi_CHF_joint$gy,
                                   byse = snoring_bmi_CHF_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(snoring_bmi_CHF_joint$gx, snoring_bmi_CHF_joint$ga), 
                               bxse = cbind(snoring_bmi_CHF_joint$gx_se, snoring_bmi_CHF_joint$ga_se),
                               by = snoring_bmi_CHF_joint$gy,
                               byse = snoring_bmi_CHF_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(snoring_bmi_CHF_joint$gx, snoring_bmi_CHF_joint$ga), 
                               bxse = cbind(snoring_bmi_CHF_joint$gx_se, snoring_bmi_CHF_joint$ga_se),
                               by = snoring_bmi_CHF_joint$gy,
                               byse = snoring_bmi_CHF_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Snoring", "BMI", "Snoring", "BMI")
table$Outcome <- c("CHF","CHF","CHF","CHF")
table <- rbind(table0, table)
table0=table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-snoring-CHF")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test <- rbind(test0, test)
test0 = test





##Angina
#extract SNP infomation from GWAS
snoring_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                   paste(joint_snoring$SNP))]
snoring_bmi_Angina_joint_gx <- ckb_snoring_summary1[ which( paste(ckb_snoring_summary1$SNP) %in% 
                                                              paste(snoring_bmi_joint_ga$SNP)),]
snoring_bmi_Angina_joint_gx <- snoring_bmi_Angina_joint_gx[!duplicated(snoring_bmi_Angina_joint_gx$SNP),]
snoring_bmi_Angina_joint_gy <- bbj_Angina_summary1[ which( paste(bbj_Angina_summary1$SNP) %in% 
                                                             paste(snoring_bmi_Angina_joint_gx$SNP)),]
snoring_bmi_Angina_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                        %in% paste(snoring_bmi_Angina_joint_gx$SNP)),]
#arrange
snoring_bmi_Angina_joint_gx <- arrange(snoring_bmi_Angina_joint_gx, SNP)
snoring_bmi_Angina_joint_gy <- arrange(snoring_bmi_Angina_joint_gy, SNP)
snoring_bmi_Angina_joint_ga <- arrange(snoring_bmi_Angina_joint_ga, SNP)

#????Á´beta
snoring_bmi_Angina_joint_gy$gy <- ifelse(snoring_bmi_Angina_joint_gy$ea == snoring_bmi_Angina_joint_gx$ea, 
                                         snoring_bmi_Angina_joint_gy$beta,
                                         -1*snoring_bmi_Angina_joint_gy$beta)
snoring_bmi_Angina_joint_ga$ga <- ifelse(snoring_bmi_Angina_joint_ga$ea == snoring_bmi_Angina_joint_gx$ea,
                                         snoring_bmi_Angina_joint_ga$beta,
                                         -1*snoring_bmi_Angina_joint_ga$beta)

#MVMR???Ý¿?
snoring_bmi_Angina_joint <- data.frame("V1"=1:57)
snoring_bmi_Angina_joint_gx$SNP <- snoring_bmi_Angina_joint_gx$SNP
snoring_bmi_Angina_joint$ea <- snoring_bmi_Angina_joint_gx$ea
snoring_bmi_Angina_joint$gx <- snoring_bmi_Angina_joint_gx$beta
snoring_bmi_Angina_joint$gx_se <- snoring_bmi_Angina_joint_gx$se
snoring_bmi_Angina_joint$gy <- snoring_bmi_Angina_joint_gy$gy
snoring_bmi_Angina_joint$gy_se <- snoring_bmi_Angina_joint_gy$se
snoring_bmi_Angina_joint$ga <- snoring_bmi_Angina_joint_ga$ga
snoring_bmi_Angina_joint$ga_se <- snoring_bmi_Angina_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(snoring_bmi_Angina_joint$gx, snoring_bmi_Angina_joint$ga), 
                           bxse = cbind(snoring_bmi_Angina_joint$gx_se, snoring_bmi_Angina_joint$ga_se),
                           by = snoring_bmi_Angina_joint$gy,
                           byse = snoring_bmi_Angina_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(snoring_bmi_Angina_joint$gx, snoring_bmi_Angina_joint$ga), 
                                   bxse = cbind(snoring_bmi_Angina_joint$gx_se, snoring_bmi_Angina_joint$ga_se),
                                   by = snoring_bmi_Angina_joint$gy,
                                   byse = snoring_bmi_Angina_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(snoring_bmi_Angina_joint$gx, snoring_bmi_Angina_joint$ga), 
                               bxse = cbind(snoring_bmi_Angina_joint$gx_se, snoring_bmi_Angina_joint$ga_se),
                               by = snoring_bmi_Angina_joint$gy,
                               byse = snoring_bmi_Angina_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(snoring_bmi_Angina_joint$gx, snoring_bmi_Angina_joint$ga), 
                               bxse = cbind(snoring_bmi_Angina_joint$gx_se, snoring_bmi_Angina_joint$ga_se),
                               by = snoring_bmi_Angina_joint$gy,
                               byse = snoring_bmi_Angina_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Snoring", "BMI", "Snoring", "BMI")
table$Outcome <- c("Angina","Angina","Angina","Angina")
table <- rbind(table0, table)
table0 = table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-snoring-Angina")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test <- rbind(test0, test)
test0 = test







##SAP
#extract SNP infomation from GWAS
snoring_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                   paste(joint_snoring$SNP))]
snoring_bmi_SAP_joint_gx <- ckb_snoring_summary1[ which( paste(ckb_snoring_summary1$SNP) %in% 
                                                           paste(snoring_bmi_joint_ga$SNP)),]
snoring_bmi_SAP_joint_gx <- snoring_bmi_SAP_joint_gx[!duplicated(snoring_bmi_SAP_joint_gx$SNP),]
snoring_bmi_SAP_joint_gy <- bbj_SAP_summary1[ which( paste(bbj_SAP_summary1$SNP) %in% 
                                                       paste(snoring_bmi_SAP_joint_gx$SNP)),]
snoring_bmi_SAP_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                     %in% paste(snoring_bmi_SAP_joint_gx$SNP)),]
#arrange
snoring_bmi_SAP_joint_gx <- arrange(snoring_bmi_SAP_joint_gx, SNP)
snoring_bmi_SAP_joint_gy <- arrange(snoring_bmi_SAP_joint_gy, SNP)
snoring_bmi_SAP_joint_ga <- arrange(snoring_bmi_SAP_joint_ga, SNP)

#????Á´beta
snoring_bmi_SAP_joint_gy$gy <- ifelse(snoring_bmi_SAP_joint_gy$ea == snoring_bmi_SAP_joint_gx$ea, 
                                      snoring_bmi_SAP_joint_gy$beta,
                                      -1*snoring_bmi_SAP_joint_gy$beta)
snoring_bmi_SAP_joint_ga$ga <- ifelse(snoring_bmi_SAP_joint_ga$ea == snoring_bmi_SAP_joint_gx$ea,
                                      snoring_bmi_SAP_joint_ga$beta,
                                      -1*snoring_bmi_SAP_joint_ga$beta)

#MVMR???Ý¿?
snoring_bmi_SAP_joint <- data.frame("V1"=1:57)
snoring_bmi_SAP_joint_gx$SNP <- snoring_bmi_SAP_joint_gx$SNP
snoring_bmi_SAP_joint$ea <- snoring_bmi_SAP_joint_gx$ea
snoring_bmi_SAP_joint$gx <- snoring_bmi_SAP_joint_gx$beta
snoring_bmi_SAP_joint$gx_se <- snoring_bmi_SAP_joint_gx$se
snoring_bmi_SAP_joint$gy <- snoring_bmi_SAP_joint_gy$gy
snoring_bmi_SAP_joint$gy_se <- snoring_bmi_SAP_joint_gy$se
snoring_bmi_SAP_joint$ga <- snoring_bmi_SAP_joint_ga$ga
snoring_bmi_SAP_joint$ga_se <- snoring_bmi_SAP_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(snoring_bmi_SAP_joint$gx, snoring_bmi_SAP_joint$ga), 
                           bxse = cbind(snoring_bmi_SAP_joint$gx_se, snoring_bmi_SAP_joint$ga_se),
                           by = snoring_bmi_SAP_joint$gy,
                           byse = snoring_bmi_SAP_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(snoring_bmi_SAP_joint$gx, snoring_bmi_SAP_joint$ga), 
                                   bxse = cbind(snoring_bmi_SAP_joint$gx_se, snoring_bmi_SAP_joint$ga_se),
                                   by = snoring_bmi_SAP_joint$gy,
                                   byse = snoring_bmi_SAP_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(snoring_bmi_SAP_joint$gx, snoring_bmi_SAP_joint$ga), 
                               bxse = cbind(snoring_bmi_SAP_joint$gx_se, snoring_bmi_SAP_joint$ga_se),
                               by = snoring_bmi_SAP_joint$gy,
                               byse = snoring_bmi_SAP_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(snoring_bmi_SAP_joint$gx, snoring_bmi_SAP_joint$ga), 
                               bxse = cbind(snoring_bmi_SAP_joint$gx_se, snoring_bmi_SAP_joint$ga_se),
                               by = snoring_bmi_SAP_joint$gy,
                               byse = snoring_bmi_SAP_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Snoring", "BMI", "Snoring", "BMI")
table$Outcome <- c("SAP","SAP","SAP","SAP")
table <- rbind(table0, table)
table0=table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-snoring-SAP")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test <- rbind(test0, test)
test0 = test






##UAP
#extract SNP infomation from GWAS
snoring_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                   paste(joint_snoring$SNP))]
snoring_bmi_UAP_joint_gx <- ckb_snoring_summary1[ which( paste(ckb_snoring_summary1$SNP) %in% 
                                                           paste(snoring_bmi_joint_ga$SNP)),]
snoring_bmi_UAP_joint_gx <- snoring_bmi_UAP_joint_gx[!duplicated(snoring_bmi_UAP_joint_gx$SNP),]
snoring_bmi_UAP_joint_gy <- bbj_UAP_summary1[ which( paste(bbj_UAP_summary1$SNP) %in% 
                                                       paste(snoring_bmi_UAP_joint_gx$SNP)),]
snoring_bmi_UAP_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                     %in% paste(snoring_bmi_UAP_joint_gx$SNP)),]
#arrange
snoring_bmi_UAP_joint_gx <- arrange(snoring_bmi_UAP_joint_gx, SNP)
snoring_bmi_UAP_joint_gy <- arrange(snoring_bmi_UAP_joint_gy, SNP)
snoring_bmi_UAP_joint_ga <- arrange(snoring_bmi_UAP_joint_ga, SNP)

#????Á´beta
snoring_bmi_UAP_joint_gy$gy <- ifelse(snoring_bmi_UAP_joint_gy$ea == snoring_bmi_UAP_joint_gx$ea, 
                                      snoring_bmi_UAP_joint_gy$beta,
                                      -1*snoring_bmi_UAP_joint_gy$beta)
snoring_bmi_UAP_joint_ga$ga <- ifelse(snoring_bmi_UAP_joint_ga$ea == snoring_bmi_UAP_joint_gx$ea,
                                      snoring_bmi_UAP_joint_ga$beta,
                                      -1*snoring_bmi_UAP_joint_ga$beta)

#MVMR???Ý¿?
snoring_bmi_UAP_joint <- data.frame("V1"=1:57)
snoring_bmi_UAP_joint_gx$SNP <- snoring_bmi_UAP_joint_gx$SNP
snoring_bmi_UAP_joint$ea <- snoring_bmi_UAP_joint_gx$ea
snoring_bmi_UAP_joint$gx <- snoring_bmi_UAP_joint_gx$beta
snoring_bmi_UAP_joint$gx_se <- snoring_bmi_UAP_joint_gx$se
snoring_bmi_UAP_joint$gy <- snoring_bmi_UAP_joint_gy$gy
snoring_bmi_UAP_joint$gy_se <- snoring_bmi_UAP_joint_gy$se
snoring_bmi_UAP_joint$ga <- snoring_bmi_UAP_joint_ga$ga
snoring_bmi_UAP_joint$ga_se <- snoring_bmi_UAP_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(snoring_bmi_UAP_joint$gx, snoring_bmi_UAP_joint$ga), 
                           bxse = cbind(snoring_bmi_UAP_joint$gx_se, snoring_bmi_UAP_joint$ga_se),
                           by = snoring_bmi_UAP_joint$gy,
                           byse = snoring_bmi_UAP_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(snoring_bmi_UAP_joint$gx, snoring_bmi_UAP_joint$ga), 
                                   bxse = cbind(snoring_bmi_UAP_joint$gx_se, snoring_bmi_UAP_joint$ga_se),
                                   by = snoring_bmi_UAP_joint$gy,
                                   byse = snoring_bmi_UAP_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(snoring_bmi_UAP_joint$gx, snoring_bmi_UAP_joint$ga), 
                               bxse = cbind(snoring_bmi_UAP_joint$gx_se, snoring_bmi_UAP_joint$ga_se),
                               by = snoring_bmi_UAP_joint$gy,
                               byse = snoring_bmi_UAP_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(snoring_bmi_UAP_joint$gx, snoring_bmi_UAP_joint$ga), 
                               bxse = cbind(snoring_bmi_UAP_joint$gx_se, snoring_bmi_UAP_joint$ga_se),
                               by = snoring_bmi_UAP_joint$gy,
                               byse = snoring_bmi_UAP_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Snoring", "BMI", "Snoring", "BMI")
table$Outcome <- c("UAP","UAP","UAP","UAP")
table <- rbind(table0, table)
table0=table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-snoring-UAP")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test <- rbind(test0, test)
test0 = test






#CAD---??Î»??-??gyÎ»????Îª×¼
# snps59 <- ckb_snoring_summary1$SNP[which(paste(ckb_snoring_summary1$SNP) %in% paste(snoring_bmi_T2D_joint_gx$SNP))]
# snps55 <- bbj_CAD_summary1$SNP[which(paste(bbj_CAD_summary1$SNP) %in% paste(snoring_bmi_CAD_joint_gx$SNP))]
# snps59[which(!snps59 %in% snps55)]
#"rs12885071" "rs11642015" "rs6504568"  "rs35560038"
#?Ãµ????Õµ?joint-snps
joint_snoring2 <- joint_snoring[which(joint_snoring$SNP!="rs12885071" & 
                                        joint_snoring$SNP!="rs11642015" &
                                        joint_snoring$SNP!="rs6504568" &
                                        joint_snoring$SNP!="rs35560038"),]


#extract SNP infomation from GWAS
snoring_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                   paste(joint_snoring2$SNP))]
snoring_bmi_CAD_joint_gx <- ckb_snoring_summary1[ which( paste(ckb_snoring_summary1$SNP) %in% 
                                                           paste(snoring_bmi_joint_ga$SNP)),]
snoring_bmi_CAD_joint_gx <- snoring_bmi_CAD_joint_gx[!duplicated(snoring_bmi_CAD_joint_gx$SNP),]
snoring_bmi_CAD_joint_gy <- bbj_CAD_summary1[ which( paste(bbj_CAD_summary1$SNP) %in% 
                                                       paste(snoring_bmi_CAD_joint_gx$SNP)),]
snoring_bmi_CAD_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                     %in% paste(snoring_bmi_CAD_joint_gx$SNP)),]

#arrange
snoring_bmi_CAD_joint_gx <- arrange(snoring_bmi_CAD_joint_gx, SNP)
snoring_bmi_CAD_joint_gy <- arrange(snoring_bmi_CAD_joint_gy, SNP)
snoring_bmi_CAD_joint_ga <- arrange(snoring_bmi_CAD_joint_ga, SNP)

#????Á´beta
snoring_bmi_CAD_joint_gy$gy <- ifelse(snoring_bmi_CAD_joint_gy$ea == snoring_bmi_CAD_joint_gx$ea, 
                                      snoring_bmi_CAD_joint_gy$beta,
                                      -1*snoring_bmi_CAD_joint_gy$beta)
snoring_bmi_CAD_joint_ga$ga <- ifelse(snoring_bmi_CAD_joint_ga$ea == snoring_bmi_CAD_joint_gx$ea,
                                      snoring_bmi_CAD_joint_ga$beta,
                                      -1*snoring_bmi_CAD_joint_ga$beta)

#MVMR???Ý¿?-×¢???Þ¸?SNP??
snoring_bmi_CAD_joint <- data.frame("V1"=1:54)
snoring_bmi_CAD_joint_gx$SNP <- snoring_bmi_CAD_joint_gx$SNP
snoring_bmi_CAD_joint$ea <- snoring_bmi_CAD_joint_gx$ea
snoring_bmi_CAD_joint$gx <- snoring_bmi_CAD_joint_gx$beta
snoring_bmi_CAD_joint$gx_se <- snoring_bmi_CAD_joint_gx$se
snoring_bmi_CAD_joint$gy <- snoring_bmi_CAD_joint_gy$gy
snoring_bmi_CAD_joint$gy_se <- snoring_bmi_CAD_joint_gy$se
snoring_bmi_CAD_joint$ga <- snoring_bmi_CAD_joint_ga$ga
snoring_bmi_CAD_joint$ga_se <- snoring_bmi_CAD_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(snoring_bmi_CAD_joint$gx, snoring_bmi_CAD_joint$ga), 
                           bxse = cbind(snoring_bmi_CAD_joint$gx_se, snoring_bmi_CAD_joint$ga_se),
                           by = snoring_bmi_CAD_joint$gy,
                           byse = snoring_bmi_CAD_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(snoring_bmi_CAD_joint$gx, snoring_bmi_CAD_joint$ga), 
                                   bxse = cbind(snoring_bmi_CAD_joint$gx_se, snoring_bmi_CAD_joint$ga_se),
                                   by = snoring_bmi_CAD_joint$gy,
                                   byse = snoring_bmi_CAD_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(snoring_bmi_CAD_joint$gx, snoring_bmi_CAD_joint$ga), 
                               bxse = cbind(snoring_bmi_CAD_joint$gx_se, snoring_bmi_CAD_joint$ga_se),
                               by = snoring_bmi_CAD_joint$gy,
                               byse = snoring_bmi_CAD_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(snoring_bmi_CAD_joint$gx, snoring_bmi_CAD_joint$ga), 
                               bxse = cbind(snoring_bmi_CAD_joint$gx_se, snoring_bmi_CAD_joint$ga_se),
                               by = snoring_bmi_CAD_joint$gy,
                               byse = snoring_bmi_CAD_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Snoring", "BMI", "Snoring", "BMI")
table$Outcome <- c("CAD","CAD","CAD","CAD")
table <- rbind(table0, table)
table0=table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-snoring-CAD")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test <- rbind(test0, test)
test0 = test








#############################
#ALL IVs-habitual
#############################

#joint snps from each exposure


joint_habitual <- as.data.frame(c(ckb_habitual_iv$SNP, ckb_bmi_iv$SNP))
colnames(joint_habitual)[1]<-'SNP'
joint_habitual <- clump_data(joint_habitual, clump_r2 = 0.001, 
                             clump_kb = 10000, pop = "EAS")

proxy_bmi <- fread("C:/Users/86132/Desktop/????/snoring2CAD_check/data/bmi_proxy.txt",h=T)


#?Ãµ????Õµ?joint-snps
joint_habitual1 <- joint_habitual[which(joint_habitual$SNP!="rs1532127" & 
                                          joint_habitual$SNP!="rs6044722" &
                                          joint_habitual$SNP!="rs4973506" &
                                          joint_habitual$SNP!="rs2744475"),]
joint_habitual2 <- rbind(joint_habitual1, proxy_bmi)
joint_habitual = joint_habitual2









#MI
#extract SNP infomation from GWAS
habitual_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                    paste(joint_habitual$SNP))]
habitual_bmi_MI_joint_gx <- ckb_habitual_summary1[ which( paste(ckb_habitual_summary1$SNP) %in% 
                                                            paste(habitual_bmi_joint_ga$SNP)),]
habitual_bmi_MI_joint_gx <- habitual_bmi_MI_joint_gx[!duplicated(habitual_bmi_MI_joint_gx$SNP),]
habitual_bmi_MI_joint_gy <- bbj_MI_summary1[ which( paste(bbj_MI_summary1$SNP) %in% 
                                                      paste(habitual_bmi_MI_joint_gx$SNP))]
habitual_bmi_MI_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                     %in% paste(habitual_bmi_MI_joint_gx$SNP)),]
#arrange
habitual_bmi_MI_joint_gx <- arrange(habitual_bmi_MI_joint_gx, SNP)
habitual_bmi_MI_joint_gy <- arrange(habitual_bmi_MI_joint_gy, SNP)
habitual_bmi_MI_joint_ga <- arrange(habitual_bmi_MI_joint_ga, SNP)

#????Á´beta
habitual_bmi_MI_joint_gy$gy <- ifelse(habitual_bmi_MI_joint_gy$ea == habitual_bmi_MI_joint_gx$ea, 
                                      habitual_bmi_MI_joint_gy$beta,
                                      -1*habitual_bmi_MI_joint_gy$beta)
habitual_bmi_MI_joint_ga$ga <- ifelse(habitual_bmi_MI_joint_ga$ea == habitual_bmi_MI_joint_gx$ea,
                                      habitual_bmi_MI_joint_ga$beta,
                                      -1*habitual_bmi_MI_joint_ga$beta)

#MVMR???Ý¿?
habitual_bmi_MI_joint <- data.frame("V1"=1:57)
habitual_bmi_MI_joint_gx$SNP <- habitual_bmi_MI_joint_gx$SNP
habitual_bmi_MI_joint$ea <- habitual_bmi_MI_joint_gx$ea
habitual_bmi_MI_joint$gx <- habitual_bmi_MI_joint_gx$beta
habitual_bmi_MI_joint$gx_se <- habitual_bmi_MI_joint_gx$se
habitual_bmi_MI_joint$gy <- habitual_bmi_MI_joint_gy$gy
habitual_bmi_MI_joint$gy_se <- habitual_bmi_MI_joint_gy$se
habitual_bmi_MI_joint$ga <- habitual_bmi_MI_joint_ga$ga
habitual_bmi_MI_joint$ga_se <- habitual_bmi_MI_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(habitual_bmi_MI_joint$gx, habitual_bmi_MI_joint$ga), 
                           bxse = cbind(habitual_bmi_MI_joint$gx_se, habitual_bmi_MI_joint$ga_se),
                           by = habitual_bmi_MI_joint$gy,
                           byse = habitual_bmi_MI_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(habitual_bmi_MI_joint$gx, habitual_bmi_MI_joint$ga), 
                                   bxse = cbind(habitual_bmi_MI_joint$gx_se, habitual_bmi_MI_joint$ga_se),
                                   by = habitual_bmi_MI_joint$gy,
                                   byse = habitual_bmi_MI_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(habitual_bmi_MI_joint$gx, habitual_bmi_MI_joint$ga), 
                               bxse = cbind(habitual_bmi_MI_joint$gx_se, habitual_bmi_MI_joint$ga_se),
                               by = habitual_bmi_MI_joint$gy,
                               byse = habitual_bmi_MI_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(habitual_bmi_MI_joint$gx, habitual_bmi_MI_joint$ga), 
                               bxse = cbind(habitual_bmi_MI_joint$gx_se, habitual_bmi_MI_joint$ga_se),
                               by = habitual_bmi_MI_joint$gy,
                               byse = habitual_bmi_MI_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Habitual", "BMI", "Habitual", "BMI")
table$Outcome <- c("MI","MI","MI","MI")
table <- rbind(table0, table)
table0=table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-habitual-MI")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test <- rbind(test0, test)
test0 = test







#CHF
#extract SNP infomation from GWAS
habitual_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                    paste(joint_habitual$SNP))]
habitual_bmi_CHF_joint_gx <- ckb_habitual_summary1[ which( paste(ckb_habitual_summary1$SNP) %in% 
                                                             paste(habitual_bmi_joint_ga$SNP)),]
habitual_bmi_CHF_joint_gx <- habitual_bmi_CHF_joint_gx[!duplicated(habitual_bmi_CHF_joint_gx$SNP),]
habitual_bmi_CHF_joint_gy <- bbj_CHF_summary1[ which( paste(bbj_CHF_summary1$SNP) %in% 
                                                        paste(habitual_bmi_CHF_joint_gx$SNP))]
habitual_bmi_CHF_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                      %in% paste(habitual_bmi_CHF_joint_gx$SNP)),]
#arrange
habitual_bmi_CHF_joint_gx <- arrange(habitual_bmi_CHF_joint_gx, SNP)
habitual_bmi_CHF_joint_gy <- arrange(habitual_bmi_CHF_joint_gy, SNP)
habitual_bmi_CHF_joint_ga <- arrange(habitual_bmi_CHF_joint_ga, SNP)

#????Á´beta
habitual_bmi_CHF_joint_gy$gy <- ifelse(habitual_bmi_CHF_joint_gy$ea == habitual_bmi_CHF_joint_gx$ea, 
                                       habitual_bmi_CHF_joint_gy$beta,
                                       -1*habitual_bmi_CHF_joint_gy$beta)
habitual_bmi_CHF_joint_ga$ga <- ifelse(habitual_bmi_CHF_joint_ga$ea == habitual_bmi_CHF_joint_gx$ea,
                                       habitual_bmi_CHF_joint_ga$beta,
                                       -1*habitual_bmi_CHF_joint_ga$beta)

#MVMR???Ý¿?
habitual_bmi_CHF_joint <- data.frame("V1"=1:57)
habitual_bmi_CHF_joint_gx$SNP <- habitual_bmi_CHF_joint_gx$SNP
habitual_bmi_CHF_joint$ea <- habitual_bmi_CHF_joint_gx$ea
habitual_bmi_CHF_joint$gx <- habitual_bmi_CHF_joint_gx$beta
habitual_bmi_CHF_joint$gx_se <- habitual_bmi_CHF_joint_gx$se
habitual_bmi_CHF_joint$gy <- habitual_bmi_CHF_joint_gy$gy
habitual_bmi_CHF_joint$gy_se <- habitual_bmi_CHF_joint_gy$se
habitual_bmi_CHF_joint$ga <- habitual_bmi_CHF_joint_ga$ga
habitual_bmi_CHF_joint$ga_se <- habitual_bmi_CHF_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(habitual_bmi_CHF_joint$gx, habitual_bmi_CHF_joint$ga), 
                           bxse = cbind(habitual_bmi_CHF_joint$gx_se, habitual_bmi_CHF_joint$ga_se),
                           by = habitual_bmi_CHF_joint$gy,
                           byse = habitual_bmi_CHF_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(habitual_bmi_CHF_joint$gx, habitual_bmi_CHF_joint$ga), 
                                   bxse = cbind(habitual_bmi_CHF_joint$gx_se, habitual_bmi_CHF_joint$ga_se),
                                   by = habitual_bmi_CHF_joint$gy,
                                   byse = habitual_bmi_CHF_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(habitual_bmi_CHF_joint$gx, habitual_bmi_CHF_joint$ga), 
                               bxse = cbind(habitual_bmi_CHF_joint$gx_se, habitual_bmi_CHF_joint$ga_se),
                               by = habitual_bmi_CHF_joint$gy,
                               byse = habitual_bmi_CHF_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(habitual_bmi_CHF_joint$gx, habitual_bmi_CHF_joint$ga), 
                               bxse = cbind(habitual_bmi_CHF_joint$gx_se, habitual_bmi_CHF_joint$ga_se),
                               by = habitual_bmi_CHF_joint$gy,
                               byse = habitual_bmi_CHF_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Habitual", "BMI", "Habitual", "BMI")
table$Outcome <- c("CHF","CHF","CHF","CHF")
table <- rbind(table0, table)
table0=table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-habitual-CHF")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test <- rbind(test0, test)
test0 = test









#Angina
#extract SNP infomation from GWAS
habitual_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                    paste(joint_habitual$SNP))]
habitual_bmi_Angina_joint_gx <- ckb_habitual_summary1[ which( paste(ckb_habitual_summary1$SNP) %in% 
                                                                paste(habitual_bmi_joint_ga$SNP)),]
habitual_bmi_Angina_joint_gx <- habitual_bmi_Angina_joint_gx[!duplicated(habitual_bmi_Angina_joint_gx$SNP),]
habitual_bmi_Angina_joint_gy <- bbj_Angina_summary1[ which( paste(bbj_Angina_summary1$SNP) %in% 
                                                              paste(habitual_bmi_Angina_joint_gx$SNP))]
habitual_bmi_Angina_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                         %in% paste(habitual_bmi_Angina_joint_gx$SNP)),]
#arrange
habitual_bmi_Angina_joint_gx <- arrange(habitual_bmi_Angina_joint_gx, SNP)
habitual_bmi_Angina_joint_gy <- arrange(habitual_bmi_Angina_joint_gy, SNP)
habitual_bmi_Angina_joint_ga <- arrange(habitual_bmi_Angina_joint_ga, SNP)

#????Á´beta
habitual_bmi_Angina_joint_gy$gy <- ifelse(habitual_bmi_Angina_joint_gy$ea == habitual_bmi_Angina_joint_gx$ea, 
                                          habitual_bmi_Angina_joint_gy$beta,
                                          -1*habitual_bmi_Angina_joint_gy$beta)
habitual_bmi_Angina_joint_ga$ga <- ifelse(habitual_bmi_Angina_joint_ga$ea == habitual_bmi_Angina_joint_gx$ea,
                                          habitual_bmi_Angina_joint_ga$beta,
                                          -1*habitual_bmi_Angina_joint_ga$beta)

#MVMR???Ý¿?
habitual_bmi_Angina_joint <- data.frame("V1"=1:57)
habitual_bmi_Angina_joint_gx$SNP <- habitual_bmi_Angina_joint_gx$SNP
habitual_bmi_Angina_joint$ea <- habitual_bmi_Angina_joint_gx$ea
habitual_bmi_Angina_joint$gx <- habitual_bmi_Angina_joint_gx$beta
habitual_bmi_Angina_joint$gx_se <- habitual_bmi_Angina_joint_gx$se
habitual_bmi_Angina_joint$gy <- habitual_bmi_Angina_joint_gy$gy
habitual_bmi_Angina_joint$gy_se <- habitual_bmi_Angina_joint_gy$se
habitual_bmi_Angina_joint$ga <- habitual_bmi_Angina_joint_ga$ga
habitual_bmi_Angina_joint$ga_se <- habitual_bmi_Angina_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(habitual_bmi_Angina_joint$gx, habitual_bmi_Angina_joint$ga), 
                           bxse = cbind(habitual_bmi_Angina_joint$gx_se, habitual_bmi_Angina_joint$ga_se),
                           by = habitual_bmi_Angina_joint$gy,
                           byse = habitual_bmi_Angina_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(habitual_bmi_Angina_joint$gx, habitual_bmi_Angina_joint$ga), 
                                   bxse = cbind(habitual_bmi_Angina_joint$gx_se, habitual_bmi_Angina_joint$ga_se),
                                   by = habitual_bmi_Angina_joint$gy,
                                   byse = habitual_bmi_Angina_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(habitual_bmi_Angina_joint$gx, habitual_bmi_Angina_joint$ga), 
                               bxse = cbind(habitual_bmi_Angina_joint$gx_se, habitual_bmi_Angina_joint$ga_se),
                               by = habitual_bmi_Angina_joint$gy,
                               byse = habitual_bmi_Angina_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(habitual_bmi_Angina_joint$gx, habitual_bmi_Angina_joint$ga), 
                               bxse = cbind(habitual_bmi_Angina_joint$gx_se, habitual_bmi_Angina_joint$ga_se),
                               by = habitual_bmi_Angina_joint$gy,
                               byse = habitual_bmi_Angina_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Habitual", "BMI", "Habitual", "BMI")
table$Outcome <- c("Angina","Angina","Angina","Angina")
table <- rbind(table0, table)
table0=table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-habitual-Angina")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test <- rbind(test0, test)
test0 = test







#SAP
#extract SNP infomation from GWAS
habitual_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                    paste(joint_habitual$SNP))]
habitual_bmi_SAP_joint_gx <- ckb_habitual_summary1[ which( paste(ckb_habitual_summary1$SNP) %in% 
                                                             paste(habitual_bmi_joint_ga$SNP)),]
habitual_bmi_SAP_joint_gx <- habitual_bmi_SAP_joint_gx[!duplicated(habitual_bmi_SAP_joint_gx$SNP),]
habitual_bmi_SAP_joint_gy <- bbj_SAP_summary1[ which( paste(bbj_SAP_summary1$SNP) %in% 
                                                        paste(habitual_bmi_SAP_joint_gx$SNP))]
habitual_bmi_SAP_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                      %in% paste(habitual_bmi_SAP_joint_gx$SNP)),]
#arrange
habitual_bmi_SAP_joint_gx <- arrange(habitual_bmi_SAP_joint_gx, SNP)
habitual_bmi_SAP_joint_gy <- arrange(habitual_bmi_SAP_joint_gy, SNP)
habitual_bmi_SAP_joint_ga <- arrange(habitual_bmi_SAP_joint_ga, SNP)

#????Á´beta
habitual_bmi_SAP_joint_gy$gy <- ifelse(habitual_bmi_SAP_joint_gy$ea == habitual_bmi_SAP_joint_gx$ea, 
                                       habitual_bmi_SAP_joint_gy$beta,
                                       -1*habitual_bmi_SAP_joint_gy$beta)
habitual_bmi_SAP_joint_ga$ga <- ifelse(habitual_bmi_SAP_joint_ga$ea == habitual_bmi_SAP_joint_gx$ea,
                                       habitual_bmi_SAP_joint_ga$beta,
                                       -1*habitual_bmi_SAP_joint_ga$beta)

#MVMR???Ý¿?
habitual_bmi_SAP_joint <- data.frame("V1"=1:57)
habitual_bmi_SAP_joint_gx$SNP <- habitual_bmi_SAP_joint_gx$SNP
habitual_bmi_SAP_joint$ea <- habitual_bmi_SAP_joint_gx$ea
habitual_bmi_SAP_joint$gx <- habitual_bmi_SAP_joint_gx$beta
habitual_bmi_SAP_joint$gx_se <- habitual_bmi_SAP_joint_gx$se
habitual_bmi_SAP_joint$gy <- habitual_bmi_SAP_joint_gy$gy
habitual_bmi_SAP_joint$gy_se <- habitual_bmi_SAP_joint_gy$se
habitual_bmi_SAP_joint$ga <- habitual_bmi_SAP_joint_ga$ga
habitual_bmi_SAP_joint$ga_se <- habitual_bmi_SAP_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(habitual_bmi_SAP_joint$gx, habitual_bmi_SAP_joint$ga), 
                           bxse = cbind(habitual_bmi_SAP_joint$gx_se, habitual_bmi_SAP_joint$ga_se),
                           by = habitual_bmi_SAP_joint$gy,
                           byse = habitual_bmi_SAP_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(habitual_bmi_SAP_joint$gx, habitual_bmi_SAP_joint$ga), 
                                   bxse = cbind(habitual_bmi_SAP_joint$gx_se, habitual_bmi_SAP_joint$ga_se),
                                   by = habitual_bmi_SAP_joint$gy,
                                   byse = habitual_bmi_SAP_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(habitual_bmi_SAP_joint$gx, habitual_bmi_SAP_joint$ga), 
                               bxse = cbind(habitual_bmi_SAP_joint$gx_se, habitual_bmi_SAP_joint$ga_se),
                               by = habitual_bmi_SAP_joint$gy,
                               byse = habitual_bmi_SAP_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(habitual_bmi_SAP_joint$gx, habitual_bmi_SAP_joint$ga), 
                               bxse = cbind(habitual_bmi_SAP_joint$gx_se, habitual_bmi_SAP_joint$ga_se),
                               by = habitual_bmi_SAP_joint$gy,
                               byse = habitual_bmi_SAP_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Habitual", "BMI", "Habitual", "BMI")
table$Outcome <- c("SAP","SAP","SAP","SAP")
table <- rbind(table0, table)
table0=table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-habitual-SAP")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test <- rbind(test0, test)
test0 = test





#UAP
#extract SNP infomation from GWAS
habitual_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                    paste(joint_habitual$SNP))]
habitual_bmi_UAP_joint_gx <- ckb_habitual_summary1[ which( paste(ckb_habitual_summary1$SNP) %in% 
                                                             paste(habitual_bmi_joint_ga$SNP)),]
habitual_bmi_UAP_joint_gx <- habitual_bmi_UAP_joint_gx[!duplicated(habitual_bmi_UAP_joint_gx$SNP),]
habitual_bmi_UAP_joint_gy <- bbj_UAP_summary1[ which( paste(bbj_UAP_summary1$SNP) %in% 
                                                        paste(habitual_bmi_UAP_joint_gx$SNP))]
habitual_bmi_UAP_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                      %in% paste(habitual_bmi_UAP_joint_gx$SNP)),]
#arrange
habitual_bmi_UAP_joint_gx <- arrange(habitual_bmi_UAP_joint_gx, SNP)
habitual_bmi_UAP_joint_gy <- arrange(habitual_bmi_UAP_joint_gy, SNP)
habitual_bmi_UAP_joint_ga <- arrange(habitual_bmi_UAP_joint_ga, SNP)

#????Á´beta
habitual_bmi_UAP_joint_gy$gy <- ifelse(habitual_bmi_UAP_joint_gy$ea == habitual_bmi_UAP_joint_gx$ea, 
                                       habitual_bmi_UAP_joint_gy$beta,
                                       -1*habitual_bmi_UAP_joint_gy$beta)
habitual_bmi_UAP_joint_ga$ga <- ifelse(habitual_bmi_UAP_joint_ga$ea == habitual_bmi_UAP_joint_gx$ea,
                                       habitual_bmi_UAP_joint_ga$beta,
                                       -1*habitual_bmi_UAP_joint_ga$beta)

#MVMR???Ý¿?
habitual_bmi_UAP_joint <- data.frame("V1"=1:57)
habitual_bmi_UAP_joint_gx$SNP <- habitual_bmi_UAP_joint_gx$SNP
habitual_bmi_UAP_joint$ea <- habitual_bmi_UAP_joint_gx$ea
habitual_bmi_UAP_joint$gx <- habitual_bmi_UAP_joint_gx$beta
habitual_bmi_UAP_joint$gx_se <- habitual_bmi_UAP_joint_gx$se
habitual_bmi_UAP_joint$gy <- habitual_bmi_UAP_joint_gy$gy
habitual_bmi_UAP_joint$gy_se <- habitual_bmi_UAP_joint_gy$se
habitual_bmi_UAP_joint$ga <- habitual_bmi_UAP_joint_ga$ga
habitual_bmi_UAP_joint$ga_se <- habitual_bmi_UAP_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(habitual_bmi_UAP_joint$gx, habitual_bmi_UAP_joint$ga), 
                           bxse = cbind(habitual_bmi_UAP_joint$gx_se, habitual_bmi_UAP_joint$ga_se),
                           by = habitual_bmi_UAP_joint$gy,
                           byse = habitual_bmi_UAP_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(habitual_bmi_UAP_joint$gx, habitual_bmi_UAP_joint$ga), 
                                   bxse = cbind(habitual_bmi_UAP_joint$gx_se, habitual_bmi_UAP_joint$ga_se),
                                   by = habitual_bmi_UAP_joint$gy,
                                   byse = habitual_bmi_UAP_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(habitual_bmi_UAP_joint$gx, habitual_bmi_UAP_joint$ga), 
                               bxse = cbind(habitual_bmi_UAP_joint$gx_se, habitual_bmi_UAP_joint$ga_se),
                               by = habitual_bmi_UAP_joint$gy,
                               byse = habitual_bmi_UAP_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(habitual_bmi_UAP_joint$gx, habitual_bmi_UAP_joint$ga), 
                               bxse = cbind(habitual_bmi_UAP_joint$gx_se, habitual_bmi_UAP_joint$ga_se),
                               by = habitual_bmi_UAP_joint$gy,
                               byse = habitual_bmi_UAP_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Habitual", "BMI", "Habitual", "BMI")
table$Outcome <- c("UAP","UAP","UAP","UAP")
table <- rbind(table0, table)
table0=table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-habitual-UAP")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test <- rbind(test0, test)
test0 = test









#CAD---??Î»??-??gyÎ»????Îª×¼
# snps58 <- ckb_habitual_summary1$SNP[which(paste(ckb_habitual_summary1$SNP) %in% paste(habitual_bmi_T2D_joint_gx$SNP))]
# snps55 <- bbj_CAD_summary1$SNP[which(paste(bbj_CAD_summary1$SNP) %in% paste(habitual_bmi_CAD_joint_gx$SNP))]
# snps58[which(!snps58 %in% snps55)]
#"rs12885071" "rs11642015" "rs6504568"  "rs35560038"
#?Ãµ????Õµ?joint-snps
joint_habitual2 <- joint_habitual[which(joint_habitual$SNP!="rs12885071" & 
                                          joint_habitual$SNP!="rs6504568" &
                                          joint_habitual$SNP!="rs35560038"),]


#extract SNP infomation from GWAS
habitual_bmi_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP) %in% 
                                                    paste(joint_habitual2$SNP))]
habitual_bmi_CAD_joint_gx <- ckb_habitual_summary1[ which( paste(ckb_habitual_summary1$SNP) %in% 
                                                             paste(habitual_bmi_joint_ga$SNP)),]
habitual_bmi_CAD_joint_gx <- habitual_bmi_CAD_joint_gx[!duplicated(habitual_bmi_CAD_joint_gx$SNP),]
habitual_bmi_CAD_joint_gy <- bbj_CAD_summary1[ which( paste(bbj_CAD_summary1$SNP) %in% 
                                                        paste(habitual_bmi_CAD_joint_gx$SNP)),]
habitual_bmi_CAD_joint_ga <- ckb_bmi_summary1[ which( paste(ckb_bmi_summary1$SNP)
                                                      %in% paste(habitual_bmi_CAD_joint_gx$SNP)),]

#arrange
habitual_bmi_CAD_joint_gx <- arrange(habitual_bmi_CAD_joint_gx, SNP)
habitual_bmi_CAD_joint_gy <- arrange(habitual_bmi_CAD_joint_gy, SNP)
habitual_bmi_CAD_joint_ga <- arrange(habitual_bmi_CAD_joint_ga, SNP)

#????Á´beta
habitual_bmi_CAD_joint_gy$gy <- ifelse(habitual_bmi_CAD_joint_gy$ea == habitual_bmi_CAD_joint_gx$ea, 
                                       habitual_bmi_CAD_joint_gy$beta,
                                       -1*habitual_bmi_CAD_joint_gy$beta)
habitual_bmi_CAD_joint_ga$ga <- ifelse(habitual_bmi_CAD_joint_ga$ea == habitual_bmi_CAD_joint_gx$ea,
                                       habitual_bmi_CAD_joint_ga$beta,
                                       -1*habitual_bmi_CAD_joint_ga$beta)

#MVMR???Ý¿?-×¢???Þ¸?SNP??
habitual_bmi_CAD_joint <- data.frame("V1"=1:54)
habitual_bmi_CAD_joint_gx$SNP <- habitual_bmi_CAD_joint_gx$SNP
habitual_bmi_CAD_joint$ea <- habitual_bmi_CAD_joint_gx$ea
habitual_bmi_CAD_joint$gx <- habitual_bmi_CAD_joint_gx$beta
habitual_bmi_CAD_joint$gx_se <- habitual_bmi_CAD_joint_gx$se
habitual_bmi_CAD_joint$gy <- habitual_bmi_CAD_joint_gy$gy
habitual_bmi_CAD_joint$gy_se <- habitual_bmi_CAD_joint_gy$se
habitual_bmi_CAD_joint$ga <- habitual_bmi_CAD_joint_ga$ga
habitual_bmi_CAD_joint$ga_se <- habitual_bmi_CAD_joint_ga$se

ivw <- mr_mvivw(mr_mvinput(bx =cbind(habitual_bmi_CAD_joint$gx, habitual_bmi_CAD_joint$ga), 
                           bxse = cbind(habitual_bmi_CAD_joint$gx_se, habitual_bmi_CAD_joint$ga_se),
                           by = habitual_bmi_CAD_joint$gy,
                           byse = habitual_bmi_CAD_joint$gy_se) )

mvmedian <- mr_mvmedian(mr_mvinput(bx =cbind(habitual_bmi_CAD_joint$gx, habitual_bmi_CAD_joint$ga), 
                                   bxse = cbind(habitual_bmi_CAD_joint$gx_se, habitual_bmi_CAD_joint$ga_se),
                                   by = habitual_bmi_CAD_joint$gy,
                                   byse = habitual_bmi_CAD_joint$gy_se), iterations = 100)

egger <- mr_mvegger(mr_mvinput(bx =cbind(habitual_bmi_CAD_joint$gx, habitual_bmi_CAD_joint$ga), 
                               bxse = cbind(habitual_bmi_CAD_joint$gx_se, habitual_bmi_CAD_joint$ga_se),
                               by = habitual_bmi_CAD_joint$gy,
                               byse = habitual_bmi_CAD_joint$gy_se), orientate = 1)

lasso <- mr_mvlasso(mr_mvinput(bx =cbind(habitual_bmi_CAD_joint$gx, habitual_bmi_CAD_joint$ga), 
                               bxse = cbind(habitual_bmi_CAD_joint$gx_se, habitual_bmi_CAD_joint$ga_se),
                               by = habitual_bmi_CAD_joint$gy,
                               byse = habitual_bmi_CAD_joint$gy_se))

table.ivw <- as.data.frame(cbind(ivw@Exposure,ivw@Outcome, ivw@Estimate, ivw@StdError, 
                                 ivw@CILower, ivw@CIUpper, ivw@Pvalue))
table.egger <- as.data.frame(cbind(egger@Exposure, egger@Outcome, egger@Estimate, 
                                   egger@StdError.Est, egger@CILower.Est, 
                                   egger@CIUpper.Est, egger@Pvalue.Est))
table <- rbind(table.ivw, table.egger)
rownames(table) <- c("MVMR-IVW.gx", "MVMR-IVW.ga", "MVMR-Egger.gx", "MVMR-Egger.ga")
colnames(table) <- c("Exposure","Outcome", "beta", "se", "beta.low", "beta.upper", "P")
table$Exposure <- c("Habitual", "BMI", "Habitual", "BMI")
table$Outcome <- c("CAD","CAD","CAD","CAD")
table <- rbind(table0, table)
table0=table

test <- as.data.frame(cbind(ivw@SNPs, ivw@Heter.Stat[1], ivw@Heter.Stat[2], 
                            egger@Intercept, egger@StdError.Int, egger@CILower.Int, 
                            egger@CIUpper.Int, egger@Pvalue.Int))

rownames(test) <- c("MVMR-test-habitual-CAD")
colnames(test) <- c("nSNPs","Heterogeneity test statistic", "Heterogeneity P", 
                    "Egger.intercept", "Egger.intercept.SE", "Egger.intercept.low",  
                    "Egger.intercept.upper", "Egger.intercept.P")
test <- rbind(test0, test)
test0 = test







#????Ïµ??

table0$OR1.5 <- ifelse(
  table0$Exposure=="BMI", exp(as.numeric(table0$beta)),
  exp(as.numeric(table0$beta) * 0.405)
)


table0$OR_L <- ifelse(
  table0$Exposure=="BMI", exp(as.numeric(table0$beta.low)),
  exp(as.numeric(table0$beta.low) * 0.405)
)


table0$OR_U <- ifelse(
  table0$Exposure=="BMI", exp(as.numeric(table0$beta.upper)),
  exp(as.numeric(table0$beta.upper) * 0.405)
)



write.csv(table0, "C:/Users/86132/Desktop/????/snoring2CAD_check/result/MVMR/20221007_snoring_bmi_CAD_MVMR_MRpackage.csv")
write.csv(test0, "C:/Users/86132/Desktop/????/snoring2CAD_check/result/MVMR/20221007_snoring_bmi_CAD_MVMR_MRpackage_test.csv")

