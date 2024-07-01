########
#By Yunqing Zhu, Canqing Yu
#Email: zhuyun_qing@126.com yucanqing@pku.edu.cn
#########

#MVMR-extract????

#????snoring SVMR?Ľ???
snoring_sv <- read.csv("C:/Users/86132/Desktop/????/snoring2CAD_check/2sampleMR/3snps_snoring_CAD.csv")
colnames(snoring_sv) <- c("Outcome", "Method", "N.SNPs", "OR", "OR_L", "OR_U", "P")
snoring_sv$ORCI <- paste0(sprintf("%0.2f",snoring_sv$OR)," (", 
                          sprintf("%0.2f",snoring_sv$OR_L), ",", 
                          sprintf("%0.2f",snoring_sv$OR_U), ")")
snoring_sv$P <- sprintf("%0.3f",snoring_sv$P)
snoring_sv$Method = ifelse(snoring_sv$Method=="Inverse variance weighted", "SVMR IVW", snoring_sv$Method )
snoring_sv$Method = ifelse(snoring_sv$Method=="Weighted median", "SVMR WM", snoring_sv$Method )
snoring_sv$Method = ifelse(snoring_sv$Method=="MR Egger", "SVMR Egger", snoring_sv$Method )
snoring_sv$Method = ifelse(snoring_sv$Method=="RAPS", "SVMR RAPS", snoring_sv$Method )

snoring_sv$Exposure <- ifelse(snoring_sv$Method=="SVMR IVW",  "Snoring", " ")

snoring_sv1 <- snoring_sv[,c(9,1,2,8,7,4,5,6)]

#???ɷֽ??ֵ?????
snoring_CAD_sv1 <- subset(snoring_sv1, Outcome=="CAD")[,-2]
snoring_MI_sv1 <- subset(snoring_sv1, Outcome=="MI")[,-2]
snoring_CHF_sv1 <- subset(snoring_sv1, Outcome=="CHF")[,-2]
snoring_Angina_sv1 <- subset(snoring_sv1, Outcome=="Angina")[,-2]
snoring_UAP_sv1 <- subset(snoring_sv1, Outcome=="UAP")[,-2]
snoring_SAP_sv1 <- subset(snoring_sv1, Outcome=="SAP")[,-2]







#????habitual SVMR?Ľ???
habitual_sv <- read.csv("C:/Users/86132/Desktop/????/snoring2CAD_check/2sampleMR/3snps_habitual_CAD.csv")
colnames(habitual_sv) <- c("Outcome", "Method", "N.SNPs", "OR", "OR_L", "OR_U", "P")
habitual_sv$ORCI <- paste0(sprintf("%0.2f",habitual_sv$OR)," (", 
                           sprintf("%0.2f",habitual_sv$OR_L), ",", 
                           sprintf("%0.2f",habitual_sv$OR_U), ")")
habitual_sv$P <- sprintf("%0.3f",habitual_sv$P)
habitual_sv$Method = ifelse(habitual_sv$Method=="Inverse variance weighted", "SVMR IVW", habitual_sv$Method )
habitual_sv$Method = ifelse(habitual_sv$Method=="Weighted median", "SVMR WM", habitual_sv$Method )
habitual_sv$Method = ifelse(habitual_sv$Method=="MR Egger", "SVMR Egger", habitual_sv$Method )
habitual_sv$Method = ifelse(habitual_sv$Method=="RAPS", "SVMR RAPS", habitual_sv$Method )

habitual_sv$Exposure <- ifelse(habitual_sv$Method=="SVMR IVW",  "Habitual snoring", " ")

habitual_sv1 <- habitual_sv[,c(9,1,2,8,7,4,5,6)]

#???ɷֽ??ֵ?????
habitual_CAD_sv1 <- subset(habitual_sv1, Outcome=="CAD")[,-2]
habitual_MI_sv1 <- subset(habitual_sv1, Outcome=="MI")[,-2]
habitual_CHF_sv1 <- subset(habitual_sv1, Outcome=="CHF")[,-2]
habitual_Angina_sv1 <- subset(habitual_sv1, Outcome=="Angina")[,-2]
habitual_UAP_sv1 <- subset(habitual_sv1, Outcome=="UAP")[,-2]
habitual_SAP_sv1 <- subset(habitual_sv1, Outcome=="SAP")[,-2]






#????BMI SVMR?Ľ???
BMI_sv <- read.csv("C:/Users/86132/Desktop/????/snoring2CAD_check/2sampleMR/57snps_BMI_CAD.csv")
colnames(BMI_sv) <- c("Outcome", "Method", "N.SNPs", "OR", "OR_L", "OR_U", "P")
BMI_sv$ORCI <- paste0(sprintf("%0.2f",BMI_sv$OR)," (", 
                      sprintf("%0.2f",BMI_sv$OR_L), ",", 
                      sprintf("%0.2f",BMI_sv$OR_U), ")")
BMI_sv$P <- signif(BMI_sv$P,3)
BMI_sv$Method = ifelse(BMI_sv$Method=="Inverse variance weighted", "SVMR IVW", BMI_sv$Method )
BMI_sv$Method = ifelse(BMI_sv$Method=="Weighted median", "SVMR WM", BMI_sv$Method )
BMI_sv$Method = ifelse(BMI_sv$Method=="MR Egger", "SVMR Egger", BMI_sv$Method )
BMI_sv$Method = ifelse(BMI_sv$Method=="RAPS", "SVMR RAPS", BMI_sv$Method )

BMI_sv$Exposure <- ifelse(BMI_sv$Method=="SVMR IVW",  "BMI", " ")

BMI_sv1 <- BMI_sv[,c(9,1,2,8,7,4,5,6)]

#???ɷֽ??ֵ?????
BMI_CAD_sv1 <- subset(BMI_sv1, Outcome=="CAD")[,-2]
BMI_MI_sv1 <- subset(BMI_sv1, Outcome=="MI")[,-2]
BMI_CHF_sv1 <- subset(BMI_sv1, Outcome=="CHF")[,-2]
BMI_Angina_sv1 <- subset(BMI_sv1, Outcome=="Angina")[,-2]
BMI_UAP_sv1 <- subset(BMI_sv1, Outcome=="UAP")[,-2]
BMI_SAP_sv1 <- subset(BMI_sv1, Outcome=="SAP")[,-2]




#????MVMR??????
snoring_bmi_mv <- read.csv("C:/Users/86132/Desktop/????/snoring2CAD_check/result/MVMR/20221007_snoring_bmi_CAD_MVMR_MRpackage.csv")
colnames(snoring_bmi_mv) <- c("Method", "Exposure", "Outcome", "beta", "se", 
                              "beta_l", "beta_u", "P", "OR", "OR_L", "OR_U")
snoring_bmi_mv$ORCI <- paste0(sprintf("%0.2f",snoring_bmi_mv$OR)," (", 
                              sprintf("%0.2f",snoring_bmi_mv$OR_L), ",", 
                              sprintf("%0.2f",snoring_bmi_mv$OR_U), ")")

snoring_bmi_mv$P <- sprintf("%0.3f",snoring_bmi_mv$P)
snoring_bmi_mv$Method <- rep(c("MVMR IVW", "MVMR IVW", "MVMR Egger", "MVMR Egger"),nrow(snoring_bmi_mv)/4)
head(snoring_bmi_mv)
snoring_bmi_mv1 <- snoring_bmi_mv[,c(2,1,3,12,8,9,10,11)]

#????MVMR????-snoring
snoring_MI_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Snoring" &
                           snoring_bmi_mv1$Outcome=="MI")[,-3]
snoring_MI_mv1$Exposure=c(" ", " ")
snoring_CHF_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Snoring" &
                            snoring_bmi_mv1$Outcome=="CHF")[,-3]
snoring_CHF_mv1$Exposure=c(" ", " ")
snoring_Angina_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Snoring" &
                               snoring_bmi_mv1$Outcome=="Angina")[,-3]
snoring_Angina_mv1$Exposure=c(" ", " ")
snoring_SAP_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Snoring" &
                            snoring_bmi_mv1$Outcome=="SAP")[,-3]
snoring_SAP_mv1$Exposure=c(" ", " ")
snoring_UAP_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Snoring" &
                            snoring_bmi_mv1$Outcome=="UAP")[,-3]
snoring_UAP_mv1$Exposure=c(" ", " ")
snoring_CAD_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Snoring" &
                            snoring_bmi_mv1$Outcome=="CAD")[,-3]
snoring_CAD_mv1$Exposure=c(" ", " ")


#????MVMR????-BMI in snoring
snoring_insnoring_bmi_mv1 <- snoring_bmi_mv1[c(1:24),]

BMI_insnoring_MI_mv1 <- subset(snoring_insnoring_bmi_mv1, snoring_insnoring_bmi_mv1$Exposure=="BMI" &
                                 snoring_insnoring_bmi_mv1$Outcome=="MI")[,-3]
BMI_insnoring_MI_mv1$Exposure=c(" ", " ")
BMI_insnoring_CHF_mv1 <- subset(snoring_insnoring_bmi_mv1, snoring_insnoring_bmi_mv1$Exposure=="BMI" &
                                  snoring_insnoring_bmi_mv1$Outcome=="CHF")[,-3]
BMI_insnoring_CHF_mv1$Exposure=c(" ", " ")
BMI_insnoring_Angina_mv1 <- subset(snoring_insnoring_bmi_mv1, snoring_insnoring_bmi_mv1$Exposure=="BMI" &
                                     snoring_insnoring_bmi_mv1$Outcome=="Angina")[,-3]
BMI_insnoring_Angina_mv1$Exposure=c(" ", " ")
BMI_insnoring_SAP_mv1 <- subset(snoring_insnoring_bmi_mv1, snoring_insnoring_bmi_mv1$Exposure=="BMI" &
                                  snoring_insnoring_bmi_mv1$Outcome=="SAP")[,-3]
BMI_insnoring_SAP_mv1$Exposure=c(" ", " ")
BMI_insnoring_UAP_mv1 <- subset(snoring_insnoring_bmi_mv1, snoring_insnoring_bmi_mv1$Exposure=="BMI" &
                                  snoring_insnoring_bmi_mv1$Outcome=="UAP")[,-3]
BMI_insnoring_UAP_mv1$Exposure=c(" ", " ")
BMI_insnoring_CAD_mv1 <- subset(snoring_insnoring_bmi_mv1, snoring_insnoring_bmi_mv1$Exposure=="BMI" &
                                  snoring_insnoring_bmi_mv1$Outcome=="CAD")[,-3]
BMI_insnoring_CAD_mv1$Exposure=c(" ", " ")




#????MVMR????-habitual
habitual_MI_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Habitual" &
                            snoring_bmi_mv1$Outcome=="MI")[,-3]
habitual_MI_mv1$Exposure=c(" ", " ")
habitual_CHF_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Habitual" &
                             snoring_bmi_mv1$Outcome=="CHF")[,-3]
habitual_CHF_mv1$Exposure=c(" ", " ")
habitual_Angina_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Habitual" &
                                snoring_bmi_mv1$Outcome=="Angina")[,-3]
habitual_Angina_mv1$Exposure=c(" ", " ")
habitual_SAP_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Habitual" &
                             snoring_bmi_mv1$Outcome=="SAP")[,-3]
habitual_SAP_mv1$Exposure=c(" ", " ")
habitual_UAP_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Habitual" &
                             snoring_bmi_mv1$Outcome=="UAP")[,-3]
habitual_UAP_mv1$Exposure=c(" ", " ")
habitual_CAD_mv1 <- subset(snoring_bmi_mv1, snoring_bmi_mv1$Exposure=="Habitual" &
                             snoring_bmi_mv1$Outcome=="CAD")[,-3]
habitual_CAD_mv1$Exposure=c(" ", " ")


#????MVMR????-BMI in habitual
snoring_inhabitual_bmi_mv1 <- snoring_bmi_mv1[c(25:48),]

BMI_inhabitual_MI_mv1 <- subset(snoring_inhabitual_bmi_mv1, snoring_inhabitual_bmi_mv1$Exposure=="BMI" &
                                  snoring_inhabitual_bmi_mv1$Outcome=="MI")[,-3]
BMI_inhabitual_MI_mv1$Exposure=c(" ", " ")
BMI_inhabitual_CHF_mv1 <- subset(snoring_inhabitual_bmi_mv1, snoring_inhabitual_bmi_mv1$Exposure=="BMI" &
                                   snoring_inhabitual_bmi_mv1$Outcome=="CHF")[,-3]
BMI_inhabitual_CHF_mv1$Exposure=c(" ", " ")
BMI_inhabitual_Angina_mv1 <- subset(snoring_inhabitual_bmi_mv1, snoring_inhabitual_bmi_mv1$Exposure=="BMI" &
                                      snoring_inhabitual_bmi_mv1$Outcome=="Angina")[,-3]
BMI_inhabitual_Angina_mv1$Exposure=c(" ", " ")
BMI_inhabitual_SAP_mv1 <- subset(snoring_inhabitual_bmi_mv1, snoring_inhabitual_bmi_mv1$Exposure=="BMI" &
                                   snoring_inhabitual_bmi_mv1$Outcome=="SAP")[,-3]
BMI_inhabitual_SAP_mv1$Exposure=c(" ", " ")
BMI_inhabitual_UAP_mv1 <- subset(snoring_inhabitual_bmi_mv1, snoring_inhabitual_bmi_mv1$Exposure=="BMI" &
                                   snoring_inhabitual_bmi_mv1$Outcome=="UAP")[,-3]
BMI_inhabitual_UAP_mv1$Exposure=c(" ", " ")
BMI_inhabitual_CAD_mv1 <- subset(snoring_inhabitual_bmi_mv1, snoring_inhabitual_bmi_mv1$Exposure=="BMI" &
                                   snoring_inhabitual_bmi_mv1$Outcome=="CAD")[,-3]
BMI_inhabitual_CAD_mv1$Exposure=c(" ", " ")




#??װ?????뵼??
setwd("C:/Users/86132/Desktop/????/snoring2CAD_check/result/MVMR_raw")

#snoringΪ??¶
#CAD
snoring_CAD <- rbind(snoring_CAD_sv1, snoring_CAD_mv1, BMI_CAD_sv1, BMI_insnoring_CAD_mv1)
write.csv(snoring_CAD, "snoring_CAD.csv", row.names = FALSE)
#MI
snoring_MI <- rbind(snoring_MI_sv1, snoring_MI_mv1, BMI_MI_sv1, BMI_insnoring_MI_mv1)
write.csv(snoring_MI, "snoring_MI.csv", row.names = FALSE)
#CHF
snoring_CHF <- rbind(snoring_CHF_sv1, snoring_CHF_mv1, BMI_CHF_sv1, BMI_insnoring_CHF_mv1)
write.csv(snoring_CHF, "snoring_CHF.csv", row.names = FALSE)
#Angina
snoring_Angina <- rbind(snoring_Angina_sv1, snoring_Angina_mv1, BMI_Angina_sv1, BMI_insnoring_Angina_mv1)
write.csv(snoring_Angina, "snoring_Angina.csv", row.names = FALSE)
#UAP
snoring_UAP <- rbind(snoring_UAP_sv1, snoring_UAP_mv1, BMI_UAP_sv1, BMI_insnoring_UAP_mv1)
write.csv(snoring_UAP, "snoring_UAP.csv", row.names = FALSE)
#SAP
snoring_SAP <- rbind(snoring_SAP_sv1, snoring_SAP_mv1, BMI_SAP_sv1, BMI_insnoring_SAP_mv1)
write.csv(snoring_SAP, "snoring_SAP.csv", row.names = FALSE)


#habitualΪ??¶
#CAD
habitual_CAD <- rbind(habitual_CAD_sv1, habitual_CAD_mv1, BMI_CAD_sv1, BMI_inhabitual_CAD_mv1)
write.csv(habitual_CAD, "habitual_CAD.csv", row.names = FALSE)
#MI
habitual_MI <- rbind(habitual_MI_sv1, habitual_MI_mv1, BMI_MI_sv1, BMI_inhabitual_MI_mv1)
write.csv(habitual_MI, "habitual_MI.csv", row.names = FALSE)
#CHF
habitual_CHF <- rbind(habitual_CHF_sv1, habitual_CHF_mv1, BMI_CHF_sv1, BMI_inhabitual_CHF_mv1)
write.csv(habitual_CHF, "habitual_CHF.csv", row.names = FALSE)
#Angina
habitual_Angina <- rbind(habitual_Angina_sv1, habitual_Angina_mv1, BMI_Angina_sv1, BMI_inhabitual_Angina_mv1)
write.csv(habitual_Angina, "habitual_Angina.csv", row.names = FALSE)
#UAP
habitual_UAP <- rbind(habitual_UAP_sv1, habitual_UAP_mv1, BMI_UAP_sv1, BMI_inhabitual_UAP_mv1)
write.csv(habitual_UAP, "habitual_UAP.csv", row.names = FALSE)
#SAP
habitual_SAP <- rbind(habitual_SAP_sv1, habitual_SAP_mv1, BMI_SAP_sv1, BMI_inhabitual_SAP_mv1)
write.csv(habitual_SAP, "habitual_SAP.csv", row.names = FALSE)

#ע?⣬Ҫ??excel???ֶ??޸Ŀ?ѧ????????