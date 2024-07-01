########
#By Yunqing Zhu, Canqing Yu
#Email: zhuyun_qing@126.com yucanqing@pku.edu.cn
#########

library(data.table)
setwd("C:/Users/86132/Desktop/复核/snoring2CAD_check/2sampleMR")

#snoring
snoring <- read.csv("3snps_snoring_CAD.csv")
colnames(snoring) <- c("Outcomes", "Method", "N.SNPs", "OR", "OR_L", "OR_U", "P")
snoring <- snoring[which(snoring$Method=="Inverse variance weighted"),]
snoring$ORCI <- paste0(round(snoring$OR,2)," (", round(snoring$OR_L,2), ",", round(snoring$OR_U,2), ")")
snoring$P <- sprintf("%0.3f",snoring$P)
snoring <- snoring[,c(1,3,8,7)]
header_snoring <- c("Snoring", " ", " ", " ")
snoring1 <- rbind(header_snoring, snoring)

#habitual snoring
habitual <- read.csv("3snps_habitual_CAD.csv")
colnames(habitual) <- c("Outcomes", "Method", "N.SNPs", "OR", "OR_L", "OR_U", "P")
habitual <- habitual[which(habitual$Method=="Inverse variance weighted"),]
habitual$ORCI <- paste0(round(habitual$OR,2)," (", round(habitual$OR_L,2), ",", round(habitual$OR_U,2), ")")
habitual$P <- sprintf("%0.3f",habitual$P)
habitual <- habitual[,c(1,3,8,7)]
header_habitual <- c("Habitual snoring", " ", " ", " ")
habitual1 <- rbind(header_habitual, habitual)

#BMI
BMI <- read.csv("57snps_BMI_CAD.csv")
colnames(BMI) <- c("Outcomes", "Method", "N.SNPs", "OR", "OR_L", "OR_U", "P")
BMI <- BMI[which(BMI$Method=="Inverse variance weighted"),]
BMI$ORCI <- paste0(round(BMI$OR,2)," (", round(BMI$OR_L,2), ",", round(BMI$OR_U,2), ")")
BMI$P <- signif(BMI$P,3)
BMI <- BMI[,c(1,3,8,7)]
header_BMI <- c("BMI", " ", " ", " ")
BMI1 <- rbind(header_BMI, BMI)

result <- rbind(snoring1, habitual1, BMI1)

write.csv(result,"C:/Users/86132/Desktop/复核/snoring2CAD_check/result/1.SVMR.csv")
