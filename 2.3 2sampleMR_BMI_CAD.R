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
library(xtable)
library(flextable)
library(officer)
library(ggplot2)
library('gsmr')
library(mr.raps)
library(rio)

#画图函数
mr_scatter_plot1 <- function (mr_results, dat) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), 
                       function(d) {
                         d <- plyr::mutate(d)
                         if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         d <- subset(d, mr_keep)
                         index <- d$beta.exposure < 0
                         d$beta.exposure[index] <- d$beta.exposure[index] * 
                           -1
                         d$beta.outcome[index] <- d$beta.outcome[index] * 
                           -1
                         mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & 
                                           id.outcome == d$id.outcome[1])
                         mrres$a <- 0
                         if ("MR Egger" %in% mrres$method) {
                           temp <- mr_egger_regression(d$beta.exposure, 
                                                       d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                       default_parameters())
                           mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                         }
                         if ("MR Egger (bootstrap)" %in% mrres$method) {
                           temp <- mr_egger_regression_bootstrap(d$beta.exposure, 
                                                                 d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                                 default_parameters())
                           mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                         }
                         ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, 
                                                                y = beta.outcome)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - 
                                                                                                                           se.outcome, ymax = beta.outcome + se.outcome), 
                                                                                                            colour = "black", width = 0, lwd=0.4) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - 
                                                                                                                                                                                           se.exposure, xmax = beta.exposure + se.exposure), 
                                                                                                                                                                            colour = "black", height = 0) + ggplot2::geom_point(size=1.3) + 
                           ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, 
                                                                           slope = b, colour = method), lwd = 0.8, show.legend = TRUE) + 
                           ggplot2::scale_colour_manual(values = c("#1f78b4", "#a6cee3","#b2df8a", "#33a02c", 
                                                                   "#fb9a99", "#e31a1c", "#fdbf6f", 
                                                                   "#ff7f00", "#cab2d6", "#6a3d9a", 
                                                                   "#ffff99", "#b15928")) + ggplot2::labs(colour = "MR Test", 
                                                                                                          x = paste("SNP effect on", d$exposure[1]), 
                                                                                                          y = paste("SNP effect on", d$outcome[1])) + theme_bw() +
                           ggplot2::theme(legend.position = "top", 
                                          legend.direction = "vertical", panel.grid.major=element_blank(),
                                          panel.grid.minor=element_blank())  +  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
                       })
  mrres
}



mr_forest_plot1 <- 
  function (singlesnp_results, exponentiate = FALSE) 
  {
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("plyr", quietly = TRUE)
    res <- plyr::dlply(singlesnp_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 2) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
                                              levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
                                              am <- grep("All", d$SNP, value = TRUE)
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 0.01
                                              d$tot[d$SNP %in% am] <- 1
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[!d$SNP %in% am]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
                                              xint <- 0
                                              if (exponentiate) {
                                                d$b <- exp(d$b)
                                                d$up <- exp(d$up)
                                                d$lo <- exp(d$lo)
                                                xint <- 1
                                              }
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = xint, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot)), size=1.3) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("black", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.5, 
                                                                                                                                                                                                                       1)) + theme_bw() + 
                                                ggplot2::theme(legend.position = "none", 
                                                               axis.text.y = ggplot2::element_text(size = 10), axis.ticks.y = ggplot2::element_line(size = 0), 
                                                               axis.title.x = ggplot2::element_text(size = 11),panel.grid.major=element_blank(),
                                                               panel.grid.minor=element_blank()) + 
                                                ggplot2::labs(y = "", x = paste0("MR effect size for", " ",
                                                                                 d$exposure[1], " on ", d$outcome[1]))
                                            })
    res
  }


mr_leaveoneout_plot1 <- 
  function (leaveoneout_results) 
  {
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("plyr", quietly = TRUE)
    res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                              "id.outcome"), function(d) {
                                                d <- plyr::mutate(d)
                                                if (sum(!grepl("All", d$SNP)) < 3) {
                                                  return(blank_plot("Insufficient number of SNPs"))
                                                }
                                                d$up <- d$b + 1.96 * d$se
                                                d$lo <- d$b - 1.96 * d$se
                                                d$tot <- 1
                                                d$tot[d$SNP != "All"] <- 0.01
                                                d$SNP <- as.character(d$SNP)
                                                nom <- d$SNP[d$SNP != "All"]
                                                nom <- nom[order(d$b)]
                                                d <- rbind(d, d[nrow(d), ])
                                                d$SNP[nrow(d) - 1] <- ""
                                                d$b[nrow(d) - 1] <- NA
                                                d$up[nrow(d) - 1] <- NA
                                                d$lo[nrow(d) - 1] <- NA
                                                d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                                ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0,linetype = "dotted") + 
                                                  ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo,xmax = up, size = as.factor(tot), colour = as.factor(tot)), height = 0) + 
                                                  ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                  ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in%  "")), colour = "grey") + 
                                                  ggplot2::scale_colour_manual(values = c("black", "red")) + 
                                                  ggplot2::scale_size_manual(values = c(0.3, 1)) + 
                                                  theme_bw() +
                                                  ggplot2::theme(legend.position = "none",axis.text.y = ggplot2::element_text(size = 8), 
                                                                 axis.ticks.y = ggplot2::element_line(size = 0), axis.title.x = ggplot2::element_text(size = 10), 
                                                                 panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
                                                  ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n", 
                                                                                   d$exposure[1], " on ", d$outcome[1]))
                                              })
    res
  }


#第一次跑要生成BMI的IV数据
#CKB-BMI
ckb_bmi_summary <- fread("C:/Users/86132/Desktop/复核/snoring2CAD_check/data/bmi.summarydata.tsmr.txt")
setnames(ckb_bmi_summary, old = c("SNP","CHR","BP","P","A1","A2","BETA","SE","A1FREQ"), 
         new = c("SNP","CHR","BP","pval.exposure","EA","OA","beta","se","EAF"))
ckb_bmi_summary$exposure <- 'BMI'
ckb_bmi_summary$N=100640


#BMI-IVs
ckb_bmi_iv0 <- subset(ckb_bmi_summary, pval.exposure<5e-08 )

export(ckb_bmi_iv0, "C:/Users/86132/Desktop/复核/snoring2CAD_check/data/bmi.summarydata.tsmr1.txt")


######
#
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

#因为CKB的SNP???5个位点在BBJ的summary data中找不到，所以排除之，找代理位点
exp_dat$SNP[which(!exp_dat$SNP %in% CAD$SNP)]
#需要在初步拼接数据后进行此???

#"rs12885071" "rs6504568"  "rs35560038" "rs6044722"  "rs2744475"
#"rs2103785"  "rs28855509"

#进一步排除反向的
#rs13047416, rs1361511, rs2510032

BMI <- BMI[which(BMI$SNP!="rs12885071" & 
                   BMI$SNP!="rs6504568" &
                   BMI$SNP!="rs35560038" &
                   BMI$SNP!="rs6044722" &
                   BMI$SNP!="rs2744475" &
                   BMI$SNP!="rs2103785" &
                   BMI$SNP!="rs28855509" & 
                   BMI$SNP!="rs13047416" &
                   BMI$SNP!="rs1361511" &
                   BMI$SNP!="rs2510032"),]

exp_dat <- clump_data(BMI, clump_r2 = 0.001, 
                      clump_kb = 10000, pop = "EAS")

nrow(exp_dat)



#CAD
CAD <- read_outcome_data(snps = exp_dat$SNP,
                         filename = "C:/Users/86132/Desktop/复核/snoring2CAD_check/data/CAD_BBJ_tsmr.txt",
                         sep = "\t",                          
                         snp_col = "SNP",
                         beta_col = "b",
                         se_col = "se",
                         chr_col = "CHR",
                         pos_col = "BP",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         eaf_col = "freq",
                         pval_col = "P",
                         samplesize_col = "n"
)

CAD$outcome = "CAD"


dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = CAD)
nrow(dat)
dat$exposure = "BMI"



#利用MR PRESSO包做水平多效性和筛选outlier
#注意修改n和种子数
# mr_presso(BetaOutcome = "beta.outcome",
#           BetaExposure = "beta.exposure",
#           SdOutcome = "se.outcome",
#           SdExposure = "se.exposure",
#           OUTLIERtest = TRUE,
#           DISTORTIONtest = TRUE,
#           data = dat,
#           NbDistribution = 1300,
#           SignifThreshold = 0.05,
#           set.seed(100))

#剔除outlier IV后的结果
#剔除的IV: 

dat1 = dat
dat1 <- dat1[-c(31,46),]
dat = dat1



#Phenoscanner
write.table(dat[1], "C:/Users/86132/Desktop/复核/snoring2CAD_check/result/BMI2CAD.txt", 
            row.names=F, col.names=F, quote=F)

#phenoscanner扫描??? 排除2个与结局强关联的位点
dat <- subset(dat, SNP!="rs476828" & SNP!="rs6265")



#MR
set.seed(1000)
res <- mr(dat,method_list = c('mr_ivw', 
                              'mr_weighted_median',
                              'mr_egger_regression'))
res.raps <- mr.raps(b_exp = dat$beta.exposure,
                    b_out = dat$beta.outcome,
                    se_exp = dat$se.exposure,
                    se_out = dat$se.outcome)
res.raps.result <- c(res$id.exposure[1], res$id.outcome[1], res$outcome[1], res$exposure[1], "RAPS", res$nsnp[1], res.raps$beta.hat, res.raps$beta.se, res.raps$beta.p.value)
res <- rbind(res, res.raps.result)
res$b <- as.numeric(res$b)
res$se <- as.numeric(res$se)
res$nsnp <- as.numeric(res$nsnp)
res$method[1] <- "Inverse variance weighted"


#BMI不用adjust, 是每增加1个标准差
res$OR.adj <- exp(res$b)
res$OR.adj_L <- exp((res$b - 1.96*res$se))
res$OR.adj_U <- exp((res$b + 1.96*res$se))

result_CAD<-cbind(res$method, res$nsnp, res$OR.adj, res$OR.adj_L, res$OR.adj_U,res$pval)

colnames(result_CAD) <- c("method","nsnp","OR.adj","OR.adj_L","OR.adj_U","P")
rownames(result_CAD) <- res$outcome 


#水平多效性检???
#多效性检验，主要是检验多个IV是否存在水平多效性，
#常用MR Egger法的截距项表示，如果该截距项???0差异很大，说明存在水平多效???
pleio <- mr_pleiotropy_test(dat)
pleio


#TESTS
het <- mr_heterogeneity(dat)
het

#MR Steiger
out <- directionality_test(dat)
out
kable(out)


#F-Statistics
R2=out$snp_r2.exposure
n=dat$samplesize.exposure[1] #注意修改样本???
k=res$nsnp[1]
F = R2/(1-R2) * (n-k-1)/k
F

test_CAD <- cbind(pleio$exposure, pleio$outcome, F, pleio$egger_intercept, 
                  pleio$se, pleio$pval, het$Q[2], het$Q_pval[2], out$snp_r2.exposure,
                  out$snp_r2.outcome, out$correct_causal_direction, out$steiger_pval)
colnames(test_CAD) <- c("Exposure","Outcome","F","Egger.inter",
                        "Egger.SE","P.inter","Q","P.heterogeneity","snp_r2.exposure",
                        "snp_r2.outcome","Direction","P.Steiger")


setwd("C:/Users/86132/Desktop/复核/snoring2CAD_check/2sampleMR/plots/BMI")
single <- mr_leaveoneout(dat, method = mr_ivw) #default就是IVW
single
plot_leave1 <- mr_leaveoneout_plot1(single)
plot_leave1
ggsave(plot_leave1[[1]], file="BMI_CAD_leaveoneout.pdf", width=9, height=9)

p1<-mr_scatter_plot1(res, dat)
p1
ggsave(p1[[1]], file="BMI_CAD_scatterplot.pdf", width=7, height=7)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw"))
p2 <- mr_forest_plot1(res_single)
p2
ggsave(p2[[1]], file="BMI_CAD_singleplot.pdf", width=9, height=9)
p3 <- mr_funnel_plot(res_single)
p3




# summary(lm(dat$beta.outcome ~ dat$beta.exposure -1, weights = 1/dat$se.outcome^2))
# 
# ivw.res <- summary(lm(dat$beta.outcome ~ dat$beta.exposure -1, weights = 1/dat$se.outcome^2))
# b <- ivw.res$coef["dat$beta.exposure", "Estimate"]
# se <- ivw.res$coef["dat$beta.exposure", "Std. Error"]/ivw.res$sigma
# pval <- 2 * pnorm(abs(b/se), lower.tail = FALSE)



#b se pval
# 
# p1<-mr_scatter_plot(res, dat)
# p1
# res_single <- mr_singlesnp(dat, all_method = c("mr_ivw"))
# res_single <- mr_singlesnp(dat)
# p2 <- mr_forest_plot(res_single)
# p2
# p3 <- mr_funnel_plot(res_single)
# p3












#MI
#从这里开始跑
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

#因为CKB的SNP???5个位点在BBJ的summary data中找不到，所以排除之，找代理位点
exp_dat$SNP[which(!exp_dat$SNP %in% MI$SNP)]
#需要在初步拼接数据后进行此???

#"rs12885071" "rs6504568"  "rs35560038" "rs6044722"  "rs2744475"
#"rs2103785"  "rs28855509"

#进一步排除反向的
#rs13047416, rs1361511, rs2510032 rs6732947

#还有3个位点找不到 继续排除
#"rs2238689" "rs1532127" "rs4973506"

BMI <- BMI[which(BMI$SNP!="rs12885071" & 
                   BMI$SNP!="rs6504568" &
                   BMI$SNP!="rs35560038" &
                   BMI$SNP!="rs6044722" &
                   BMI$SNP!="rs2744475" &
                   BMI$SNP!="rs2103785" &
                   BMI$SNP!="rs28855509" & 
                   BMI$SNP!="rs13047416" &
                   BMI$SNP!="rs1361511" &
                   BMI$SNP!="rs2510032" & 
                   BMI$SNP!="rs2238689" &
                   BMI$SNP!="rs1532127" &
                   BMI$SNP!="rs4973506" &
                   BMI$SNP!="rs6732947" & 
                   BMI$SNP!="rs62136856"),]

exp_dat <- clump_data(BMI, clump_r2 = 0.001, 
                      clump_kb = 10000, pop = "EAS")

nrow(exp_dat)



MI <- read_outcome_data(snps = exp_dat$SNP,
                        filename = "C:/Users/86132/Desktop/复核/snoring2CAD_check/data/MI_BBJ_tsmr.txt",
                        sep = "\t",                          
                        snp_col = "SNP",
                        beta_col = "b",
                        se_col = "se",
                        chr_col = "CHR",
                        pos_col = "BP",
                        effect_allele_col = "A2",
                        other_allele_col = "A1",
                        eaf_col = "freq",
                        pval_col = "P",
                        samplesize_col = "n"
)

MI$outcome = "MI"

dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = MI)
nrow(dat)
dat$exposure = "BMI"





#利用MR PRESSO包做水平多效性和筛选outlier
#注意修改n和种子数
# mr_presso(BetaOutcome = "beta.outcome",
#           BetaExposure = "beta.exposure",
#           SdOutcome = "se.outcome",
#           SdExposure = "se.exposure",
#           OUTLIERtest = TRUE,
#           DISTORTIONtest = TRUE,
#           data = dat,
#           NbDistribution = 1300,
#           SignifThreshold = 0.05,
#           set.seed(100))

#剔除outlier IV后的结果

dat1 = dat
dat1 <- dat1[-c(4,30),]
dat = dat1




#phenoscanner
write.table(dat[1], "C:/Users/86132/Desktop/复核/snoring2CAD_check/result/BMI2MI.txt", 
            row.names=F, col.names=F, quote=F)

#phenoscanner扫描??? 排除2个与结局强关联的位点
dat <- subset(dat, SNP!="rs476828" & SNP!="rs6265")





set.seed(1000)
res <- mr(dat,method_list = c('mr_ivw',
                              'mr_weighted_median',
                              'mr_egger_regression'))
res.raps <- mr.raps(b_exp = dat$beta.exposure,
                    b_out = dat$beta.outcome,
                    se_exp = dat$se.exposure,
                    se_out = dat$se.outcome)
res.raps.result <- c(res$id.exposure[1], res$id.outcome[1], res$outcome[1], res$exposure[1], "RAPS", res$nsnp[1], res.raps$beta.hat, res.raps$beta.se, res.raps$beta.p.value)
res <- rbind(res, res.raps.result)
res$b <- as.numeric(res$b)
res$se <- as.numeric(res$se)
res$nsnp <- as.numeric(res$nsnp)
res$method[1] <- "Inverse variance weighted"


#导出结果 ln(1.5)，打鼾风险变为原来的1.5倍时...
res$OR.adj <- exp(res$b)
res$OR.adj_L <- exp((res$b - 1.96*res$se))
res$OR.adj_U <- exp((res$b + 1.96*res$se))

result_MI<-cbind(res$method, res$nsnp, res$OR.adj, res$OR.adj_L, res$OR.adj_U,res$pval)

colnames(result_MI) <- c("method","nsnp","OR.adj","OR.adj_L","OR.adj_U","P")
rownames(result_MI) <- res$outcome 


#水平多效性检???
#多效性检验，主要是检验多个IV是否存在水平多效性，
#常用MR Egger法的截距项表示，如果该截距项???0差异很大，说明存在水平多效???
pleio <- mr_pleiotropy_test(dat)
pleio


#TESTS
het <- mr_heterogeneity(dat)
het

#MR Steiger
out <- directionality_test(dat)
out
kable(out)


#F-Statistics
R2=out$snp_r2.exposure
n=dat$samplesize.exposure[1] #注意修改样本???
k=res$nsnp[1]
F = R2/(1-R2) * (n-k-1)/k
F

test_MI <- cbind(pleio$exposure, pleio$outcome, F, pleio$egger_intercept, 
                 pleio$se, pleio$pval, het$Q[2], het$Q_pval[2], out$snp_r2.exposure,
                 out$snp_r2.outcome, out$correct_causal_direction, out$steiger_pval)
colnames(test_MI) <- c("Exposure","Outcome","F","Egger.inter",
                       "Egger.SE","P.inter","Q","P.heterogeneity","snp_r2.exposure",
                       "snp_r2.outcome","Direction","P.Steiger")


setwd("C:/Users/86132/Desktop/复核/snoring2CAD_check/2sampleMR/plot/BMI")
single <- mr_leaveoneout(dat, method = mr_ivw) #default就是IVW
single
plot_leave1 <- mr_leaveoneout_plot1(single)
plot_leave1
ggsave(plot_leave1[[1]], file="BMI_MI_leaveoneout.pdf", width=9, height=9)

p1<-mr_scatter_plot1(res, dat)
p1
ggsave(p1[[1]], file="BMI_MI_scatterplot.pdf", width=7, height=7)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw"))
p2 <- mr_forest_plot1(res_single)
p2
ggsave(p2[[1]], file="BMI_MI_singleplot.pdf", width=9, height=9)
p3 <- mr_funnel_plot(res_single)
p3








#CHF
#从这里开始跑
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

#因为CKB的SNP???5个位点在BBJ的summary data中找不到，所以排除之，找代理位点
exp_dat$SNP[which(!exp_dat$SNP %in% CHF$SNP)]
#需要在初步拼接数据后进行此???

#"rs12885071" "rs6504568"  "rs35560038" "rs6044722"  "rs2744475"
#"rs2103785"  "rs28855509"

#进一步排除反向的
#rs13047416, rs1361511, rs2510032 rs6732947

#还有3个位点找不到 继续排除
#"rs2238689" "rs1532127" "rs4973506"

BMI <- BMI[which(BMI$SNP!="rs12885071" & 
                   BMI$SNP!="rs6504568" &
                   BMI$SNP!="rs35560038" &
                   BMI$SNP!="rs6044722" &
                   BMI$SNP!="rs2744475" &
                   BMI$SNP!="rs2103785" &
                   BMI$SNP!="rs28855509" & 
                   BMI$SNP!="rs13047416" &
                   BMI$SNP!="rs1361511" &
                   BMI$SNP!="rs2510032" & 
                   BMI$SNP!="rs2238689" &
                   BMI$SNP!="rs1532127" &
                   BMI$SNP!="rs4973506" &
                   BMI$SNP!="rs6732947" & 
                   BMI$SNP!="rs62136856"),]

exp_dat <- clump_data(BMI, clump_r2 = 0.001, 
                      clump_kb = 10000, pop = "EAS")

nrow(exp_dat)

CHF <- read_outcome_data(snps = exp_dat$SNP,
                         filename = "C:/Users/86132/Desktop/复核/snoring2CAD_check/data/CHF_BBJ_tsmr.txt",
                         sep = "\t",                          
                         snp_col = "SNP",
                         beta_col = "b",
                         se_col = "se",
                         chr_col = "CHR",
                         pos_col = "BP",
                         effect_allele_col = "A2",
                         other_allele_col = "A1",
                         eaf_col = "freq",
                         pval_col = "P",
                         samplesize_col = "n"
)

CHF$outcome = "CHF"

dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = CHF)
nrow(dat)
dat$exposure = "BMI"



#利用MR PRESSO包做水平多效性和筛选outlier
#注意修改n和种子数
# mr_presso(BetaOutcome = "beta.outcome",
#           BetaExposure = "beta.exposure",
#           SdOutcome = "se.outcome",
#           SdExposure = "se.exposure",
#           OUTLIERtest = TRUE,
#           DISTORTIONtest = TRUE,
#           data = dat,
#           NbDistribution = 1300,
#           SignifThreshold = 0.05,
#           set.seed(100))


dat <- dat[-c(7,14),]


#phenoscanner
write.table(dat[1], "C:/Users/86132/Desktop/复核/snoring2CAD_check/result/BMI2CHF.txt", 
            row.names=F, col.names=F, quote=F)

#phenoscanner扫描??? 排除2个与结局强关联的位点
dat <- subset(dat, SNP!="rs476828" & SNP!="rs6265")



#MR
set.seed(1000)
res <- mr(dat,method_list = c('mr_ivw',
                              'mr_weighted_median',
                              'mr_egger_regression'))
res.raps <- mr.raps(b_exp = dat$beta.exposure,
                    b_out = dat$beta.outcome,
                    se_exp = dat$se.exposure,
                    se_out = dat$se.outcome)
res.raps.result <- c(res$id.exposure[1], res$id.outcome[1], res$outcome[1], res$exposure[1], "RAPS", res$nsnp[1], res.raps$beta.hat, res.raps$beta.se, res.raps$beta.p.value)
res <- rbind(res, res.raps.result)
res$b <- as.numeric(res$b)
res$se <- as.numeric(res$se)
res$nsnp <- as.numeric(res$nsnp)
res$method[1] <- "Inverse variance weighted"


#导出结果 ln(1.5)，打鼾风险变为原来的1.5倍时...
res$OR.adj <- exp(res$b)
res$OR.adj_L <- exp((res$b - 1.96*res$se))
res$OR.adj_U <- exp((res$b + 1.96*res$se))

result_CHF<-cbind(res$method, res$nsnp, res$OR.adj, res$OR.adj_L, res$OR.adj_U,res$pval)

colnames(result_CHF) <- c("method","nsnp","OR.adj","OR.adj_L","OR.adj_U","P")
rownames(result_CHF) <- res$outcome 


#水平多效性检???
#多效性检验，主要是检验多个IV是否存在水平多效性，
#常用MR Egger法的截距项表示，如果该截距项???0差异很大，说明存在水平多效???
pleio <- mr_pleiotropy_test(dat)
pleio


#TESTS
het <- mr_heterogeneity(dat)
het

#MR Steiger
out <- directionality_test(dat)
out
kable(out)


#F-Statistics
R2=out$snp_r2.exposure
n=dat$samplesize.exposure[1] #注意修改样本???
k=res$nsnp[1]
F = R2/(1-R2) * (n-k-1)/k
F

test_CHF <- cbind(pleio$exposure, pleio$outcome, F, pleio$egger_intercept, 
                  pleio$se, pleio$pval, het$Q[2], het$Q_pval[2], out$snp_r2.exposure,
                  out$snp_r2.outcome, out$correct_causal_direction, out$steiger_pval)
colnames(test_CHF) <- c("Exposure","Outcome","F","Egger.inter",
                        "Egger.SE","P.inter","Q","P.heterogeneity","snp_r2.exposure",
                        "snp_r2.outcome","Direction","P.Steiger")


setwd("C:/Users/86132/Desktop/复核/snoring2CAD_check/2sampleMR/plots/BMI")
single <- mr_leaveoneout(dat, method = mr_ivw) #default就是IVW
single
plot_leave1 <- mr_leaveoneout_plot1(single)
plot_leave1
ggsave(plot_leave1[[1]], file="BMI_CHF_leaveoneout.pdf", width=9, height=9)

p1<-mr_scatter_plot1(res, dat)
p1
ggsave(p1[[1]], file="BMI_CHF_scatterplot.pdf", width=7, height=7)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw"))
p2 <- mr_forest_plot1(res_single)
p2
ggsave(p2[[1]], file="BMI_CHF_singleplot.pdf", width=9, height=9)
p3 <- mr_funnel_plot(res_single)
p3








#Angina
Angina <- read_outcome_data(snps = exp_dat$SNP,
                            filename = "C:/Users/86132/Desktop/复核/snoring2CAD_check/data/Angina_BBJ_tsmr.txt",
                            sep = "\t",                          
                            snp_col = "SNP",
                            beta_col = "b",
                            se_col = "se",
                            chr_col = "CHR",
                            pos_col = "BP",
                            effect_allele_col = "ALT",
                            other_allele_col = "REF",
                            eaf_col = "freq",
                            pval_col = "P",
                            samplesize_col = "n"
)

Angina$outcome = "Angina"

dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = Angina)
nrow(dat)
dat$exposure = "BMI"



#利用MR PRESSO包做水平多效性和筛选outlier
#注意修改n和种子数
# mr_presso(BetaOutcome = "beta.outcome",
#           BetaExposure = "beta.exposure",
#           SdOutcome = "se.outcome",
#           SdExposure = "se.exposure",
#           OUTLIERtest = TRUE,
#           DISTORTIONtest = TRUE,
#           data = dat,
#           NbDistribution = 1300,
#           SignifThreshold = 0.05,
#           set.seed(100))

#剔除outlier IV后的结果

dat1 = dat
dat1 <- dat1[-30,]
dat = dat1
#dat1 <- dat1[-c(11,22),]


#phenoscanner
#write.table(dat[1], "C:/Users/86132/Desktop/复核/snoring2CAD_check/result/BMI2CHF.txt", 
#           row.names=F, col.names=F, quote=F)

#phenoscanner扫描??? 排除2个与结局强关联的位点
dat <- subset(dat, SNP!="rs476828" & SNP!="rs6265")



#MR
set.seed(1000)
res <- mr(dat,method_list = c('mr_ivw',
                              'mr_weighted_median',
                              'mr_egger_regression'))
res.raps <- mr.raps(b_exp = dat$beta.exposure,
                    b_out = dat$beta.outcome,
                    se_exp = dat$se.exposure,
                    se_out = dat$se.outcome)
res.raps.result <- c(res$id.exposure[1], res$id.outcome[1], res$outcome[1], res$exposure[1], "RAPS", res$nsnp[1], res.raps$beta.hat, res.raps$beta.se, res.raps$beta.p.value)
res <- rbind(res, res.raps.result)
res$b <- as.numeric(res$b)
res$se <- as.numeric(res$se)
res$nsnp <- as.numeric(res$nsnp)
res$method[1] <- "Inverse variance weighted"



res$OR.adj <- exp(res$b)
res$OR.adj_L <- exp((res$b - 1.96*res$se))
res$OR.adj_U <- exp((res$b + 1.96*res$se))

result_Angina<-cbind(res$method, res$nsnp, res$OR.adj, res$OR.adj_L, res$OR.adj_U,res$pval)

colnames(result_Angina) <- c("method","nsnp","OR.adj","OR.adj_L","OR.adj_U","P")
rownames(result_Angina) <- res$outcome 



#水平多效性检???
#多效性检验，主要是检验多个IV是否存在水平多效性，
#常用MR Egger法的截距项表示，如果该截距项???0差异很大，说明存在水平多效???
pleio <- mr_pleiotropy_test(dat)
pleio

#TESTS
het <- mr_heterogeneity(dat)
het

#MR Steiger
out <- directionality_test(dat)
out
kable(out)


#F-Statistics
R2=out$snp_r2.exposure
n=dat$samplesize.exposure[1] #注意修改样本???
k=res$nsnp[1]
F = R2/(1-R2) * (n-k-1)/k
F

test_Angina <- cbind(pleio$exposure, pleio$outcome, F, pleio$egger_intercept, 
                     pleio$se, pleio$pval, het$Q[2], het$Q_pval[2], out$snp_r2.exposure,
                     out$snp_r2.outcome, out$correct_causal_direction, out$steiger_pval)
colnames(test_Angina) <- c("Exposure","Outcome","F","Egger.inter",
                           "Egger.SE","P.inter","Q","P.heterogeneity","snp_r2.exposure",
                           "snp_r2.outcome","Direction","P.Steiger")


setwd("C:/Users/86132/Desktop/复核/snoring2CAD_check/2sampleMR/plotsBMI")
single <- mr_leaveoneout(dat, method = mr_ivw) #default就是IVW
single
plot_leave1 <- mr_leaveoneout_plot1(single)
plot_leave1
ggsave(plot_leave1[[1]], file="BMI_Angina_leaveoneout.pdf", width=9, height=9)

p1<-mr_scatter_plot1(res, dat)
p1
ggsave(p1[[1]], file="BMI_Angina_scatterplot.pdf", width=7, height=7)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw"))
p2 <- mr_forest_plot1(res_single)
p2
ggsave(p2[[1]], file="BMI_Angina_singleplot.pdf", width=9, height=9)
p3 <- mr_funnel_plot(res_single)
p3









#UAP
UAP <- read_outcome_data(snps = exp_dat$SNP,
                         filename = "C:/Users/86132/Desktop/复核/snoring2CAD_check/data/UAP_BBJ_tsmr.txt",
                         sep = "\t",                          
                         snp_col = "SNP",
                         beta_col = "b",
                         se_col = "se",
                         chr_col = "CHR",
                         pos_col = "BP",
                         effect_allele_col = "ALT",
                         other_allele_col = "REF",
                         eaf_col = "freq",
                         pval_col = "P",
                         samplesize_col = "n"
)

UAP$outcome = "UAP"

dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = UAP)
nrow(dat)
dat$exposure = "BMI"



#利用MR PRESSO包做水平多效性和筛选outlier
#注意修改n和种子数
# mr_presso(BetaOutcome = "beta.outcome",
#           BetaExposure = "beta.exposure",
#           SdOutcome = "se.outcome",
#           SdExposure = "se.exposure",
#           OUTLIERtest = TRUE,
#           DISTORTIONtest = TRUE,
#           data = dat,
#           NbDistribution = 1300,
#           SignifThreshold = 0.05,
#           set.seed(100))

#剔除outlier IV后的结果

dat1 = dat
dat1 <- dat1[-30,]
dat = dat1

#phenoscanner
#write.table(dat[1], "C:/Users/86132/Desktop/复核/snoring2CAD_check/result/BMI2CHF.txt", 
#           row.names=F, col.names=F, quote=F)

#phenoscanner扫描??? 排除2个与结局强关联的位点
dat <- subset(dat, SNP!="rs476828" & SNP!="rs6265")




#MR
set.seed(1000)
res <- mr(dat,method_list = c('mr_ivw_mre',
                              'mr_weighted_median',
                              'mr_egger_regression'))
res.raps <- mr.raps(b_exp = dat$beta.exposure,
                    b_out = dat$beta.outcome,
                    se_exp = dat$se.exposure,
                    se_out = dat$se.outcome)
res.raps.result <- c(res$id.exposure[1], res$id.outcome[1], res$outcome[1], res$exposure[1], "RAPS", res$nsnp[1], res.raps$beta.hat, res.raps$beta.se, res.raps$beta.p.value)
res <- rbind(res, res.raps.result)
res$b <- as.numeric(res$b)
res$se <- as.numeric(res$se)
res$nsnp <- as.numeric(res$nsnp)
res$method[1] <- "Inverse variance weighted"


res$OR.adj <- exp(res$b)
res$OR.adj_L <- exp((res$b - 1.96*res$se))
res$OR.adj_U <- exp((res$b + 1.96*res$se))

result_UAP<-cbind(res$method, res$nsnp, res$OR.adj, res$OR.adj_L, res$OR.adj_U,res$pval)

colnames(result_UAP) <- c("method","nsnp","OR.adj","OR.adj_L","OR.adj_U","P")
rownames(result_UAP) <- res$outcome 


#水平多效性检???
#多效性检验，主要是检验多个IV是否存在水平多效性，
#常用MR Egger法的截距项表示，如果该截距项???0差异很大，说明存在水平多效???
pleio <- mr_pleiotropy_test(dat)
pleio

#TESTS
het <- mr_heterogeneity(dat)
het

#MR Steiger
out <- directionality_test(dat)
out
kable(out)


#F-Statistics
R2=out$snp_r2.exposure
n=dat$samplesize.exposure[1] #注意修改样本???
k=res$nsnp[1]
F = R2/(1-R2) * (n-k-1)/k
F

test_UAP <- cbind(pleio$exposure, pleio$outcome, F, pleio$egger_intercept, 
                  pleio$se, pleio$pval, het$Q[2], het$Q_pval[2], out$snp_r2.exposure,
                  out$snp_r2.outcome, out$correct_causal_direction, out$steiger_pval)
colnames(test_UAP) <- c("Exposure","Outcome","F","Egger.inter",
                        "Egger.SE","P.inter","Q","P.heterogeneity","snp_r2.exposure",
                        "snp_r2.outcome","Direction","P.Steiger")


setwd("C:/Users/86132/Desktop/复核/snoring2CAD_check/2sampleMR/plots/BMI")
single <- mr_leaveoneout(dat, method = mr_ivw) #default就是IVW
single
plot_leave1 <- mr_leaveoneout_plot1(single)
plot_leave1
ggsave(plot_leave1[[1]], file="BMI_UAP_leaveoneout.pdf", width=9, height=9)

p1<-mr_scatter_plot1(res, dat)
p1
ggsave(p1[[1]], file="BMI_UAP_scatterplot.pdf", width=7, height=7)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw"))
p2 <- mr_forest_plot1(res_single)
p2
ggsave(p2[[1]], file="BMI_UAP_singleplot.pdf", width=9, height=9)
p3 <- mr_funnel_plot(res_single)
p3








#SAP
SAP <- read_outcome_data(snps = exp_dat$SNP,
                         filename = "C:/Users/86132/Desktop/复核/snoring2CAD_check/data/SAP_BBJ_tsmr.txt",
                         sep = "\t",                          
                         snp_col = "SNP",
                         beta_col = "b",
                         se_col = "se",
                         chr_col = "CHR",
                         pos_col = "BP",
                         effect_allele_col = "ALT",
                         other_allele_col = "REF",
                         eaf_col = "freq",
                         pval_col = "P",
                         samplesize_col = "n"
)

SAP$outcome = "SAP"

dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = SAP)
nrow(dat)
dat$exposure = "BMI"



#利用MR PRESSO包做水平多效性和筛选outlier
#注意修改n和种子数
# mr_presso(BetaOutcome = "beta.outcome",
#           BetaExposure = "beta.exposure",
#           SdOutcome = "se.outcome",
#           SdExposure = "se.exposure",
#           OUTLIERtest = TRUE,
#           DISTORTIONtest = TRUE,
#           data = dat,
#           NbDistribution = 1300,
#           SignifThreshold = 0.05,
#           set.seed(100))

#剔除outlier IV后的结果

dat1 = dat
dat1 <- dat1[-c(13,30),]
dat = dat1
#dat1 <- dat1[-c(11,22),]



#phenoscanner
#write.table(dat[1], "C:/Users/86132/Desktop/复核/snoring2CAD_check/result/BMI2CHF.txt", 
#           row.names=F, col.names=F, quote=F)

#phenoscanner扫描??? 排除2个与结局强关联的位点
dat <- subset(dat, SNP!="rs476828" & SNP!="rs6265")









#MR
set.seed(1000)
res <- mr(dat,method_list = c('mr_ivw',
                              'mr_weighted_median',
                              'mr_egger_regression'))
res.raps <- mr.raps(b_exp = dat$beta.exposure,
                    b_out = dat$beta.outcome,
                    se_exp = dat$se.exposure,
                    se_out = dat$se.outcome)
res.raps.result <- c(res$id.exposure[1], res$id.outcome[1], res$outcome[1], res$exposure[1], "RAPS", res$nsnp[1], res.raps$beta.hat, res.raps$beta.se, res.raps$beta.p.value)
res <- rbind(res, res.raps.result)
res$b <- as.numeric(res$b)
res$se <- as.numeric(res$se)
res$nsnp <- as.numeric(res$nsnp)
res$method[1] <- "Inverse variance weighted"


res$OR.adj <- exp(res$b)
res$OR.adj_L <- exp((res$b - 1.96*res$se))
res$OR.adj_U <- exp((res$b + 1.96*res$se))

result_SAP<-cbind(res$method, res$nsnp, res$OR.adj, res$OR.adj_L, res$OR.adj_U,res$pval)

colnames(result_SAP) <- c("method","nsnp","OR.adj","OR.adj_L","OR.adj_U","P")
rownames(result_SAP) <- res$outcome 


#水平多效性检???
#多效性检验，主要是检验多个IV是否存在水平多效性，
#常用MR Egger法的截距项表示，如果该截距项???0差异很大，说明存在水平多效???
pleio <- mr_pleiotropy_test(dat)
pleio

#TESTS
het <- mr_heterogeneity(dat)
het

#MR Steiger
out <- directionality_test(dat)
out
kable(out)


#F-Statistics
R2=out$snp_r2.exposure
n=dat$samplesize.exposure[1] #注意修改样本???
k=res$nsnp[1]
F = R2/(1-R2) * (n-k-1)/k
F

test_SAP <- cbind(pleio$exposure, pleio$outcome, F, pleio$egger_intercept, 
                  pleio$se, pleio$pval, het$Q[2], het$Q_pval[2], out$snp_r2.exposure,
                  out$snp_r2.outcome, out$correct_causal_direction, out$steiger_pval)
colnames(test_SAP) <- c("Exposure","Outcome","F","Egger.inter",
                        "Egger.SE","P.inter","Q","P.heterogeneity","snp_r2.exposure",
                        "snp_r2.outcome","Direction","P.Steiger")


setwd("C:/Users/86132/Desktop/复核/snoring2CAD_check/2sampleMR/plots/BMI")
single <- mr_leaveoneout(dat, method = mr_ivw) #default就是IVW
single
plot_leave1 <- mr_leaveoneout_plot1(single)
plot_leave1
ggsave(plot_leave1[[1]], file="BMI_SAP_leaveoneout.pdf", width=9, height=9)

p1<-mr_scatter_plot1(res, dat)
p1
ggsave(p1[[1]], file="BMI_SAP_scatterplot.pdf", width=7, height=7)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw"))
p2 <- mr_forest_plot1(res_single)
p2
ggsave(p2[[1]], file="BMI_SAP_singleplot.pdf", width=9, height=9)
p3 <- mr_funnel_plot(res_single)
p3


#合并结果
result <- rbind(result_CAD, result_MI, result_CHF, result_Angina, result_SAP, result_UAP)
test <- rbind(test_CAD, test_MI, test_CHF, test_Angina, test_SAP, test_UAP)

setwd("C:/Users/86132/Desktop/复核/snoring2CAD_check/2sampleMR")
write.csv(result, "57snps_BMI_CAD.csv")
write.csv(test, "57snps_BMI_CAD_test.csv")