########
#By Yunqing Zhu, Canqing Yu
#Email: zhuyun_qing@126.com yucanqing@pku.edu.cn
#########

#3.MVMR ????12??ɭ??ͼ
library(data.table)
library(readxl)
library(ggplot2)
library(readr)
library(tidyverse)
library(car)
library(forestplot)

setwd("C:/Users/86132/Desktop/复核/snoring2CAD_check/result/MVMR_raw")

styles <- fpShapesGp(
  lines = list(
    gpar(col = "black"), 
    gpar(col = "black"),    
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "dodgerblue"),    
    gpar(col = "dodgerblue"),
    gpar(col = "black"),    
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "black"),
    gpar(col = "dodgerblue"),    
    gpar(col = "dodgerblue")
  ),
  box = list(
    gpar(fill = "black"),  
    gpar(fill = "black"),    
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "dodgerblue"),    
    gpar(fill = "dodgerblue"),
    gpar(fill = "black"),    
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "black"),
    gpar(fill = "dodgerblue"),    
    gpar(fill = "dodgerblue")
  )
)

#snoring_CAD
dat <- read.csv("snoring_CAD.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                align=c("l","l","l","l"),
                mean=forest$OR,
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.8, 1.0,1.2, 1.4),
                clip=c(0.8,1.4),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p

#snoring_MI
dat <- read.csv("snoring_MI.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                mean=forest$OR,
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.8, 1.0, 1.4,1.8),
                clip=c(0.8,1.8),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p


#snoring_CHF
dat <- read.csv("snoring_CHF.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                mean=forest$OR,
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.8, 1.0,1.4,1.8),
                clip=c(0.8,1.8),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p



#snoring_Angina
dat <- read.csv("snoring_Angina.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                mean=forest$OR,
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.8, 1.0,1.2, 1.4),
                clip=c(0.8,1.4),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p



#snoring_UAP
dat <- read.csv("snoring_UAP.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                mean=forest$OR,
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.8, 1.0,1.4,1.8),
                clip=c(0.8,1.8),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p





#snoring_SAP
dat <- read.csv("snoring_SAP.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                mean=forest$OR,
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.8, 1.0,1.2, 1.4),
                clip=c(0.8,1.4),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p






#habitual_CAD
dat <- read.csv("habitual_CAD.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                mean=forest$OR,
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.8, 1.0,1.2, 1.4),
                clip=c(0.8,1.4),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p

#habitual_MI
dat <- read.csv("habitual_MI.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                mean=forest$OR,
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.8, 1.0,1.4, 1.8),
                clip=c(0.8,1.8),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p


#habitual_CHF
dat <- read.csv("habitual_CHF.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                mean=forest$OR,
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.8, 1.0,1.4,1.8),
                clip=c(0.8,1.8),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p



#habitual_Angina
dat <- read.csv("habitual_Angina.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                mean=forest$OR,
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.8, 1.0,1.2, 1.4),
                clip=c(0.8,1.4),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p



#habitual_UAP
dat <- read.csv("habitual_UAP.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                mean=forest$OR,
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.8, 1.0,1.4,1.8),
                clip=c(0.8,1.8),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p





#habitual_SAP
dat <- read.csv("habitual_SAP.csv",h=T)
dat$P <- ifelse(dat$P<0.001, "<0.001", sprintf("%0.3f", dat$P))
header <- c("Exposure", "Method", "OR (95%CI)","P", " ", " "," ")
forest <- rbind(header,dat)

p <- forestplot(labeltext=as.matrix(forest[,c(1,2,3,4)]),
                mean=forest$OR,
                
                lower=forest$OR_L,
                upper=forest$OR_U,
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
                xticks=c(0.6,0.8, 1.0,1.2, 1.4),
                clip=c(0.6,1.4),
                ci.vertices.height=0.2,
                zero=1,
                cex.lab=1.5,
                graphwidth=unit(50,'mm'),
                boxsize=0.3,
                lineheight=unit(7,'mm'),
                colgap=unit(5,'mm'),
                lwd.zero=1,
                lwd.ci=1.8,
                lwd.xaxis=2,
                xlog=FALSE,
                lty.ci=1,
                shapes_gp=styles,
                col=fpColors(box = 'black',summary = 'black',lines = 'black',zero = 'grey'),
                xlab="Odds ratio",
                txt_gp = fpTxtGp(label = gpar(cex=1),
                                 ticks = gpar(cex=1),
                                 xlab = gpar(cex=1)),
                hrzl_lines = gpar(col="black", lwd=2),
                graph.pos=3) 
p