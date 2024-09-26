############################## Functions and libraries ##############################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
#library(ConsensusClusterPlus) # R3.xx and R4.xx have different versions of ConsensusClusterPlus
#library(preprocessCore)
library(openxlsx)
#library(decoderr) 
#library(stringr)
library(survival)
library(survminer)
#library(dplyr)
#library(tidyverse)
#library(tidytidbits)
library(survivalAnalysis)

# load functions
file.sources <- list.files("./R/",pattern="*.R")
file.sources <- paste("./R/", file.sources, sep="")
sapply(file.sources, source)

# load subtype info
### This is a combined subtype object with
### subtypeColList, subtypeGeneList, subtypeList and schemaList
load("cmbSubtypes.230823.RData")

############################## Parsing ##############################
# load data
dataSet <- readRDS("RASH_ACCEPT.rds")
sampSub <- 1:nrow(dataSet$sampInfo)

studyList <- c("ACCEPT","RASH")
studyColList <- c("#91d1c2","#b09c8f")

armList <- c("Gemcitabine/Afatinib","Gemcitabine","No_further","Gemcitabine/Erlotinib","FOLFIRINOX")
armColList <- c("#9f79ee","#9acd32","#ffda1f","#e36ca3","#52BBD3")


# Suppl Fig 6a  ---------------------------------------------------------------------------------------------------------------------------------------------
# heatmap

ColSideColors <- data.frame(PurIST_score = colorRampPalette(c('snow2', 'black'))(length(dataSet$sampInfo$PurIST.prob))[rank(dataSet$sampInfo$PurIST.prob)],
                            PurIST = c("orange","blue")[as.factor(dataSet$sampInfo$PurIST)], 
                            Arm = armColList[as.factor(dataSet$clinicalDat$arm_all)],
                            Study = studyColList[as.factor(dataSet$clinicalDat$study)],
                            stringsAsFactors = FALSE)

### purist
TSPgenes <- subtypeGeneList[[1]]
featSet <- match(TSPgenes$geneSymbol,dataSet$featInfo$SYMBOL)
TSPgeneMat <- dataSet$ex[featSet, sampSub]
sampOrd <- order(dataSet$sampInfo$PurIST.prob[sampSub], decreasing = TRUE)
TSPgeneMat <- TSPgeneMat[, sampOrd]
TSPgeneMat <- apply(TSPgeneMat, 2, rank)

### set parameters
mainLab <- "PurIST TSP"
x <- TSPgeneMat
ColSideColors <- ColSideColors[sampOrd,]
ColSideColorsSize <- 3
geneInfo <- TSPgenes

heatmap.3(main= mainLab
          ,x = as.matrix(x)
          ,Rowv = FALSE
          ,Colv = FALSE
          ,dendrogram = "none" ,trace="none" ,scale="row"
          ,labRow = geneInfo$geneSymbol
          ,labCol = colnames(x)
          ,col = colorRampPalette(c("darkblue","blue","white","red","darkred"))(99)
          ,ColSideColors = as.matrix(ColSideColors)
          ,ColSideColorsSize = ColSideColorsSize
          ,RowSideColors = t(as.matrix(geneInfo$Color))
          ,RowSideColorsSize=1
          ,lwid = c(1,5)
          ,lhei = c(1,3)
          ,margins =c(15,20)
          ,keysize=0.5
          ,cex.lab=0.4
          ,cexRow = 1.1
          ,cexCol = 0.6
)
legend(xy.coords(x=.7,y=.98),
       legend=c("Study",studyList[c(2,1)],"Arm",armList[c(3,4,5,1,2)],
                "PurIST","Basal-like","Classical",
                "PurIST_score", "High","Low"),
       fill=c("white",studyColList[c(2,1)],
              "white",armColList[c(3,4,5,1,2)],
              "white","orange","blue",
              "white","black","snow2"),
       border=FALSE, bty="n",
       y.intersp = 0.85 , cex = 0.8)

# survival plots ---------------------------------------------------------------------------------------------------------------------------------------------

survDat <- dataSet$clinicalDat_Parsed

# Suppl Fig 6h
s="Accept"
survDatTmp <- survDat[which(survDat$study==s),]
factorList <- levels(factor(survDatTmp$arm_all))
survDatObj <- survDatTmp
names(survDatObj)[which(names(survDatObj) %in% "arm_all")] <- "obj"
km_fit <- survfit(Surv(survDatObj[,"os_time"], survDatObj[, "os"]) ~ obj, data = survDatObj, type = "kaplan-meier")
df <- data.frame(surv_median(km_fit))
textTmp <- list()
for (i in 1:nrow(df)) {
  textTmp[[i]] <- paste(gsub("obj=","",df$strata[i]), ":",
                        round(df$median[i],2)," (",
                        round(df$lower[i],2),",",
                        round(df$upper[i],2),")",sep="")
}
textPrt <- as.character(paste(unlist(textTmp), sep="", collapse="\n"))

ggsurv <- ggsurvplot(km_fit, 
                     data = survDatObj, 
                     break.time.by = 6,
                     title = "", 
                     #subtitle = textPrt,
                     conf.int = F, 
                     risk.table = T,
                     risk.table.height = 0.35,
                     surv.median.line = "hv", 
                     size = 0.7, censor.size=5, 
                     palette = armColList[c(1,2)],
                     legend.title="ACCEPT", 
                     legend.labs=factorList, # legend
                     xlab = "", ylab="Survival probability", xlim = c(0,30),
                     pval = T, pval.size = 9, pval.coord = c(22,0.75),  # p-values
                     font.main = c(23, "plain", "black"), 
                     font.legend = c(19, "plain"),
                     font.x = c(23, "plain", "black"),
                     font.y = c(23, "plain", "black"),
                     font.tickslab = c(23, "plain", "black"), ) + 
  guides(colour = guide_legend(nrow = 2)) 
ggsurv$table <- ggrisktable(km_fit, 
                            data = survDatObj, 
                            break.time.by = 6,
                            color = "strata", 
                            palette = armColList[c(1,2)],
                            y.text = F,  
                            fontsize = 8,
                            xlab = "Time (months)", ylab = "",xlim = c(0,30),
                            tables.theme = theme_survminer(font.main = 23, font.tickslab = 23, font.x = 23))
print(ggsurv)

# Fig 6g
s="Accept"
survDatTmp <- survDat[which(survDat$study==s),]
factorList <- levels(factor(survDatTmp$arm_all))
survDatObj <- survDatTmp
survDatObj$arm_subtype <- paste(survDatObj$arm_all, 
                                survDatObj[,"PurIST"], 
                                sep = "_")
survDatObj$arm_subtype <- factor(survDatObj$arm_subtype, 
                                 levels = c("Gemcitabine/Afatinib_Basal-like","Gemcitabine/Afatinib_Classical",
                                            "Gemcitabine_Basal-like","Gemcitabine_Classical"))
names(survDatObj)[which(names(survDatObj) %in% "arm_subtype")] <- "obj"
km_fit <- survfit(Surv(survDatObj[,"os_time"], survDatObj[, "os"]) ~ obj, data = survDatObj, type = "kaplan-meier")
df <- data.frame(surv_median(km_fit))
textTmp <- list()
for (i in 1:nrow(df)) {
  textTmp[[i]] <- paste(gsub("obj=","",df$strata[i]), ":",
                        round(df$median[i],2)," (",
                        round(df$lower[i],2),",",
                        round(df$upper[i],2),")",sep="")
}
textPrt <- as.character(paste(unlist(textTmp), sep="", collapse="\n"))

ggsurv <- ggsurvplot(km_fit, 
                     data = survDatObj, 
                     break.time.by = 6,
                     title = "", 
                     #subtitle = textPrt,
                     conf.int = F, 
                     risk.table = T,
                     risk.table.height = 0.35,
                     surv.median.line = "hv", 
                     size = 0.7, censor.size=5, 
                     palette = paletteer::paletteer_d("ggthemes::excel_Headlines")[c(2,5,3,1)],
                     legend.title="ACCEPT", 
                     legend.labs=c("G/A:Basal-like","G/A:Classical","G:Basal-like","G:Classical"), # legend
                     xlab = "", ylab="Survival probability", xlim = c(0,30),
                     pval = T, pval.size = 9, pval.coord = c(22,0.75),  # p-values
                     font.main = c(23, "plain", "black"), 
                     font.legend = c(19, "plain"),
                     font.x = c(23, "plain", "black"),
                     font.y = c(23, "plain", "black"),
                     font.tickslab = c(23, "plain", "black") ) + 
  guides(colour = guide_legend(nrow = 2)) 
ggsurv$table <- ggrisktable(km_fit, 
                            data = survDatObj, 
                            break.time.by = 6,
                            color = "strata", 
                            palette = paletteer::paletteer_d("ggthemes::excel_Headlines")[c(2,5,3,1)],
                            y.text = F,  
                            fontsize = 8,
                            xlab = "Time (months)", ylab = "", xlim = c(0,30),
                            tables.theme = theme_survminer(font.main = 23, font.tickslab = 23, font.x = 23))
print(ggsurv)

# Suppl Fig 6b
s="Rash"
survDatTmp <- survDat[which(survDat$study==s & survDat$arm_all != "No_further"),]
factorList <- levels(factor(survDatTmp$arm_all))
survDatObj <- survDatTmp
names(survDatObj)[which(names(survDatObj) %in% "arm_all")] <- "obj"
km_fit <- survfit(Surv(survDatObj[,"os_time"], survDatObj[, "os"]) ~ obj, data = survDatObj, type = "kaplan-meier")
df <- data.frame(surv_median(km_fit))
textTmp <- list()
for (i in 1:nrow(df)) {
  textTmp[[i]] <- paste(gsub("obj=","",df$strata[i]), ":",
                        round(df$median[i],2)," (",
                        round(df$lower[i],2),",",
                        round(df$upper[i],2),")",sep="")
}
textPrt <- as.character(paste(unlist(textTmp), sep="", collapse="\n"))

ggsurv <- ggsurvplot(km_fit, 
                     data = survDatObj, 
                     break.time.by = 6,
                     title = "", 
                     #subtitle = textPrt,
                     conf.int = F, 
                     risk.table = T,
                     risk.table.height = 0.35,
                     surv.median.line = "hv", 
                     size = 0.7, censor.size=5, 
                     palette = armColList[c(4,5)],
                     legend.title="RASH", 
                     legend.labs=factorList, # legend
                     xlab = "", ylab="Survival probability", xlim = c(0,30),
                     pval = T, pval.size = 9, pval.coord = c(22,0.75),  # p-values
                     font.main = c(23, "plain", "black"), 
                     font.legend = c(19, "plain"),
                     font.x = c(23, "plain", "black"),
                     font.y = c(23, "plain", "black"),
                     font.tickslab = c(23, "plain", "black"), ) + 
  guides(colour = guide_legend(nrow = 2)) 
ggsurv$table <- ggrisktable(km_fit, 
                            data = survDatObj, 
                            break.time.by = 6,
                            color = "strata", 
                            palette = armColList[c(4,5)],
                            y.text = F,  
                            fontsize = 8,
                            xlab = "Time (months)", ylab = "", xlim = c(0,30),
                            tables.theme = theme_survminer(font.main = 23, font.tickslab = 23, font.x = 23))
print(ggsurv)

# Suppl Fig 6f 
s="Rash"
survDatTmp <- survDat[which(survDat$study==s & survDat$arm_all != "No_further"),]
factorList <- levels(factor(survDatTmp$arm_all))
survDatObj <- survDatTmp
survDatObj$arm_subtype <- paste(survDatObj$arm_all, 
                                survDatObj[,"PurIST"], 
                                sep = "_")
survDatObj$arm_subtype <- factor(survDatObj$arm_subtype, 
                                 levels = c("Gemcitabine/Erlotinib_Basal-like","Gemcitabine/Erlotinib_Classical",
                                 "FOLFIRINOX_Basal-like","FOLFIRINOX_Classical"))
names(survDatObj)[which(names(survDatObj) %in% "arm_subtype")] <- "obj"
km_fit <- survfit(Surv(survDatObj[,"os_time"], survDatObj[, "os"]) ~ obj, data = survDatObj, type = "kaplan-meier")
df <- data.frame(surv_median(km_fit))
textTmp <- list()
for (i in 1:nrow(df)) {
  textTmp[[i]] <- paste(gsub("obj=","",df$strata[i]), ":",
                        round(df$median[i],2)," (",
                        round(df$lower[i],2),",",
                        round(df$upper[i],2),")",sep="")
}
textPrt <- as.character(paste(unlist(textTmp), sep="", collapse="\n"))

ggsurv <- ggsurvplot(km_fit, 
                     data = survDatObj, 
                     break.time.by = 6,
                     title = "", 
                     #subtitle = textPrt,
                     conf.int = F, 
                     risk.table = T,
                     risk.table.height = 0.35,
                     surv.median.line = "hv", 
                     size = 0.7, censor.size=5, 
                     palette = paletteer::paletteer_d("ggthemes::excel_Headlines")[c(2,5,3,1)],
                     #lty = c(1,2,1,2),
                     legend.title="RASH", 
                     legend.labs=c("G/E:Basal-like","G/E:Classical","FFX:Basal-like","FFX:Classical"), # legend
                     xlab = "", ylab="Survival probability", xlim = c(0,30),
                     pval = T, pval.size = 9, pval.coord = c(22,0.75),  # p-values
                     font.main = c(23, "plain", "black"), 
                     font.legend = c(19, "plain"),
                     font.x = c(23, "plain", "black"),
                     font.y = c(23, "plain", "black"),
                     font.tickslab = c(23, "plain", "black"), ) + 
  guides(colour = guide_legend(nrow = 2)) 
ggsurv$table <- ggrisktable(km_fit, 
                            data = survDatObj, 
                            break.time.by = 6,
                            color = "strata", 
                            palette = paletteer::paletteer_d("ggthemes::excel_Headlines")[c(2,5,3,1)],
                            y.text = F,  
                            fontsize = 8,
                            xlab = "Time (months)", ylab = "", xlim = c(0,30),
                            tables.theme = theme_survminer(font.main = 23, font.tickslab = 23, font.x = 23))
print(ggsurv)



# Suppl Fig 6d&e  ---------------------------------------------------------------------------------------------------------------------------------------------
# boxplot

clinicalDat <- dataSet$clinicalDat_Parsed

datBox <- data.frame(PurIST = clinicalDat$PurIST.prob,
                     Arm = clinicalDat$arm_all,
                     Study = clinicalDat$study,
                     Rash = clinicalDat$rash_grade,
                     stringsAsFactors = FALSE)

datBoxTmp <- datBox[which(datBox$Study %in% c("Rash")), ]
datBoxTmp$Rash[which(datBoxTmp$Rash != 0)] <- 1
datBoxTmp$Rash <- factor(datBoxTmp$Rash, levels = c(1,0))
levels(datBoxTmp$Rash) <- c("Rash","No rash")
pVal <- round(wilcox.test(PurIST~Rash, data = datBoxTmp)$p.value,3)
boxplot(PurIST~Rash, data = datBoxTmp,
        outline = F, #las=2,
        border = c("salmon","grey"),
        lwd =1, col = "white",
        ylim = c(0,1), 
        xlab = "",ylab = "PurIST\nbasal-like probability",
        cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.1,
        main=paste("RASH trial\np=",pVal,sep=""))

# RASH
datBoxTmp <- datBox[which(datBox$Arm %in% c("Gemcitabine/Erlotinib","FOLFIRINOX")), ]
datBoxTmp$Arm <- factor(datBoxTmp$Arm)
pVal <- round(wilcox.test(PurIST~Arm, data = datBoxTmp)$p.value,3)
boxplot(PurIST~Arm, data = datBoxTmp,
        outline = F, #las=2,
        border = armColList[c(4,5)],
        lwd =1, col = "white",
        ylim = c(0,1), names = c("G/E","FFX"),
        xlab = "",ylab = "PurIST\nbasal-like probability",
        cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.1,
        main=paste("RASH trial\np=",pVal,sep=""))

# Suppl Fig. 6c  ----------------------------------------------------------------------------------------
# barplot

clinicalDat <- dataSet$clinicalDat_Parsed

# RASH
datTmp <- clinicalDat[which(clinicalDat$study == "Rash" ),
                  c("arm_all","resp3","PurIST")]
datTmp$arm_all <- droplevels(datTmp$arm_all)

datTmp <- clinicalDat[which(clinicalDat$arm_all %in% c("Gemcitabine/Erlotinib","FOLFIRINOX") ),
                      c("arm_all","resp3","PurIST")]
datTmp$arm_all <- droplevels(datTmp$arm_all)
tableTmp <- as.data.frame(table(datTmp))

# ACCEPT
datTmp <- clinicalDat[which(clinicalDat$study == "Accept" ),
                      c("arm_all","resp3","PurIST")]
datTmp$arm_all <- droplevels(datTmp$arm_all)
tableTmp <- as.data.frame(table(datTmp))

dev.off()


