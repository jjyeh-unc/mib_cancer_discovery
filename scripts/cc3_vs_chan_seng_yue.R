# set dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

############################################################
#################### Functions and libraries

# load libraries
library(stringr)
library(openxlsx)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)

# load functions
file.sources <- list.files("./R/",pattern="*.R")
file.sources <- paste("./R/", file.sources, sep="")
sapply(file.sources, source)

# load subtype info
load("cmbSubtypes.230823.RData")

# function
Plot_heatmap_CC <- function(tmp, schemaTmp, ColSideColors, mainLab=schemaTmp) {
  i <- which(schemaList %in% schemaTmp)
  genesTmp <- subtypeGeneList[[i]]
  idxGene <- match(tmp$featInfo,genesTmp$geneSymbol)
  geneInfo <- genesTmp[idxGene,]
  
  if(ncol(as.matrix(ColSideColors))<=15) {
    ColSideColorsSize <- ncol(as.matrix(ColSideColors))
  } else {
    ColSideColorsSize <- 15
  }
  
  heatmap.3(main= mainLab
            ,x = as.matrix(tmp$ex)
            ,Rowv = tmp$Rowv
            ,Colv = tmp$Colv
            ,dendrogram = "both" ,trace="none" ,scale="row"
            ,labRow = tmp$featInfo
            ,labCol = tmp$sampID
            ,col = colorRampPalette(c("darkblue","blue","white","red","darkred"))(99)
            ,ColSideColors = as.matrix(ColSideColors)
            ,ColSideColorsSize = ColSideColorsSize
            ,RowSideColors = t(as.matrix(geneInfo$Color))
            ,RowSideColorsSize=1
            ,lwid = c(1,5)
            ,lhei = c(1,3)
            ,margins =c(15,18)
            ,keysize=0.5
            ,cex.lab=0.4
            ,cexRow = 0.5
            ,cexCol = 0.6
  )
}
############################## Analyze ################################

sampInfo <- read.csv("CC3.csv")
tmpSplit <- data.frame(str_split_fixed(sampInfo$X, "_", 2))
sampInfo$Line <-tmpSplit$X1
newClusterCall <- data.frame(Line = unique(sampInfo$Line),
                             Cluster = "",
                             stringsAsFactors = FALSE)
for (i in 1:nrow(newClusterCall)) {
  tmpDat <- sampInfo[which(sampInfo$Line %in% newClusterCall$Line[i]),]
  newClusterCall$Cluster[i] <- names(which.max(table(tmpDat$CC2.c3.consensusClass)))
}
newClusterCall <- newClusterCall[order(newClusterCall$Cluster),]


# load RNA data ---------------------------------------------------------------------------------
datRNA <- read.csv("MIB_macthed_RNAseq.csv")
rownames(datRNA) <- datRNA$X
datRNA <- datRNA[,-1]
dataSet <- list()
sampInfo <- data.frame(sampleID = colnames(datRNA), stringsAsFactors = FALSE)
tmpSplit <- data.frame(str_split_fixed(sampInfo$sampleID, "_", 2))
sampInfo$Line <- tmpSplit$X1
sampInfo$Rep <- tmpSplit$X2
sampInfo$KinaseCluster <- newClusterCall$Cluster[match(sampInfo$Line, newClusterCall$Line)]
dataSet$sampInfo <- sampInfo
dataSet$ex <- datRNA
dataSet$featInfo <- data.frame(SYMBOL = rownames(datRNA), stringsAsFactors = FALSE)
dataSet$metadata$log.transformed <- TRUE

# purist
#dataSet <- Call_PurIST(dataSet)

# cluster 
sampSub <- 1:nrow(dataSet$sampInfo)
## Chan-Seng-Yue
tmpSchema <- "Chan-Seng-Yue"
tmpK <- 5
tmpRowScale <- FALSE
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 4, tmpK, tmpRowScale, TRUE, FALSE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("Classical B","Basal-like B","Hybrid","Classical A","Basal-like A")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, 5, tmpRowScale)
table(dataSet$sampInfo[,c("Chan-Seng-Yue.unscaled","KinaseCluster")])
write.xlsx(dataSet$sampInfo, "EGFR_paper_k3_sampInfo.xlsx")

pdf("Chan-Seng-Yue_vs_CC3.heatmap.pdf")

# heatmap
lineCol <- c("P319T1" = "#1F77B4",
  "P710T1" = "#FF7F0E",
  "P203T1" = "#2CA02C",
  "P616T1" = "#D62728",
  "P109T1" = "#9467BD",
  "P411T1" = "#8C564B",
  "P125T2" = "#E377C2",
  "P119T1" = "#7F7F7F",
  "P225T1" = "#BCBD22",
  "PancT6" = "#17BECF",
  "P508T1" = "#AEC7E8",
  "P902T1B" = "#FFBB78",
  "P806T1" = "#98DF8A",
  "P713T1" = "#FF9896")

# with PurIST
ColSideColors <- data.frame(#Line = c(brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 8, name = "Dark2")[c(4,7,8)])[factor(dataSet$sampInfo$Line)],  
                            Line = lineCol[match(dataSet$sampInfo$Line, names(lineCol))],
                            KinaseCluster = c("#B15928","#6A3D9A","#f4cc03")[as.numeric(dataSet$sampInfo$KinaseCluster)],
                            #PurIST = c("orange","blue")[as.factor(dataSet$sampInfo$PurIST)], 
                            Chan_Seng_Yue = Get_subtype_color(dataSet$sampInfo$`Chan-Seng-Yue.unscaled`,"Chan-Seng-Yue"),
                            stringsAsFactors = FALSE)
Plot_heatmap_CC(dataSet$`Chan-Seng-Yue.unscaled`, "Chan-Seng-Yue", ColSideColors, "Kinome-matched RNAseq samples")
legend(xy.coords(x=.75,y=.98),
       legend=c("Kinase Cluster", "1","2","3",
                "Line",names(lineCol)),
       fill=c("white",c("#B15928","#6A3D9A","#f4cc03"),
              "white","brown","orange","grey","blue3","skyblue",
              "white",lineCol),
       border=FALSE, bty="n",
       y.intersp = 0.85 , cex = 0.8)

dev.off()