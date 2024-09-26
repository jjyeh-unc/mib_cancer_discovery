Relabel_subtype <- function(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale) {
  
  tmpRowScaleLab <- c("unscaled", "scaled")[(tmpRowScale + 1)]
  #tmpMainLab <- paste(tmpSchema,".",tmpRowScaleLab,".","K",tmpK,sep="")
  tmpMainLab <- paste(tmpSchema,".",tmpRowScaleLab,sep="")
  
  tmp <- dataSet[[tmpMainLab]]
  sampSub <- tmp$sampSub
  
  # find schema #
  i <- which(schemaList %in% tmp$schema)
  
  # assign subtype label
  tmp$subtypeCall.cluster <- tmp$subtypeCall
  tmp$subtypeCall <- tmpSubtypeLab[tmp$subtypeCall]
  
  # save subtype to sampInfo
  dataSet$sampInfo$tmpCluster <- NA
  dataSet$sampInfo$tmpCluster[sampSub] <- tmp$subtypeCall
  names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- tmp$mainLab
  
  # color samples
  ColSideColors <- c(subtypeColList[[i]])[match(tmp$subtypeCall, subtypeList[[i]])]
  ColSideColors[which(!(ColSideColors %in% c(subtypeColList[[i]])))] <- "grey"
  ColSideColors[which(tmp$subtypeCall %in% "Absent")] <- "black"
  ColSideColors[which(tmp$subtypeCall %in% "Mixed")] <- "grey"
  tmp$ColSideColors <- ColSideColors
  
  # plot heatmap
  Plot_heatmap_CC(tmp, tmpSchema, ColSideColors, 
                  paste(tmp$mainLab,".","K",tmpK,sep="") )
  
  # assign tmp
  dataSet[[tmpMainLab]] <- tmp
  
  return(dataSet)
}