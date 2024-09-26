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
library(tidyverse)
#library(tidytidbits)
library(survivalAnalysis)
library(stringr)
library(patchwork)
library(forestplot)

# load functions
file.sources <- list.files("./R/",pattern="*.R")
file.sources <- paste("./R/", file.sources, sep="")
sapply(file.sources, source)

# load subtype info
### This is a combined subtype object with
### subtypeColList, subtypeGeneList, subtypeList and schemaList
#load("cmbSubtypes.220922.RData")
#print("Subtype schemas available for use:")
#print(schemaList)

############################## Parsing ##############################
# load data
dataSet <- readRDS("RASH_ACCEPT.rds")
survDat <- dataSet$clinicalDat_Parsed
survDat.okane <- readRDS("OKane.Salmon_Gencode.rds")
survDat.okane0 <- survDat.okane

# Fig. 6h  ---------------------------------------------------------------------------------------------------------------------------------------------
# forest plot

# initialize
statsHR <- data.frame(matrix(ncol = 5))
names(statsHR) <- c("Group","HR","Lower 95% CI","Upper 95% CI","P-value")

# ACCEPT -----------------------------------------

# Basal-like
survDatTmp0 <- survDat[which(survDat$study=="Accept" & survDat$PurIST == "Basal-like"),c("arm_all","os","os_time")]
survDatObj <- survDatTmp0 %>%
  mutate(arm_all = case_when(
    arm_all == "Gemcitabine" ~ "G",
    arm_all == "Gemcitabine/Afatinib" ~ "G/A",
    TRUE ~ arm_all  # Keep other values unchanged
  ))
survDatObj$arm_all <- factor(survDatObj$arm_all,
                             levels = c("G","G/A"))

fit <- coxph(Surv(os_time, os) ~ arm_all, data = survDatObj)
summary_fit <- summary(fit)
hazard_ratios <- summary_fit$coefficients[, "exp(coef)"]
conf_int <- summary_fit$conf.int
lower_bound <- conf_int[, "lower .95"]
upper_bound <- conf_int[, "upper .95"]
p_value <- summary_fit$coefficients[, 'Pr(>|z|)']
statsTmp <- c("ACCEPT basal-like G/A vs G",
              hazard_ratios,lower_bound,upper_bound,p_value)
statsHR[1,] <- statsTmp

# Classical
survDatTmp0 <- survDat[which(survDat$study=="Accept" & survDat$PurIST == "Classical"),c("arm_all","os","os_time")]
survDatObj <- survDatTmp0 %>%
  mutate(arm_all = case_when(
    arm_all == "Gemcitabine" ~ "G",
    arm_all == "Gemcitabine/Afatinib" ~ "G/A",
    TRUE ~ arm_all  # Keep other values unchanged
  ))
survDatObj$arm_all <- factor(survDatObj$arm_all,
                             levels = c("G","G/A"))

fit <- coxph(Surv(os_time, os) ~ arm_all, data = survDatObj)
summary_fit <- summary(fit)
hazard_ratios <- summary_fit$coefficients[, "exp(coef)"]
conf_int <- summary_fit$conf.int
lower_bound <- conf_int[, "lower .95"]
upper_bound <- conf_int[, "upper .95"]
p_value <- summary_fit$coefficients[, 'Pr(>|z|)']
statsTmp <- c("ACCEPT classical G/A vs G",
              hazard_ratios,lower_bound,upper_bound,p_value)
statsHR[2,] <- statsTmp

# RASH -----------------------------------------

# Basal-like
survDatTmp0 <- survDat[which(survDat$study=="Rash" & survDat$PurIST == "Basal-like"),c("arm_all","os","os_time")]
survDatObj <- survDatTmp0 %>%
  mutate(arm_all = case_when(
    arm_all %in% c("FOLFIRINOX") ~ "G/E_to_FFX",
    arm_all == "Gemcitabine/Erlotinib" ~ "G/E",
    TRUE ~ arm_all  # Keep other values unchanged
  )) %>%
  filter(arm_all != "No_further")
survDatObj$arm_all <- factor(survDatObj$arm_all,
                             levels = c("G/E_to_FFX","G/E"))

fit <- coxph(Surv(os_time, os) ~ arm_all, data = survDatObj)
summary_fit <- summary(fit)
hazard_ratios <- summary_fit$coefficients[, "exp(coef)"]
conf_int <- summary_fit$conf.int
lower_bound <- conf_int[, "lower .95"]
upper_bound <- conf_int[, "upper .95"]
p_value <- summary_fit$coefficients[, 'Pr(>|z|)']
statsTmp <- c("RASH basal-like G/E vs G/E_to_FFX",
              hazard_ratios,lower_bound,upper_bound,p_value)
statsHR[3,] <- statsTmp

# Classical
survDatTmp0 <- survDat[which(survDat$study=="Rash" & survDat$PurIST == "Classical"),c("arm_all","os","os_time")]
survDatObj <- survDatTmp0 %>%
  mutate(arm_all = case_when(
    arm_all %in% c("FOLFIRINOX") ~ "G/E_to_FFX",
    arm_all == "Gemcitabine/Erlotinib" ~ "G/E",
    TRUE ~ arm_all  # Keep other values unchanged
  )) %>%
  filter(arm_all != "No_further")
survDatObj$arm_all <- factor(survDatObj$arm_all,
                             levels = c("G/E_to_FFX","G/E"))

fit <- coxph(Surv(os_time, os) ~ arm_all, data = survDatObj)
summary_fit <- summary(fit)
hazard_ratios <- summary_fit$coefficients[, "exp(coef)"]
conf_int <- summary_fit$conf.int
lower_bound <- conf_int[, "lower .95"]
upper_bound <- conf_int[, "upper .95"]
p_value <- summary_fit$coefficients[, 'Pr(>|z|)']
statsTmp <- c("RASH classical G/E vs G/E_to_FFX",
              hazard_ratios,lower_bound,upper_bound,p_value)
statsHR[4,] <- statsTmp


# COMPASS -------------------------------------------

# Basal-like
survDat.okane <- survDat.okane0
survDat.okane <- survDat.okane$sampInfo[which( (survDat.okane$sampInfo$PurIST %in% "Basal-like") &
                                                 !(survDat.okane$sampInfo$Study.ID %in% "COMP-0021") &
                                                 !(survDat.okane$sampInfo$First.Tx %in% c("none","GA/experimental","Gemcitabine"))), ]
survDat.okane$First.Tx[which(survDat.okane$First.Tx %in% "GA")] <- "GnP"
survDat.okane <- data.frame(arm_all = paste("COMPASS",survDat.okane$First.Tx,sep = "_"),
                            os = survDat.okane$Alive,
                            os_time = as.numeric(survDat.okane$Survival.days/30),
                            stringsAsFactors = FALSE)
survDat.okane$os[which(survDat.okane$os %in% "Dead")] <- 1
survDat.okane$os[which(survDat.okane$os %in% "Alive")] <- 0
survDat.okane$os <- as.numeric(survDat.okane$os)
survDat.okane$arm_all <- factor(survDat.okane$arm_all,
                                levels = c("COMPASS_FFX","COMPASS_GnP"))
survDatObj <- survDat.okane

fit <- coxph(Surv(os_time, os) ~ arm_all, data = survDatObj)
summary_fit <- summary(fit)
hazard_ratios <- summary_fit$coefficients[, "exp(coef)"]
conf_int <- summary_fit$conf.int
lower_bound <- conf_int[, "lower .95"]
upper_bound <- conf_int[, "upper .95"]
p_value <- summary_fit$coefficients[, 'Pr(>|z|)']
statsTmp <- c("COMPASS basal-like GnP vs FFX",
              hazard_ratios[1],lower_bound[1],upper_bound[1],p_value[1])
statsHR[5,] <- statsTmp

# Classical
survDat.okane <- survDat.okane0
survDat.okane <- survDat.okane$sampInfo[which( (survDat.okane$sampInfo$PurIST %in% "Classical") &
                                                 !(survDat.okane$sampInfo$Study.ID %in% "COMP-0021") &
                                                 !(is.na(survDat.okane$sampInfo$First.Tx)) &
                                                 !(survDat.okane$sampInfo$First.Tx %in% c("none","GA/experimental","cisplatin/gemcitabine","Gemcitabine"))), ]
survDat.okane$First.Tx[which(survDat.okane$First.Tx %in% "GA")] <- "GnP"
survDat.okane <- data.frame(arm_all = paste("COMPASS",survDat.okane$First.Tx,sep = "_"),
                            os = survDat.okane$Alive,
                            os_time = as.numeric(survDat.okane$Survival.days/30),
                            stringsAsFactors = FALSE)
survDat.okane$os[which(survDat.okane$os %in% "Dead")] <- 1
survDat.okane$os[which(survDat.okane$os %in% "Alive")] <- 0
survDat.okane$os <- as.numeric(survDat.okane$os)
survDat.okane$arm_all <- factor(survDat.okane$arm_all,
                                levels = c("COMPASS_FFX","COMPASS_GnP"))
survDatObj <- survDat.okane

fit <- coxph(Surv(os_time, os) ~ arm_all, data = survDatObj)
summary_fit <- summary(fit)
hazard_ratios <- summary_fit$coefficients[, "exp(coef)"]
conf_int <- summary_fit$conf.int
lower_bound <- conf_int[, "lower .95"]
upper_bound <- conf_int[, "upper .95"]
p_value <- summary_fit$coefficients[, 'Pr(>|z|)']
statsTmp <- c("COMPASS classical GnP vs FFX",
              hazard_ratios[1],lower_bound[1],upper_bound[1],p_value[1])
statsHR[6,] <- statsTmp

# forest plot
splitTmp <- data.frame(str_split_fixed(statsHR$Group, " ", 3))
names(splitTmp) <- c("Study","Subtype","Arms")
statsPlot <- data.frame(Group = statsHR$Group,
                        stringsAsFactors = FALSE)
statsPlot$HR <- log2(as.numeric(statsHR$HR))
statsPlot$CI_lower <- log2(as.numeric(statsHR$`Lower 95% CI`))
statsPlot$CI_upper <- log2(as.numeric(statsHR$`Upper 95% CI`))
statsPlot$Pvalue <- round(as.numeric(statsHR$`P-value`),3)
statsPlot$Subtype <- splitTmp$Subtype
statsPlot$HR <- as.numeric(statsPlot$HR)
statsPlot$CI_lower <- as.numeric(statsPlot$CI_lower)
statsPlot$CI_upper <- as.numeric(statsPlot$CI_upper)

# ggplot
statsPlot$Group <- factor(statsPlot$Group, levels = statsPlot$Group[order(statsPlot$HR, decreasing = T)])
ggplot(statsPlot, aes(x = HR, y = Group, color = Subtype)) +
  scale_color_manual(values = c("basal-like" = "orange", "classical" = "blue")) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
  theme_minimal() +
  labs(x = "Log2 HR", y = "Treatment comprisons by study and subtype", title = "") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),  # Larger x-axis label size
    axis.title.y = element_text(size = 16, face = "bold"),  # Larger y-axis label size
    axis.text.x = element_text(size = 14),  # Larger x-axis text size
    axis.text.y = element_text(size = 14)   # Larger y-axis text size
  )

# Suppl Fig f,g,i,j  ---------------------------------------------------------------------------------------------------------------------------------------------
# survival plot 

survDat <- dataSet$clinicalDat_Parsed

# RASH -----------------------------------------

## basal + classical
survDatTmp0 <- survDat[which((survDat$arm_all %in% c("Gemcitabine/Erlotinib") )),c("arm_all","os","os_time")]
survDat.okane <- survDat.okane0
survDat.okane <- survDat.okane$sampInfo[which(  !(survDat.okane$sampInfo$Study.ID %in% "COMP-0021") &
                                                  !(survDat.okane$sampInfo$First.Tx %in% c("none","GA/experimental","cisplatin/gemcitabine"))), ]
survDat.okane$First.Tx[which(survDat.okane$First.Tx %in% "GA")] <- "GnP"
survDat.okane <- data.frame(arm_all = paste("COMPASS",survDat.okane$First.Tx,sep = "_"),
                            #First.Tx = survDat.okane$First.Tx,
                            os = survDat.okane$Alive,
                            os_time = as.numeric(survDat.okane$Survival.days/30),
                            stringsAsFactors = FALSE)
survDat.okane$os[which(survDat.okane$os %in% "Dead")] <- 1
survDat.okane$os[which(survDat.okane$os %in% "Alive")] <- 0
survDatTmp <- rbind(survDatTmp0, survDat.okane)

survDatObj <- survDatTmp
survDatObj$os <- as.numeric(survDatObj$os)
survDatObj$arm_all <- factor(survDatObj$arm_all,
                             levels = c("Gemcitabine/Erlotinib", "COMPASS_FFX", "COMPASS_GnP", "COMPASS_Gemcitabine"))
km_fit <- survfit(Surv(survDatObj[,"os_time"], survDatObj[, "os"]) ~ arm_all, data = survDatObj, type = "kaplan-meier")
p <- ggsurvplot(km_fit, 
                data = survDatObj,
                conf.int = F, 
                legend.title="",
                legend = c(0.75, 0.75),
                font.legend = 16, 
                font.title = 16,        
                font.x = 20,  
                font.y = 20, 
                font.tickslab = 20,
                risk.table = T,
                risk.table.pos = "out",
                risk.table.fontsize = 5,
                risk.table.height = 0.3,
                surv.median.line = "hv", 
                size = 0.5, 
                censor.size=3, 
                pval = T,
                pval.size = 8,
                pval.coord = c(24, 0.4),
                xlim = c(0,40),
                xlab = "Time (months)", 
                break.time.by = 12,
                palette = c("#e36ca3","cyan3","gold3","yellowgreen"),
                linetype = c("solid","dashed","dashed","dashed"),
                legend.labs=c("RASH_G/E","COMPASS_mFFX","COMPASS_GnP","COMPASS_G"),
                title = "RASH vs COMPASS: all")
p

## Basal-like
survDatTmp0 <- survDat[which((survDat$arm_all %in% c("Gemcitabine/Erlotinib")) & survDat$PurIST %in% "Basal-like" ),c("arm_all","os","os_time")]
survDat.okane <- survDat.okane0
survDat.okane <- survDat.okane$sampInfo[which( (survDat.okane$sampInfo$PurIST %in% "Basal-like") &
                                                 !(survDat.okane$sampInfo$Study.ID %in% "COMP-0021") &
                                                 !(survDat.okane$sampInfo$First.Tx %in% c("none","GA/experimental"))), ]
survDat.okane$First.Tx[which(survDat.okane$First.Tx %in% "GA")] <- "GnP"
survDat.okane <- data.frame(arm_all = paste("COMPASS",survDat.okane$First.Tx,sep = "_"),
                            #First.Tx = survDat.okane$First.Tx,
                            os = survDat.okane$Alive,
                            os_time = as.numeric(survDat.okane$Survival.days/30),
                            stringsAsFactors = FALSE)
survDat.okane$os[which(survDat.okane$os %in% "Dead")] <- 1
survDat.okane$os[which(survDat.okane$os %in% "Alive")] <- 0
survDatTmp <- rbind(survDatTmp0, survDat.okane)

survDatObj <- survDatTmp
survDatObj$os <- as.numeric(survDatObj$os)
survDatObj$arm_all <- factor(survDatObj$arm_all,
                             levels = c("Gemcitabine/Erlotinib", "COMPASS_FFX", "COMPASS_GnP", "COMPASS_Gemcitabine"))
km_fit <- survfit(Surv(survDatObj[,"os_time"], survDatObj[, "os"]) ~ arm_all, data = survDatObj, type = "kaplan-meier")
p <- ggsurvplot(km_fit, 
                data = survDatObj,
                conf.int = F, 
                legend.title="",
                legend = c(0.75, 0.75),
                font.legend = 16, 
                font.title = 16,         
                font.x = 20,  
                font.y = 20, 
                font.tickslab = 20,
                risk.table = T,
                risk.table.pos = "out",
                risk.table.fontsize = 5,
                risk.table.height = 0.3,
                surv.median.line = "hv", 
                size = 0.5, 
                censor.size=3, 
                pval = T,
                pval.size = 8,
                pval.coord = c(24, 0.4),
                xlim = c(0,40),
                xlab = "Time (months)", 
                break.time.by = 12,
                palette = c("#e36ca3","cyan3","gold3","yellowgreen"),
                linetype = c("solid","dashed","dashed","dashed"),
                legend.labs=c("RASH_G/E","COMPASS_mFFX","COMPASS_GnP","COMPASS_G"),
                title = "RASH vs COMPASS: basal-like") 
p

## Classical
survDatTmp0 <- survDat[which((survDat$arm_all %in% c("Gemcitabine/Erlotinib") ) & survDat$PurIST == "Classical"),c("arm_all","os","os_time")]
survDat.okane <- survDat.okane0
survDat.okane <- survDat.okane$sampInfo[which( (survDat.okane$sampInfo$PurIST %in% "Classical") &
                                                 !(survDat.okane$sampInfo$Study.ID %in% "COMP-0021") &
                                                 !(survDat.okane$sampInfo$First.Tx %in% c("none","GA/experimental","cisplatin/gemcitabine"))), ]
survDat.okane$First.Tx[which(survDat.okane$First.Tx %in% "GA")] <- "GnP"
survDat.okane <- data.frame(arm_all = paste("COMPASS",survDat.okane$First.Tx,sep = "_"),
                            os = survDat.okane$Alive,
                            os_time = as.numeric(survDat.okane$Survival.days/30),
                            stringsAsFactors = FALSE)
survDat.okane$os[which(survDat.okane$os %in% "Dead")] <- 1
survDat.okane$os[which(survDat.okane$os %in% "Alive")] <- 0
survDatTmp <- rbind(survDatTmp0, survDat.okane)

survDatObj <- survDatTmp
survDatObj$os <- as.numeric(survDatObj$os)
survDatObj$arm_all <- factor(survDatObj$arm_all,
                             levels = c("Gemcitabine/Erlotinib", "COMPASS_FFX", "COMPASS_GnP", "COMPASS_Gemcitabine"))
km_fit <- survfit(Surv(survDatObj[,"os_time"], survDatObj[, "os"]) ~ arm_all, data = survDatObj, type = "kaplan-meier")
p <- ggsurvplot(km_fit, 
                data = survDatObj,
                conf.int = F, 
                legend.title="",
                legend = c(0.75, 0.75),
                font.legend = 16, 
                font.title = 16,           
                font.x = 20,  
                font.y = 20, 
                font.tickslab = 20,
                risk.table = T,
                risk.table.pos = "out",
                risk.table.fontsize = 5,
                risk.table.height = 0.3,
                surv.median.line = "hv", 
                size = 0.5, 
                censor.size=3, 
                pval = T,
                pval.size = 8,
                pval.coord = c(24, 0.4),
                xlim = c(0,40),
                xlab = "Time (months)", 
                break.time.by = 12,
                palette = c("#e36ca3","cyan3","gold3","yellowgreen"),
                linetype = c("solid","dashed","dashed","dashed"),
                legend.labs=c("RASH_G/E","COMPASS_mFFX","COMPASS_GnP","COMPASS_G"),
                title = "RASH + OKane: classical") 

p

# Accept -----------------------------------------

## all
survDatTmp0 <- survDat[which(survDat$study=="Accept" ),c("arm_all","os","os_time")]

survDat.okane <- survDat.okane0
survDat.okane <- survDat.okane$sampInfo[which(!(survDat.okane$sampInfo$Study.ID %in% "COMP-0021") &
                                                !is.na(survDat.okane$sampInfo$First.Tx) &
                                                !(survDat.okane$sampInfo$First.Tx %in% c("none","GA/experimental","cisplatin/gemcitabine"))), ]
survDat.okane$First.Tx[which(survDat.okane$First.Tx %in% "GA")] <- "GnP"
survDat.okane <- data.frame(arm_all = paste("COMPASS",survDat.okane$First.Tx,sep = "_"),
                            os = survDat.okane$Alive,
                            os_time = as.numeric(survDat.okane$Survival.days/30),
                            stringsAsFactors = FALSE)
survDat.okane$os[which(survDat.okane$os %in% "Dead")] <- 1
survDat.okane$os[which(survDat.okane$os %in% "Alive")] <- 0
survDatTmp <- rbind(survDatTmp0, survDat.okane)

survDatObj <- survDatTmp
survDatObj$os <- as.numeric(survDatObj$os)
survDatObj$arm_all <- factor(survDatObj$arm_all,
                             levels = c("Gemcitabine/Afatinib", "Gemcitabine", "COMPASS_FFX", "COMPASS_GnP", "COMPASS_Gemcitabine"))
km_fit <- survfit(Surv(survDatObj[,"os_time"], survDatObj[, "os"]) ~ arm_all, data = survDatObj, type = "kaplan-meier")
p <- ggsurvplot(km_fit, 
                data = survDatObj,
                conf.int = F, 
                legend.title="",
                legend = c(0.75, 0.75),
                font.legend = 16, 
                font.title = 16,          
                font.x = 20,  
                font.y = 20, 
                font.tickslab = 20,
                risk.table = T,
                risk.table.pos = "out",
                risk.table.fontsize = 5,
                risk.table.height = 0.3,
                surv.median.line = "hv", 
                size = 0.5, 
                censor.size=3, 
                pval = T,
                pval.size = 8,
                pval.coord = c(24, 0.4),
                xlim = c(0,40),
                xlab = "Time (months)", 
                break.time.by = 12,
                palette = c("#9f79ee","#9acd32","cyan3","gold3","yellowgreen"),
                linetype = c("solid","solid","dashed","dashed","dashed"),
                legend.labs=c("ACCEPT_G/A","ACCEPT_G","COMPASS_mFFX","COMPASS_GnP","COMPASS_G"),
                title = "ACCEPT vs COMPASS: all") 

p

## Basal-like ----------------------------------------
survDatTmp0 <- survDat[which(survDat$study=="Accept" & survDat$PurIST == "Basal-like"),c("arm_all","os","os_time")]
survDat.okane <- survDat.okane0
survDat.okane <- survDat.okane$sampInfo[which( (survDat.okane$sampInfo$PurIST %in% "Basal-like") &
                                                 !(survDat.okane$sampInfo$Study.ID %in% "COMP-0021") &
                                                 !is.na(survDat.okane$sampInfo$First.Tx) &
                                                 !(survDat.okane$sampInfo$First.Tx %in% c("none","GA/experimental","cisplatin/gemcitabine"))), ]
survDat.okane$First.Tx[which(survDat.okane$First.Tx %in% "GA")] <- "GnP"
survDat.okane <- data.frame(arm_all = paste("COMPASS",survDat.okane$First.Tx,sep = "_"),
                            os = survDat.okane$Alive,
                            os_time = as.numeric(survDat.okane$Survival.days/30),
                            stringsAsFactors = FALSE)
survDat.okane$os[which(survDat.okane$os %in% "Dead")] <- 1
survDat.okane$os[which(survDat.okane$os %in% "Alive")] <- 0
survDatTmp <- rbind(survDatTmp0, survDat.okane)

survDatObj <- survDatTmp
survDatObj$os <- as.numeric(survDatObj$os)
survDatObj$arm_all <- factor(survDatObj$arm_all,
                             levels = c("Gemcitabine/Afatinib", "Gemcitabine", "COMPASS_FFX", "COMPASS_GnP", "COMPASS_Gemcitabine"))

km_fit <- survfit(Surv(survDatObj[,"os_time"], survDatObj[, "os"]) ~ arm_all, data = survDatObj, type = "kaplan-meier")
p <- ggsurvplot(km_fit, 
                data = survDatObj,
                conf.int = F, 
                legend.title="",
                legend = c(0.75, 0.75),
                font.legend = 16, 
                font.title = 16,              
                font.x = 20,  
                font.y = 20, 
                font.tickslab = 20,
                risk.table = T,
                risk.table.pos = "out",
                risk.table.fontsize = 5,
                risk.table.height = 0.3,
                surv.median.line = "hv", 
                size = 0.5, 
                censor.size=3, 
                pval = T,
                pval.size = 8,
                pval.coord = c(24, 0.4),
                xlim = c(0,40),
                xlab = "Time (months)", 
                break.time.by = 12,
                palette = c("#9f79ee","#9acd32","cyan3","gold3","yellowgreen"),
                linetype = c("solid","solid","dashed","dashed","dashed"),
                legend.labs=c("ACCEPT_G/A","ACCEPT_G","COMPASS_mFFX","COMPASS_GnP","COMPASS_G"),
                title = "ACCEPT vs COMPASS: basal-like") 
p


## Classical
survDatTmp0 <- survDat[which(survDat$study=="Accept" & survDat$PurIST == "Classical"),c("arm_all","os","os_time")]

survDat.okane <- survDat.okane0
survDat.okane <- survDat.okane$sampInfo[which( (survDat.okane$sampInfo$PurIST %in% "Classical") &
                                                 !(survDat.okane$sampInfo$Study.ID %in% "COMP-0021") &
                                                 !is.na(survDat.okane$sampInfo$First.Tx) &
                                                 !(survDat.okane$sampInfo$First.Tx %in% c("none","GA/experimental","cisplatin/gemcitabine"))), ]
survDat.okane$First.Tx[which(survDat.okane$First.Tx %in% "GA")] <- "GnP"
survDat.okane <- data.frame(arm_all = paste("COMPASS",survDat.okane$First.Tx,sep = "_"),
                            os = survDat.okane$Alive,
                            os_time = as.numeric(survDat.okane$Survival.days/30),
                            stringsAsFactors = FALSE)
survDat.okane$os[which(survDat.okane$os %in% "Dead")] <- 1
survDat.okane$os[which(survDat.okane$os %in% "Alive")] <- 0
survDatTmp <- rbind(survDatTmp0, survDat.okane)

survDatObj <- survDatTmp
survDatObj$os <- as.numeric(survDatObj$os)
survDatObj$arm_all <- factor(survDatObj$arm_all,
                             levels = c("Gemcitabine/Afatinib", "Gemcitabine", "COMPASS_FFX", "COMPASS_GnP", "COMPASS_Gemcitabine"))
km_fit <- survfit(Surv(survDatObj[,"os_time"], survDatObj[, "os"]) ~ arm_all, data = survDatObj, type = "kaplan-meier")
p <- ggsurvplot(km_fit, 
                data = survDatObj,
                conf.int = F, 
                legend.title="",
                legend = c(0.75, 0.75),
                font.legend = 16, 
                font.title = 16,             
                font.x = 20,  
                font.y = 20, 
                font.tickslab = 20,
                risk.table = T,
                risk.table.pos = "out",
                risk.table.fontsize = 5,
                risk.table.height = 0.3,
                surv.median.line = "hv", 
                size = 0.5, 
                censor.size=3, 
                pval = T,
                pval.size = 8,
                pval.coord = c(24, 0.4),
                xlim = c(0,40),
                xlab = "Time (months)", 
                break.time.by = 12,
                palette = c("#9f79ee","#9acd32","cyan3","gold3","yellowgreen"),
                linetype = c("solid","solid","dashed","dashed","dashed"),
                legend.labs=c("ACCEPT_G/A","ACCEPT_G","COMPASS_mFFX","COMPASS_GnP","COMPASS_G"),
                title = "ACCEPT vs COMPASS: classical") 
#p$table <- p$table + layer_scales(p$plot)$x
p

