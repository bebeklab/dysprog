#!/usr/bin/Rscript
suppressMessages(library(ggplot2))
suppressMessages(library(NMF))
suppressMessages(library(RColorBrewer))
suppressMessages(library(GGally))
suppressMessages(library(survival))

## NMF Clustering ##
args = commandArgs(trailingOnly = TRUE)

if(length(args) < 3){
  cat(sprintf(" !Missing Files! Usage: <FSM Score File> <Clinical File> <Output Prefix> <Event Days Threshold (Optional, Default 10000)>\n"))
  quit(save = "no", status = 1)
}else{
  score.file <- args[1]
  clinical.file <- args[2]
  out.file <- args[3]
  thr <- ifelse(!is.na(args[4]), args[4], 3500)
}

log <- file(sprintf("%s-log.txt", out.file), open = "a")
cat(sprintf("\n-- NMF Clustering (%s) --\n", Sys.time()), file = log)
cat(sprintf("Score File: %s\nClinical File: %s\nOutput Prefix: %s\n", score.file, clinical.file, out.file), file = log)

cat(sprintf(" %% Reading FSM Score File\n"))
data <- as.matrix(read.table(score.file, header = TRUE, sep = "\t", check.names = FALSE))
if(all(is.na(data))){
  cat(sprintf(" !Empty FSM Score File\n"))
  quit(save = "no", status = 1)
}
data.random <- randomize(data)

cat(sprintf(" %% Running Column Factorizations\n\tSeed: Random\n\tMethod nsNMF\n\tnRun: 150\n"))
result <- nmf(data, 2:6, method = "nsNMF", seed = "random", nrun = 150, .options = "vp20")
save(result, file = sprintf("%s.nmf", out.file))

result.random <- nmf(data.random, 2:6, method = "nsNMF", seed = "random", nrun = 150, .options = "vp20")
save(result.random, file = sprintf("%s-random.nmf", out.file))

cat(sprintf(" %% Plotting Cluster Metrics\n"))
pdf(file = sprintf("%s-metrics.pdf", out.file), width = 16, height = 14)
plot(result, result.random)
invisible(dev.off())
cat(sprintf(" %% Results Saved on Data Objects: %s.nmf %s-random.nmf %s-metrics.pdf\n", out.file, out.file, out.file), file = log)

## Plotting Results ##
Clinical <- read.table(file = clinical.file, header = TRUE, sep = "\t", check.names = FALSE)
if(all(is.na(Clinical))){
  cat(sprintf(" !Empty Clinical File"))
  cat(sprintf(" !Empty Clinical File"), file = log)
  quit(save = "no", status = 1)
}

sortedIdx <- match(colnames(data), table = Clinical$`Patient ID`)
Clinical <- Clinical[sortedIdx, ]
rowIdx <- which(Clinical$`Followup Days` > thr)
if(length(rowIdx) > 0){
  Clinical$`Followup Days`[rowIdx] <- thr
  Clinical$`Event Status`[rowIdx] <- 0
}

if(nrow(Clinical) == 0){
  cat(sprintf(" !Error! Clinical File Sort Failed, Lines: 70 - 80\n"))
}

long.colors <- c("cadetblue", "cadetblue1", "chartreuse", 
                 "chartreuse4", "coral", "coral4", "darkblue",
                 "darkcyan","darkgreen", "darkred", "firebrick1", 
                 "khaki", "ivory4", "hotpink4", "seagreen4")

## Fixes Empty Lines ##
indices <- grep(colnames(Clinical), pattern = 'Patient|Barcode|Days|Event', perl = TRUE)
fixIdx <- which(!(1:ncol(Clinical) %in% indices) == TRUE)

for(c in fixIdx){
  for(r in 1:nrow(Clinical)){
    if(is.na(Clinical[r, c]) || Clinical[r, c] == ""){
      Clinical[r, c] <- NA
    }
  }
  Clinical[, c] <- factor(Clinical[, c])
}

colors.annotation <- list()
for(c in 1:ncol(Clinical[-indices])){
  num <- length(unique(as.matrix(Clinical[,-indices])[, c]))
  if(num > 8){
    colors.annotation[[c]] <- sample(long.colors, size = num, replace = FALSE)  
  }else{
    if(sample(0:1,1)){
      colors.annotation[[c]] <- brewer.pal(num, name = "Paired")    
    }else{
      colors.annotation[[c]] <- brewer.pal(num, name = "Set3")  
    }
  }
}

cat(sprintf(" %% Plotting Consensusmap\n"))
pdf(file = sprintf("NMF_Consensusmap_%s.pdf", out.file), width = 16, height = 12)
consensusmap(result, annCol = Clinical[, -indices], annColors = colors.annotation, 
             tracks = "consensus", fontsize = 6)
invisible(dev.off())

pdf(file = sprintf("NMF_Consensusmap_%s-R4.pdf", out.file), width = 16, height = 12)
consensusmap(result$fit$`4`, annCol = Clinical[, -indices], annColors = colors.annotation, 
             tracks = "consensus", fontsize = 6)
invisible(dev.off())

cat(sprintf(" %% Plotting Survival Curves\n"))
for(i in 1:length(result$fit)){
  res.predicted <- predict(result$fit[[i]], what = "consensus")
  localClinical <- cbind(Clinical, matrix(0, ncol = 1, nrow = nrow(Clinical)))
  colnames(localClinical)[ncol(localClinical)] <- "Group ID"
  localClinical$`Group ID` <- res.predicted
  
  ## Removes Patients With Survival Value < 30 (statistically irrelevant) ##
  rowIdx <- which(localClinical[,2] < 30)
  localClinical <- localClinical[-rowIdx, ]
  
  survival <- Surv(time = localClinical$`Followup Days`, event = localClinical$`Event Status`)
  Subtypes <- localClinical$`Group ID`
  data.surv <- survfit(survival ~ Subtypes)
  survival_result <- survdiff(survival ~ localClinical$`Group ID`, rho = 0)
  p.val <- 1 - pchisq(survival_result$chisq, length(survival_result$n) - 1)
  strata <- names(data.surv$strata)
  if(length(strata) > 8){
    surv.colors <- sample(long.colors, size = length(strata), replace = FALSE)
  }else if(length(strata) >= 3){
    surv.colors <- brewer.pal(length(strata), "Dark2")  
  }else{
    surv.colors <- c("dodgerblue", "firebrick")
  }
  plot <- ggsurv(data.surv, CI = FALSE, surv.col = surv.colors, cens.col = "black", back.white = FALSE) + 
    ggtitle("Kaplan-Meier Plot (GBM-MA)") + theme_bw() +
    geom_text(data = NULL, label = sprintf("p-value < %.2e", p.val), y = 0.9, x = (thr - 200), size = 8) +
    theme(plot.title = element_text(size = 24), axis.text = element_text(size = 18), axis.title = element_text(size = 18, face = "bold"),legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18))
  ggsave(file = sprintf("NMF_Survival_%s_Rank%d.pdf", out.file, length(strata)), plot = plot, width = 12, height = 10)
  
  pdf(file = sprintf("NMF_Loadings_%s_Rank%d.pdf", out.file, length(strata)), width = 12, height = 6)
  layout(cbind(1, 2))
  basismap(result$fit[[i]])
  coefmap(result$fit[[i]])
  invisible(dev.off())
}

cat(sprintf(" %%All Processes Completed.\n\n"))
cat(sprintf(" %%NMF Run Completed.\n\n"), file = log)
close(log)
quit(save = "no", status = 0)
