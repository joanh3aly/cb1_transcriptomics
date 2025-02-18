setwd('C:/Users/joanhealy24/OneDrive - Royal College of Surgeons in Ireland/Documents/CB1/transcriptomics_assignment')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("edgeR")
BiocManager::install("Rsubread")


library(Rsubread)
library(limma)
library(edgeR)
library(clusterProfiler)

load("C:/Users/joanhealy24/OneDrive - Royal College of Surgeons in Ireland/Documents/CB1/transcriptomics_assignment/gene.Rdata")

barplot(as.matrix(t(cbind(propmap[,"NumMapped"]*1e-6, (propmap[,"NumTotal"]-propmap[,"NumMapped"])*1e-6))),ylab="Number of Reads (millions)", 
          names = Targets[1:6,1], las=3, cex.names=0.3, cex.axis=0.7, density = c(80,0), legend = c("Number of Mapped Reads", "Number of Unmapped Reads"), 
          args.legend = list(x = "topright",cex = 0.5, bty = "n"))

barplot(colSums(Counts)*1e-6, names = Targets[1:6,1], ylab="Library size (millions)", las = 2, cex.names = 0.3)
legend("topleft", legend = c("1.5","2.5","3.5"), pch = 15, cex = 0.7)


load("C:/Users/joanhealy24/OneDrive - Royal College of Surgeons in Ireland/Documents/CB1/transcriptomics_assignment/norm.Rdata")
head(y)
plotMDS(y, labels=Targets[1:6,1], gene.selection="common", prior.count=5)

library(limma)
library(edgeR)

design <- model.matrix(~ 0+factor(c(rep(1,3),rep(2,3))))
colnames(design) <- c("DMSO", "APR246")
contrast.matrix <- makeContrasts(APR246-DMSO, levels=design)

# v <- voom(y, design)
# fit <- lmFit(v, design)
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit3 <- treat(fit2, lfc=log2(1.2))
# topTreat(fit3)
# dim(topTreat(fit3, number=1000, lfc=1.2, p.value=0.05))

v <- voom(y, design)
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)
topTable(fit3)

dim(topTable(fit3, number=1000, lfc=1.2, p.value=0.05))


library(clusterProfiler)
fit3
de <- topTreat(fit3, number=1000, p.value=0.05)[,1]
yy <- enrichKEGG(de, organism = "hsa", pvalueCutoff=0.5)
head(yy)

