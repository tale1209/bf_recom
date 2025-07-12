library("Hapi")

rm(list = ls()) 
setwd("path/to/wd")

missing <- 5
impute <- 5     

gmtDa <- read.table("bubble_matrix_filter", header = TRUE)
gmtFrame <- hapiFrameSelection(gmt = gmtDa, n = missing)     
imputedFrame <- hapiImupte(gmt = gmtFrame, nSPT = impute, allowNA = 0)     
draftHap <- hapiPhase(gmt = imputedFrame)    

cvCluster <- hapiCVCluster(draftHap = draftHap, minDistance = 1000000, cvlink = 2)
filter <- c()
for (i in 1:nrow(cvCluster)) {filter <- c(filter, which (rownames(draftHap) >= cvCluster$left[i] & rownames(draftHap) <= cvCluster$right[i]))}
if (length(filter) > 0) {imputedFrame <- imputedFrame[-filter, ];draftHap <- hapiPhase(imputedFrame)}

finalDraft <- hapiBlockMPR(draftHap = draftHap, gmtFrame = gmtFrame, cvlink = 2, smallBlock = 0)
consensusHap <- hapiAssemble(draftHap = finalDraft, gmt = gmtDa)
hap1 <- sum(consensusHap$hap1==0)
hap2 <- sum(consensusHap$hap1==1)
hap7 <- sum(consensusHap$hap1==7)
max(hap1, hap2)/sum(hap1, hap2)

hap <- consensusHap[,1:2]
sample <- gmtDa[,1:length(gmtDa)]         
cvOutput <- hapiIdentifyCV(hap = hap, gmt = sample)
cvOutput_result <- as.data.frame(cvOutput)
hap_result <- as.data.frame(hap)
write.table(cvOutput_result,file="crossover")
write.table(hap_result,file="haplotype")
