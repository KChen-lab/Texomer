### R code from vignette source 'TitanCNA.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: TitanCNA.Rnw:75-80
###################################################
library(TitanCNA)
infile <- system.file("extdata", "test_alleleCounts_chr2.txt", 
                      package = "TitanCNA")
data <- loadAlleleCounts(infile)
names(data)


###################################################
### code chunk number 2: TitanCNA.Rnw:88-93
###################################################
numClusters <- 2
params <- loadDefaultParameters(copyNumber = 5, 
                                numberClonalClusters = numClusters,
                                symmetric = TRUE, data = data)
params


###################################################
### code chunk number 3: TitanCNA.Rnw:107-109 (eval = FALSE)
###################################################
## params$ploidyParams$phi_0 <- 2 # for diploid or
## params$ploidyParams$phi_0 <- 4 # for tetraploid/ployploid


###################################################
### code chunk number 4: TitanCNA.Rnw:114-117
###################################################
K <- length(params$genotypeParams$alphaKHyper)
params$genotypeParams$alphaKHyper <- rep(500, K)
params$ploidyParams$phi_0 <- 1.5


###################################################
### code chunk number 5: TitanCNA.Rnw:125-131
###################################################
tumWig <- system.file("extdata", "test_tum_chr2.wig", package = "TitanCNA")
normWig <- system.file("extdata", "test_norm_chr2.wig", package = "TitanCNA")
gc <- system.file("extdata", "gc_chr2.wig", package = "TitanCNA")
map <- system.file("extdata", "map_chr2.wig", package = "TitanCNA")
cnData <- correctReadDepth(tumWig, normWig, gc, map)
head(cnData)


###################################################
### code chunk number 6: TitanCNA.Rnw:138-140
###################################################
logR <- getPositionOverlap(data$chr, data$posn, cnData)
data$logR <- log(2^logR)  #transform the log ratio to natural logs


###################################################
### code chunk number 7: TitanCNA.Rnw:150-152
###################################################
data <- filterData(data, c(1:22, "X", "Y"), minDepth = 10, maxDepth = 200, 
                   positionList = NULL, centromere = NULL, centromere.flankLength = 10000)


###################################################
### code chunk number 8: TitanCNA.Rnw:159-161 (eval = FALSE)
###################################################
## library(doMC)
## registerDoMC(cores = 4) #use 4 cores on a single machine


###################################################
### code chunk number 9: TitanCNA.Rnw:176-186
###################################################
convergeParams <- runEMclonalCN(data, gParams = params$genotypeParams, 
                                nParams = params$normalParams, 
                                pParams = params$ploidyParams, 
                                sParams = params$cellPrevParams, 
                                maxiter = 3, maxiterUpdate = 50, 
                                useOutlierState = FALSE, txnExpLen = 1e15, 
                                txnZstrength = 5e5, 
                                normalEstimateMethod = "map", 
                                estimateS = TRUE, estimatePloidy = TRUE)
names(convergeParams)


###################################################
### code chunk number 10: TitanCNA.Rnw:206-208
###################################################
optimalPath <- viterbiClonalCN(data, convergeParams)
head(optimalPath)


###################################################
### code chunk number 11: TitanCNA.Rnw:219-223
###################################################
results <- outputTitanResults(data, convergeParams, optimalPath,
                              filename = NULL, posteriorProbs = FALSE,
                              subcloneProfiles = TRUE)
head(results)


###################################################
### code chunk number 12: TitanCNA.Rnw:250-254 (eval = FALSE)
###################################################
## corrResults <- removeEmptyClusters(convergeParams, results, proportionThreshold = 0.001, 
## 																	 proportionThresholdClonal = 0.3)
## convergeParams <- corrResults$convergeParams
## results <- corrResults$results


###################################################
### code chunk number 13: TitanCNA.Rnw:259-261 (eval = FALSE)
###################################################
## outparam <- paste("test_cluster02_params.txt", sep = "")
## outputModelParameters(convergeParams, results, outparam, S_Dbw.scale = 1)


###################################################
### code chunk number 14: TitanCNA.Rnw:281-282 (eval = FALSE)
###################################################
## outputModelParameters(convergeParams, results, outparam, S_Dbw.scale = 10)


###################################################
### code chunk number 15: TitanCNA.Rnw:292-298
###################################################
ploidy <- tail(convergeParams$phi, 1)
ploidy
normal <- tail(convergeParams$n, 1)
normal
plotCNlogRByChr(results, chr = 2, ploidy = ploidy, normal = normal, 
		ylim = c(-2, 2), cex = 0.25, xlab = "", main = "Chr 2")


###################################################
### code chunk number 16: TitanCNA.Rnw:311-313
###################################################
plotAllelicRatio(results, chr = 2, ylim = c(0, 1), cex = 0.25, 
                 xlab = "", main = "Chr 2")


###################################################
### code chunk number 17: TitanCNA.Rnw:327-332
###################################################
norm <- tail(convergeParams$n, 1) 
norm # estimated normal contamination
1 - convergeParams$s[, ncol(convergeParams$s)] # estimated cellular prevalence 
plotClonalFrequency(results, chr = 2, normal = norm, ylim = c(0, 1), 
                    cex = 0.25, xlab = "", main = "Chr 2")


###################################################
### code chunk number 18: TitanCNA.Rnw:341-342
###################################################
plotSubcloneProfiles(results, chr = 2, cex = 1, spacing = 2, main = "Chr 2")


###################################################
### code chunk number 19: TitanCNA.Rnw:351-352
###################################################
toLatex(sessionInfo())


