### R code from vignette source 'copynumber.rnw'

###################################################
### code chunk number 1: copynumber.rnw:77-78
###################################################
library(copynumber)


###################################################
### code chunk number 2: copynumber.rnw:83-84
###################################################
data(lymphoma)


###################################################
### code chunk number 3: copynumber.rnw:87-89
###################################################
sub.lymphoma <- subsetData(data=lymphoma,sample=1:3) 
sub.lymphoma[1:10,]


###################################################
### code chunk number 4: copynumber.rnw:93-94
###################################################
lymph.wins <- winsorize(data=sub.lymphoma,verbose=FALSE)


###################################################
### code chunk number 5: copynumber.rnw:97-98
###################################################
lymph.wins[1:10,]


###################################################
### code chunk number 6: copynumber.rnw:101-103
###################################################
wins.res <- winsorize(data=sub.lymphoma,return.outliers=TRUE,verbose=FALSE)
wins.res$wins.outliers[1:10,]


###################################################
### code chunk number 7: copynumber.rnw:108-109
###################################################
single.seg <- pcf(data=lymph.wins,gamma=12,verbose=FALSE)


###################################################
### code chunk number 8: copynumber.rnw:113-114
###################################################
head(single.seg)


###################################################
### code chunk number 9: genomeFile
###################################################
png("plotGenome.png", width = 900, height = 500)


###################################################
### code chunk number 10: plotGenome
###################################################
plotGenome(data=sub.lymphoma,segments=single.seg,sample=1,cex=3) 


###################################################
### code chunk number 11: genomeClose
###################################################
null <- dev.off()


###################################################
### code chunk number 12: sampleFile
###################################################
png("plotSample.png", width = 900, height = 700)


###################################################
### code chunk number 13: plotSample
###################################################
plotSample(data=sub.lymphoma,segments=single.seg,layout=c(5,5),sample=1,cex=3)


###################################################
### code chunk number 14: sampleClose
###################################################
null <- dev.off()


###################################################
### code chunk number 15: copynumber.rnw:165-167
###################################################
multi.seg <- multipcf(data=lymph.wins,verbose=FALSE)
head(multi.seg)


###################################################
### code chunk number 16: chromFile
###################################################
png("chromPlot.png", width = 480, height=480)


###################################################
### code chunk number 17: chromPlot
###################################################
plotChrom(data=lymph.wins,segments=multi.seg,layout=c(3,1),chrom=1)


###################################################
### code chunk number 18: chromClose
###################################################
null <- dev.off()


###################################################
### code chunk number 19: copynumber.rnw:196-198
###################################################
data(logR)
data(BAF)


###################################################
### code chunk number 20: copynumber.rnw:201-202
###################################################
logR.wins <- winsorize(logR,verbose=FALSE)


###################################################
### code chunk number 21: copynumber.rnw:205-207
###################################################
allele.seg <- aspcf(logR.wins,BAF,verbose=FALSE)
head(allele.seg) 


###################################################
### code chunk number 22: aspcfFile
###################################################
png("plotAllele.png", width = 800, height = 500)


###################################################
### code chunk number 23: chromPlot
###################################################
plotAllele(logR,BAF,allele.seg,sample=1,chrom=c(1:4),layout=c(2,2))


###################################################
### code chunk number 24: aspcfClose
###################################################
null <- dev.off()


###################################################
### code chunk number 25: copynumber.rnw:237-238
###################################################
lymphoma.res <- pcf(data=lymphoma,gamma=12,verbose=FALSE)


###################################################
### code chunk number 26: freqFile
###################################################
png("freqPlot.png", width = 800, height = 400)


###################################################
### code chunk number 27: freqPlot
###################################################
plotFreq(segments=lymphoma.res,thres.gain=0.2,thres.loss=-0.1)


###################################################
### code chunk number 28: freqClose
###################################################
null <- dev.off()


###################################################
### code chunk number 29: copynumber.rnw:269-275
###################################################
chr.from <- c(2,12,4)
pos.from <- c(168754669,847879349,121809306)
chr.to <- c(14,21,17)
pos.to <- c(6147539,301955563,12364465)
cl <- c(1,1,2)
arcs <- cbind(chr.from,pos.from,chr.to,pos.to,cl)  


###################################################
### code chunk number 30: circleFile
###################################################
png("circlePlot.png", width = 600, height = 600)


###################################################
### code chunk number 31: circlePlot
###################################################
plotCircle(segments=lymphoma.res,thres.gain=0.15,arcs=arcs) 


###################################################
### code chunk number 32: circleClose
###################################################
null <- dev.off()


###################################################
### code chunk number 33: heatFile
###################################################
png("plotHeatmap.png", width = 800, height = 400)


###################################################
### code chunk number 34: plotHeatmap
###################################################
plotHeatmap(segments=lymphoma.res,upper.lim=0.3)


###################################################
### code chunk number 35: heatClose
###################################################
null <- dev.off()


###################################################
### code chunk number 36: aberrationFile
###################################################
png("plotAberration.png", width = 800, height = 400)


###################################################
### code chunk number 37: plotAberration
###################################################
plotAberration(segments=lymphoma.res,thres.gain=0.2)


###################################################
### code chunk number 38: aberrationClose
###################################################
null <- dev.off()


###################################################
### code chunk number 39: gammaFile
###################################################
png("plotGamma.png", width = 800, height = 600)


###################################################
### code chunk number 40: plotGamma
###################################################
data(micma)
plotGamma(micma,chrom=17,cex=3)


###################################################
### code chunk number 41: gammaClose
###################################################
null <- dev.off()


