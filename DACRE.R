# TODO: Add comment
# 
# Author: FWang9
###############################################################################
#R version needs >= 3.2.0
#Rscript segCN.R parameters
#the order of parameters is germline input,somaticinput,RNAinput,script path,outpath
#the outpath is optional and others are neceaasry
args<-commandArgs(T)
SNPinput=args[1]
somaticinput=args[2]
RNAinput=args[3]
path=args[4]
outpath=args[5]
optindex=args[6]
library(doMC)
library(bbmle,lib.loc=paste(path,"/library",sep=""))
library(emdbook,lib.loc=paste(path,"/library",sep=""))
library(copynumber,lib.loc=paste(path,"/library",sep=""))
library(TitanCNA,lib.loc=paste(path,"/library",sep=""))
library(facets,lib.loc=paste(path,"/library",sep=""))
library(mixtools,lib.loc=paste(path,"/library",sep=""))
library(sfsmisc)
source(paste(path,"/ASCAT.R",sep=""))
source(paste(path,"/DNAfunction.R",sep=""))
source(paste(path,"/RNAfunction.R",sep=""))
source(paste(path,"/sequenza.R",sep=""))
#setwd(outpath)
#dir.create("temp")
temppath=paste(outpath,"/temp",sep="")
caseid=strsplit(SNPinput,split="/")[[1]][length(strsplit(SNPinput,split="/")[[1]])]
if (length(args)==6){
	if (optindex==1){
		DNAout=try(DNArun(SNPinput=SNPinput,somaticinput=somaticinput,sample=caseid,temppath=temppath),silent=TRUE)
	}else{
		DNAout=try(DNArun1(SNPinput=SNPinput,somaticinput=somaticinput,sample=caseid,temppath=temppath),silent=TRUE)
	}
	if (!is.null(names(DNAout))){
		DNAout=Heterogeneity(DNAout)
		heterogeneity=DNAout$heterogeneity
		DNApurity=DNAout$alpha
		ploidy=DNAout$ploidy
		summaryres=rbind(c("DNApurity",DNApurity),c("Heterogeneity",heterogeneity),c("Ploidy",ploidy))
		if (RNAinput=="NA"){
			segment=DNAout$segment
			somatic=DNAout$somatic
			somatic=data.frame(chr=somatic$chr,pos=somatic$pos,ref=somatic$ref,alt=somatic$alt,refNumN=somatic$refNumN,altNumN=somatic$altNumN,refNumT=somatic$refNumT,altNumT=somatic$altNumT,altD=somatic$SMCN)
			write.table(summaryres,paste(outpath,"/output.summaryres.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
			write.table(segment,paste(outpath,"/output.segment.txt",sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
			write.table(somatic,paste(outpath,"/output.mutation.txt",sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
		}else{
			print.noquote("Inferring at RNA level!")
			RNAout=RNArun(SNPinput=SNPinput,RNAinput=RNAinput,DNAout=DNAout)
			segment=RNAout$segment
			RNApurity=RNAout$alpha
			germline=RNAout$germline
			summaryres=rbind(summaryres,c("RNApurity",RNApurity))
			resout=c()
    		resout=wildtype(data=germline,segment=segment,type="germline",resout=resout)
    		if ("somatic" %in% names(RNAout)){
    			somatic=RNAout$somatic
    			resout=wildtype(data=somatic,segment=segment,type="somatic",resout=resout)
    		}
			res=DACRE(resout=resout)
			segment=data.frame(chr=segment$chr,start=segment$startpos,end=segment$endpos,Dmajor=segment$nMajor,Dminor=segment$nMinor,Rmajor=segment$Rmajor,Rminor=segment$Rminor,RTEL=segment$RTCN,BayesP=segment$BayesP)
			write.table(summaryres,paste(outpath,"/output.summaryres.txt",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
			write.table(segment,paste(outpath,"/output.segment.txt",sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
			write.table(res,paste(outpath,"/output.mutation.txt",sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
		}
		
	}else{
		print.noquote(paste(sample,": No output!",sep=""))
	}
}
