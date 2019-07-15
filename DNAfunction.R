# TODO: Add comment
#
# Author: FWang9
###############################################################################


############Using ASCAT
############segmentation and estimate DNA copy number based on germline mutation
############the formate of inputdata were chr, position, refallele counts in tumor, altAllele counts in tumor, refallele counts in normal, altAllele counts in normal
############tab seprate and without headline
############load ASCAT function
############install bbmle and emdbook package
##############################################################################################################
####part1.run ASCAT to get segmentation and allelic copy number in DNA level
###ASCAT in: get the inputdata of ASCAT
ASCATin<-function(data,sample){
	Me=median(data$TlogR)
	M=max(data$TlogR)
	m=min(data$TlogR)
	data$Tlog[data$TlogR>=Me]=(data$TlogR[data$TlogR>=Me]-Me)/(M-Me)
	data$Tlog[data$TlogR<Me]=(data$TlogR[data$TlogR<Me]-Me)/(Me-m)
	data$Nlog=1
	tumor_logR=data.frame(SNP=paste("SNP",c(1:dim(data)[1]),sep=""),chr=as.character(data$chr),pos=data$pos,sample=data$Tlog)
	tumor_BAF=data.frame(SNP=paste("SNP",c(1:dim(data)[1]),sep=""),chr=as.character(data$chr),pos=data$pos,sample=data$tfrac)
	normal_logR=data.frame(SNP=paste("SNP",c(1:dim(data)[1]),sep=""),chr=as.character(data$chr),pos=data$pos,sample=data$Nlog)
	normal_BAF=data.frame(SNP=paste("SNP",c(1:dim(data)[1]),sep=""),chr=as.character(data$chr),pos=data$pos,sample=data$nfrac)
	colnames(tumor_logR)=c("SNP","chr","pos",sample)
	colnames(tumor_BAF)=c("SNP","chr","pos",sample)
	colnames(normal_logR)=c("SNP","chr","pos",sample)
	colnames(normal_BAF)=c("SNP","chr","pos",sample)
	tumor_logRname=paste(sample,"_tumorlogR.txt",sep="")
	tumor_BAFname=paste(sample,"_tumorBAF.txt",sep="")
	normal_logRname=paste(sample,"_normallogR.txt",sep="")
	normal_BAFname=paste(sample,"_normalBAF.txt",sep="")
	write.table(tumor_logR,tumor_logRname,sep="\t",col.names=TRUE,row.names = FALSE,quote = FALSE)
	write.table(tumor_BAF,tumor_BAFname,sep="\t",col.names=TRUE,row.names = FALSE,quote = FALSE)
	write.table(normal_logR,normal_logRname,sep="\t",col.names=TRUE,row.names = FALSE,quote = FALSE)
	write.table(normal_BAF,normal_BAFname,sep="\t",col.names=TRUE,row.names = FALSE,quote = FALSE)
	ASCATdata=c(tumor_logRname,tumor_BAFname,normal_logRname,normal_BAFname)
	return (ASCATdata)
}

###ASCATout: get the output of ASCAT
ASCATout<-function(ASCATdata,sample,chromosome){
	tumor_logRname=ASCATdata[1]
	tumor_BAFname=ASCATdata[2]
	normal_logRname=ASCATdata[3]
	normal_BAFname=ASCATdata[4]
	ascat.bc = ascat.loadData(tumor_logRname,tumor_BAFname,normal_logRname,normal_BAFname,chrs=chromosome)
	ascat.bc = ascat.aspcf(ascat.bc)
	ascat.output = ascat.runAscat(ascat.bc)
	alpha=ascat.output$aberrantcellfraction
	segment=ascat.output$segments
	ACN=data.frame(nMajor=segment$nMajor,nMinor=segment$nMinor)
	segment$nMajor=apply(ACN,1,max)
	segment$nMinor=apply(ACN,1,min)
	segment=data.frame(chr=segment$chr,startpos=segment$startpos,endpos=segment$endpos,nMajor=segment$nMajor,nMinor=segment$nMinor)
	alpha=ascat.output$aberrantcellfraction
	ploidy=ascat.output$ploidy
	DNAout=list(alpha=alpha,segment=segment,ploidy=ploidy)
	#save(DNAout,file=paste(sample,".somaticCN.Rdata",sep=""))
	return(DNAout)
}
####part 1 end
#############################################################################################################
#############################################################################################################
####part 2: estimate the copy number of somatic mutation and iterative optimization
#MCN: get the integer copy number of somatic mutation
MCN<-function(alpha,segment,somatic){
	for (j in 1:dim(segment)[1]){
		tempsomatic=somatic[somatic$chr==as.character(segment$chr[j])&somatic$pos>=segment$startpos[j]&somatic$pos<=segment$endpos[j],]
		tempsomatic$SMCN=tempsomatic$tfrac*(2*(1-alpha)+alpha*(segment$nMajor[j]+segment$nMinor[j]))/alpha
		if (j==1){
			somaticnew=tempsomatic
		}else{
			somaticnew=rbind(somaticnew,tempsomatic)
		}
	}
	somaticnew$SACN=round(somaticnew$SMCN)
	return(somaticnew)
}
#MutPlot: plot the copy number of somatic mutation
MutPlot<-function(segment,summaryres,mutation,plotname="SNV.CNV.pdf"){
	segment=read.csv(segment,sep="\t")
	summaryres=read.csv(summaryres,sep="\t",header=F)
	data=read.csv(mutation,sep="\t")
	purity=round(summaryres[1,2],2)
	Heterogeneity=round(summaryres[2,2],2)
	ploidy=round(summaryres[3,2],2)
	pdf(plotname,width=20,height=7)
	d=sum(as.numeric(segment$end)-as.numeric(segment$start))
	nMajor=max(segment$Dmajor,segment$Dminor,data$altD)
	plot(0,0,col="white",xlim=c(0,d),ylim=c(0,nMajor),xlab="",ylab="Absolute copy number",main=paste("Purity = ",purity,", Heterogeneity = ",Heterogeneity,", Ploidy = ",ploidy,sep=""),axes=FALSE)
	chro=as.character(unique(segment$chr))
	axisindex=0
	for (j in 1:length(chro)){
		subseg=segment[segment$chr==chro[j],]
		subd=sum(subseg$end-subseg$start)
		subsomatic=data[data$chr==chro[j]&data$type=="Somatic",]
		for (k in 1:dim(subseg)[1]){
			subsom=subsomatic[subsomatic$pos>=subseg$start[k]&subsomatic$pos<=subseg$end[k],]
			subsom$index=subsom$pos-subseg$start[k]
			#points(subsom$index+axisindex,subsom$SACN,pch=20,col=rgb(0,0,1,alpha=0.7),cex=1.2)
			points(subsom$index+axisindex,subsom$altD,pch=20,col=rgb(0,0,1,alpha=0.7),cex=1.2)
			if (subseg$Dmajor[k]==subseg$Dminor[k]){
				x=c(axisindex,axisindex+subseg$end[k]-subseg$start[k])
				y1=rep(subseg$Dmajor[k]-0.1,length(x))
				lines(x,y1,col="purple",lwd=2)
				y2=rep(subseg$Dmajor[k]+0.1,length(x))
				lines(x,y2,col="purple",lwd=2)
				z=c(subseg$Dmajor[k]-0.1,subseg$Dmajor[k]+0.1)
				x1=c(axisindex,axisindex)
				lines(x1,z,col="purple",lwd=2)
				x2=c(axisindex+subseg$end[k]-subseg$start[k],axisindex+subseg$end[k]-subseg$start[k])
				lines(x2,z,col="purple",lwd=2)
			}else{
				x=c(axisindex,axisindex+subseg$end[k]-subseg$start[k])
				y1=rep(subseg$Dmajor[k]-0.1,length(x))
				lines(x,y1,col="red",lwd=2)
				y2=rep(subseg$Dmajor[k]+0.1,length(x))
				lines(x,y2,col="red",lwd=2)
				z=c(subseg$Dmajor[k]-0.1,subseg$Dmajor[k]+0.1)
				x1=c(axisindex,axisindex)
				lines(x1,z,col="red",lwd=2)
				x2=c(axisindex+subseg$end[k]-subseg$start[k],axisindex+subseg$end[k]-subseg$start[k])
				lines(x2,z,col="red",lwd=2)
				y11=rep(subseg$Dminor[k]-0.1,length(x))
				lines(x,y11,col="green",lwd=2)
				y21=rep(subseg$Dminor[k]+0.1,length(x))
				lines(x,y21,col="green",lwd=2)
				z=c(subseg$Dminor[k]-0.1,subseg$Dminor[k]+0.1)
				lines(x1,z,col="green",lwd=2)
				lines(x2,z,col="green",lwd=2)
			}
			axisindex=axisindex+subseg$end[k]-subseg$start[k]+1
		}
		abline(v=axisindex-1,col="gray")
		text((axisindex-1-sum(subseg$end-subseg$start)+axisindex-1)/2,nMajor,substr(chro[j],4,nchar(chro[j])))
	}
	axis(side=2)
	dev.off()
}

#SNPplot<-function(segment,data,plotname,index,alpha){
#	pdf(plotname,width=20,height=7)
#	d=sum(segment$endpos-segment$startpos)
#	nMajor=max(segment$nMajor,segment$nMinor,data$SMCN)
#	plot(0,0,col="white",xlim=c(0,d),ylim=c(0,nMajor),xlab="",ylab="Absolute copy number",main="",axes=FALSE)
#	chro=unique(segment$chr)
#	axisindex=0
#	for (j in 1:length(chro)){
#		subseg=segment[segment$chr==chro[j],]
#		subd=sum(subseg$endpos-subseg$startpos)
#		subsomatic=data[data$chr==chro[j],]
#		for (k in 1:dim(subseg)[1]){
#			subsom=subsomatic[subsomatic$pos>=subseg$startpos[k]&subsomatic$pos<=subseg$endpos[k],]
#			TCN=subseg$nMajor[k]+subseg$nMinor[k]
#			subsom$SMCN=(subsom$tfrac*(2*(1-alpha)+alpha*TCN)-(1-alpha))/alpha
#			subsom$SACN=round(subsom$SMCN)
#			subsom$index=subsom$pos-subseg$startpos[k]
#			if (index==1){
#				points(subsom$index+axisindex,subsom$SACN,pch=20,col=rgb(0,0,1,alpha=0.7),cex=1.2)
#			}
#			else{
#				points(subsom$index+axisindex,subsom$SMCN,pch=20,col=rgb(0,0,1,alpha=0.7),cex=1.2)
#			}
#			if (subseg$nMajor[k]==subseg$nMinor[k]){
#				x=c(axisindex,axisindex+subseg$endpos[k]-subseg$startpos[k])
#				y1=rep(subseg$nMajor[k]-0.1,length(x))
#				lines(x,y1,col="purple",lwd=2)
#				y2=rep(subseg$nMajor[k]+0.1,length(x))
#				lines(x,y2,col="purple",lwd=2)
#				z=c(subseg$nMajor[k]-0.1,subseg$nMajor[k]+0.1)
#				x1=c(axisindex,axisindex)
#				lines(x1,z,col="purple",lwd=2)
#				x2=c(axisindex+subseg$endpos[k]-subseg$startpos[k],axisindex+subseg$endpos[k]-subseg$startpos[k])
#				lines(x2,z,col="purple",lwd=2)
#			}else{
#				x=c(axisindex,axisindex+subseg$endpos[k]-subseg$startpos[k])
#				y1=rep(subseg$nMajor[k]-0.1,length(x))
#				lines(x,y1,col="red",lwd=2)
#				y2=rep(subseg$nMajor[k]+0.1,length(x))
#				lines(x,y2,col="red",lwd=2)
#				z=c(subseg$nMajor[k]-0.1,subseg$nMajor[k]+0.1)
#				x1=c(axisindex,axisindex)
#				lines(x1,z,col="red",lwd=2)
#				x2=c(axisindex+subseg$endpos[k]-subseg$startpos[k],axisindex+subseg$endpos[k]-subseg$startpos[k])
#				lines(x2,z,col="red",lwd=2)
#				y11=rep(subseg$nMinor[k]-0.1,length(x))
#				lines(x,y11,col="green",lwd=2)
#				y21=rep(subseg$nMinor[k]+0.1,length(x))
#				lines(x,y21,col="green",lwd=2)
#				z=c(subseg$nMinor[k]-0.1,subseg$nMinor[k]+0.1)
#				lines(x1,z,col="green",lwd=2)
#				lines(x2,z,col="green",lwd=2)
#			}
#			axisindex=axisindex+subseg$endpos[k]-subseg$startpos[k]+1
#		}
#		abline(v=axisindex-1,col="gray")
#		text((axisindex-1-sum(subseg$endpos-subseg$startpos)+axisindex-1)/2,nMajor,substr(chro[j],4,nchar(chro[j])))
#	}
#	axis(side=2)
#	dev.off()
#}
#randomBinom: random sampling based on binomial distribution
randomBinom<-function(y,times){
	return(rbinom(times,y[1],y[2]))
}
#ransomSample: random sample based on times
randomSample<-function(data,times){
	if (dim(data)[1]>0){
		random=apply(data,1,randomBinom,times=times)
		return(t(random))
	}
}
#RMCN: get the allele copy number of random sampling
RMCN<-function(R,somatic,alpha,segment){
	randomsomatic=data.frame(chr=somatic$chr,pos=somatic$pos,refNumT=somatic$refNumT,altNumT=R)
	randomsomatic$tfrac=randomsomatic$altNumT/(randomsomatic$altNumT+randomsomatic$refNumT)
	randomres=MCN(alpha,segment,somatic=randomsomatic)
	return(randomres$SMCN)
}
#difCN: calculated the difference between allelic copy number and integer copy number
difCN<-function(x){
	dif=abs(x-round(x))
	return(dif)
}
####calculated the fraction of integer copy number = 0
#subclone: calculated the number of somatic mutation with copy number = 0
subclone<-function(x){
	return(sum(round(x)==0))
}
###somaticPurity: calculated the tumor purity based on the copy number of somatic mutation
somaticPurity<-function(somatic,segment,DNAalpha){
	alpha_up=c()
	for (j in 1:dim(segment)[1]){
		subdata=somatic[somatic$chr==as.character(segment$chr[j])&somatic$pos>=segment$startpos[j]&somatic$pos<=segment$endpos[j],]
		if (dim(subdata)[1]>0){
			TCN=segment$nMajor[j]+segment$nMinor[j]
			alpha1=2*subdata$tfrac/(2*subdata$tfrac+subdata$SACN-subdata$tfrac*TCN)
			alpha_up=c(alpha_up,alpha1)
		}

	}
	alpha_up=alpha_up[!is.na(alpha_up)&alpha_up!=Inf]
	if (length(alpha_up)<=1){
		if (length(alpha_up)==0){
			alpha1=DNAalpha
			}else{
				if (alpha_up>=0 & alpha_up <=1){
					alpha1=alpha_up
				}else{
					alpha1=DNAalpha
				}
			}

	}else{
		alpha1=density(alpha_up)$x[which.max(density(alpha_up)$y)]
		if (alpha1 > 1 | alpha1 < 0){
			alpha_up=alpha_up[alpha_up >=0 & alpha_up <= 1]
			if (length(alpha_up)>0){
				if (length(alpha_up)==1){
					alpha1=alpha_up
				}else{
					alpha1=density(alpha_up)$x[which.max(density(alpha_up)$y)]
				}
			}else{
				alpha1=DNAalpha
			}

		}
	}

	return(alpha1)
}
##CNopt: determine the optimal copy number for each segment among mang candiate CNs
CNpeak<-function(CN){
	peak=c()
	peakvalue=c()
	if (max(CN)<=0){
		CN=CN+abs(min(CN))
	}
	CN=CN[CN>0]
	if (length(CN)>3){
		k=1
		x=density(CN)$x
		y=density(CN)$y
		for (i in 2:(length(y)-1)){
			if (y[i]>=y[i-1]&y[i]>=y[i+1]){
				peak[k]=x[i]
				peakvalue[k]=y[i]
				k=k+1
			}
		}
		peak=peak[peakvalue> (max(peakvalue)/10)]
	}else{
		peak=unique(CN)
	}
	return(peak)
}
CNopt<-function(data,expBAF){
	dbetabinom(data$y,f1,data$N,theta,log=FALSE)
}
##upACN: update the copy number of segment and ploidy based on germline mutation
upACN<-function(data,segment,alpha){
	ploidy=c()
	for (j in 1:dim(segment)[1]){
		subseg=segment[j,]
		subdata=data[data$chr==as.character(segment$chr[j])&data$pos>=segment$startpos[j]&data$pos<=segment$endpo[j],]
		if (dim(subdata)[1]>=10){
			TCN=round((2*subdata$Ratio-2*(1-alpha))/alpha)
			canTCN=round(CNpeak(TCN))
			lCN=Inf
			sCN=c()
			if (length(canTCN)!=0){
				for (i in 1:length(canTCN)){
					nMajor=round((2*subdata$BAF*(1-alpha)+alpha*subdata$BAF*canTCN[i]-(1-alpha))/alpha)
					canMajor=round(CNpeak(nMajor))
					if (length(canMajor)>0){
						for (k in 1:length(canMajor)){
							expBAF=((1-alpha)+alpha*canMajor[k])/(2*(1-alpha)+alpha*canTCN[i])
							if (expBAF <= 1 & expBAF >=0){
								lCN1=-sum(dbinom(round(subdata$BAF*(subdata$refNumT+subdata$altNumT)),(subdata$refNumT+subdata$altNumT),expBAF,log=TRUE))
								if (lCN1 < lCN){
									sCN=c(canMajor[k],canTCN[i])
									lCN=lCN1
								}
							}

						}
					}
				}
				if (length(sCN)!=0){
					subseg$nMajor=sCN[1]
					if (sCN[2]-sCN[1]<0){
						subseg$nMinor=0
					}else{
						subseg$nMinor=sCN[2]-sCN[1]
					}
				}
			}
		}
		if (j==1){
			upsegment=subseg
		}else{
			upsegment=rbind(upsegment,subseg)
		}
		ploidy=c(ploidy,rep(subseg$nMajor+subseg$nMinor,length=dim(subdata)[1]))
	}
	upres=list(upsegment=upsegment,ploidy=mean(ploidy))
	return(upres)
}

###iterOPT: output optimal result based on Iterative optimization
iterOPT<-function(SNP,somatic,segment,realDIF,realsub,randomdata,times,PD,PS,DNAalpha,cutoff){
	out=list()
	out1=list()
	out2=list()
	run=1
	PD2=PD
	PS2=PS
	while(run < cutoff){
		alpha1=somaticPurity(somatic=somatic,segment=segment,DNAalpha=DNAalpha)
		res1=upACN(data=SNP,segment=segment,alpha=alpha1)
		ploidy1=res1$ploidy
		segment1=res1$upsegment
		somaticres1=MCN(alpha=alpha1,segment=segment1,somatic=somatic)
		somaticres1=somaticres1[match(paste(somatic$chr,somatic$pos,sep=":"),paste(somaticres1$chr,somaticres1$pos,sep=":")),]
		realDIF1=sum(abs(somaticres1$SMCN-somaticres1$SACN))/dim(somaticres1)[1]
		registerDoMC(cores = 4)
		randomout1=foreach (j = 1:dim(randomdata)[2], .combine=cbind) %dopar% RMCN(R=randomdata[,j],somatic=somaticres1,alpha=alpha1,segment=segment1)
		randomDIF1=apply(randomout1,2,difCN)
		aveDIF1=apply(randomDIF1,2,sum)/dim(somaticres1)[1]
		realDIF1=sum(abs(somaticres1$SMCN-somaticres1$SACN))/dim(somaticres1)[1]
		realsub1=sum(somaticres1$SACN==0)/dim(somaticres1)[1]
		if (realDIF1<realDIF | realsub1<realsub){
			out$segment=segment1
			out$somatic=somaticres1
			out$alpha=alpha1
			out$ploidy=ploidy1
		}
		randomsub1=apply(randomout1,2,subclone)
		PD1=1-length(aveDIF1[aveDIF1>=realDIF1])/length(aveDIF1)
		PS1=length(randomsub1[randomsub1<sum(somaticres1$SACN==0)])/times
		if (PD1 < 0.05 & PS1 < 0.05){
			out1$segment=segment1
			out1$somatic=somaticres1
			out1$alpha=alpha1
			out1$ploidy=ploidy1
			print.noquote(paste("iteration times = ",run,sep=""))
			p=PD1+PS1
			print.noquote(paste("p value = ",p,sep=""))
			return(out1)
			break
		}
		if ((PD1 + PS1) < (PD2+PS2)){
			out2$segment=segment1
			out2$somatic=somaticres1
			out2$alpha=alpha1
			out2$ploidy=ploidy1
			PD2=PD1
			PS2=PS1
		}
		realDIF=realDIF1
		realsub=realsub1
		somatic=somaticres1
		segment=segment1
		if(dim(somatic)[1]>=100){
			if (abs(alpha1-DNAalpha)<0.01){
				if (length(out2)!=0){
					print.noquote(paste("iteration times = ",run,sep=""))
					p=PD1+PS1
					print.noquote(paste("p value = ",p,sep=""))
					return(out2)
				}else{
					print.noquote(paste("iteration times = ",run,sep=""))
					p=PD1+PS1
					print.noquote(paste("p value = ",p,sep=""))
					return(out)
				}
			}else{
				DNAalpha=alpha1
				run=run+1
			}
		}else if (dim(somatic)[1]<100&dim(somatic)[1]>=20){
			if (abs(alpha1-DNAalpha)<0.001){
				if (length(out2)!=0){
					print.noquote(paste("iteration times = ",run,sep=""))
					p=PD1+PS1
					print.noquote(paste("p value = ",p,sep=""))
					return(out2)
				}else{
					print.noquote(paste("iteration times = ",run,sep=""))
					p=PD1+PS1
					print.noquote(paste("p value = ",p,sep=""))
					return(out)
				}
			}else{
				DNAalpha=alpha1
				run=run+1
			}
		}else{
			if (abs(alpha1-DNAalpha)<0.0001){
				if (length(out2)!=0){
					print.noquote(paste("iteration times = ",run,sep=""))
					p=PD1+PS1
					print.noquote(paste("p value = ",p,sep=""))
					return(out2)
				}else{
					print.noquote(paste("iteration times = ",run,sep=""))
					p=PD1+PS1
					print.noquote(paste("p value = ",p,sep=""))
					return(out)
				}
			}else{
				DNAalpha=alpha1
				run=run+1
			}
		}
	}
	if (run==cutoff){
		if (length(out2)!=0){
			print.noquote(paste("iteration times = ",run,sep=""))
			p=PD1+PS1
			print.noquote(paste("p value = ",p,sep=""))
			return(out2)
		}else{
			print.noquote(paste("iteration times = ",run,sep=""))
			p=PD1+PS1
			print.noquote(paste("p value = ",p,sep=""))
			return(out)
		}

	}
}
sequenzRun<-function(data,sample,chromosome){
	outdata=data.frame(chromosome=data$chr,position=data$pos,base.ref=data$ref,depth.normal=data$Nsum,depth.tumor=data$Tsum,depth.ratio=(data$Nsum)/(data$Tsum),Af=1-data$tfrac,Bf=data$tfrac)
	outdata$zygosity.normal[data$nfrac>=0.2&data$nfrac<=0.8]="het"
	outdata$zygosity.normal[data$nfrac<0.2|data$nfrac>0.8]="hom"
	outdata$GC.percent=50
	outdata$good.reads=outdata$depth.tumor
	genotype=data.frame(ref=data$ref[outdata$zygosity.normal=="het"],alt=data$alt[outdata$zygosity.normal=="het"])
	genotype=t(apply(genotype,1,sort))
	AB.normal=paste(genotype[,1],genotype[,2],sep="")
	outdata$AB.normal[outdata$zygosity.normal=="het"]=AB.normal
	genotype=data.frame(ref=data$ref[outdata$zygosity.normal=="hom"],alt=data$alt[outdata$zygosity.normal=="hom"])
	subdata=data[outdata$zygosity.normal=="hom",]
	if (dim(subdata)[1]>0){
		subdata$genotype[subdata$nfrac<0.2]=as.character(subdata$ref[subdata$nfrac<0.2])
		subdata$genotype[subdata$nfrac>0.8]=as.character(subdata$alt[subdata$nfrac>0.8])
		outdata$AB.normal[outdata$zygosity.normal=="hom"]=subdata$genotype
	}
	AB.tumor <- rep(".", length(outdata$AB.normal))
	outdata$AB.tumor=AB.tumor
	strand <- AB.tumor
	outdata$tumor.strand=strand
	normal.pos <- outdata$zygosity.normal == "hom" & outdata$AB.tumor == "."
	if (sum(normal.pos)>0){
		outdata <-outdata[outdata$depth.ratio > 0 & !is.infinite(outdata$depth.ratio) & !normal.pos, ]
	}
	outname=paste(sample,".seqz",sep="")
	write.table(outdata, outname, col.names = TRUE, row.names = FALSE, sep = "\t")
	seqz.data<-read.seqz(outname,gz=FALSE)
	##GC correction and normalization depth ratio
	print.noquote("GC correction")
	gc.stats<-gc.sample.stats(seqz.data)
	#gc.stats <- gc.norm(x = seqz.data$depth.ratio,gc = seqz.data$GC.percent)
	gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
	seqz.data$adjusted.ratio <- seqz.data$depth.ratio / gc.vect[as.character(seqz.data$GC.percent)]
	##extract information of sequenza input
	print.noquote("extract sequenza input")
	seqz.data$chr=seqz.data$chromosome
	test <- sequenza.extract(seqz.data,gc.stats,chroso=intersect(chromosome,seqz.data$chromosome))
	##infer tumor purity and ploidy
	print.noquote("sequenza infering")
	CP.example <- sequenza.fit(test,mc.cores = getOption("mc.cores", 1L),female=FALSE)
	cint <- get.ci(CP.example)
	alpha <- cint$max.cellularity
	ploidy <- cint$max.ploidy
	avg.depth.ratio <- mean(test$gc$adj[, 2])
	seg.tab <- na.exclude(do.call(rbind, test$segments))
	print.noquote("integrating result")
	cn.alleles <- baf.bayes(Bf = seg.tab$Bf, depth.ratio = seg.tab$depth.ratio,cellularity = alpha, ploidy = ploidy,avg.depth.ratio = avg.depth.ratio)
	seg.tab <- cbind(seg.tab, cn.alleles)
	segment=data.frame(chr=seg.tab$chromosome,startpos=seg.tab$start.pos,endpos=seg.tab$end.pos,nMajor=seg.tab$A,nMinor=seg.tab$B)
	index=which(segment$nMajor==0&segment$nMinor==0)
	rowindex=setdiff(c(1:dim(segment)[1]),index)
	segment=segment[rowindex,]
	DNAout=list(alpha=alpha,ploidy=ploidy,segment=segment)
	return(DNAout)
}

###DNArun: get final result at the DNA level
outputTitanSegments <- function(results, id, convergeParams, filename = NULL, igvfilename = NULL){
	# get all possible states in this set of results
	stateTable <- unique(results[, c("TITANstate", "TITANcall")])
	rownames(stateTable) <- stateTable[, 1]
	rleResults <- t(sapply(unique(results$Chr), function(x){
						ind <- results$Chr == x
						r <- rle(results$TITANstate[ind])
					}))
	rleLengths <- unlist(rleResults[, "lengths"])
	rleValues <- unlist(rleResults[, "values"])
	numSegs <- length(rleLengths)

	# convert allelic ratio to symmetric ratios #
	results$AllelicRatio <- apply(cbind(results$AllelicRatio, 1-results$AllelicRatio), 1, max, na.rm = TRUE)

	segs <- as.data.frame(matrix(NA, ncol = 14, nrow = numSegs,
					dimnames = list(c(), c("Sample", "Chromosome", "Start_Position.bp.", "End_Position.bp.",
									"Length.snp.", "Median_Ratio", "Median_logR", "TITAN_state", "TITAN_call", "Copy_Number",
									"MinorCN", "MajorCN", "Clonal_Cluster", "Cellular_Frequency"))))
	segs$Sample <- id
	colNames <- c("Chr", "Position", "TITANstate", "AllelicRatio", "LogRatio")
	prevInd <- 0
	for (j in 1:numSegs){
		start <- prevInd + 1
		end <- prevInd + rleLengths[j]
		segDF <- results[start:end, ]
		prevInd <- end
		numR <- nrow(segDF)
		segs[j, "Chromosome"] <- as.character(segDF[1, "Chr"])
		segs[j, "Start_Position.bp."] <- segDF[1, "Position"]
		segs[j, "TITAN_state"] <- rleValues[j]
		segs[j, "TITAN_call"] <- segDF[1, "TITANcall"]#stateTable[as.character(rleValues[j]), 2]
		segs[j, "Copy_Number"] <- segDF[1, "CopyNumber"]
		segs[j, "Median_Ratio"] <- round(median(segDF$AllelicRatio, na.rm = TRUE), digits = 6)
		segs[j, "Median_logR"] <- round(median(segDF$LogRatio, na.rm = TRUE), digits = 6)
		segs[j, "MinorCN"] <- getMajorMinorCN(rleValues[j], convergeParams$symmetric)$majorCN
		segs[j, "MajorCN"] <- getMajorMinorCN(rleValues[j], convergeParams$symmetric)$minorCN
		segs[j, "Clonal_Cluster"] <- segDF[1, "ClonalCluster"]
		segs[j, "Cellular_Frequency"] <- segDF[1, "CellularPrevalence"]
		if (segDF[1, "Chr"] == segDF[numR, "Chr"]){
			segs[j, "End_Position.bp."] <- segDF[numR, "Position"]
			segs[j, "Length.snp."] <- numR
		}else{ # segDF contains 2 different chromosomes
			print(j)
		}
	}
	if (!is.null(filename)){
		# write out detailed segment file #
		write.table(segs, file = filename, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	}
	# write out IGV seg file #
	if (!is.null(igvfilename)){
		igv <- segs[, c("Sample", "Chromosome", "Start_Position.bp.",
						"End_Position.bp.", "Length.snp.", "Median_logR")]
		colnames(igv) <- c("sample", "chr", "start", "end", "num.snps", "median.logR")
		write.table(igv, file = igvfilename, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	}
	return(segs)
}
getMajorMinorCN <- function(state, symmetric = TRUE){
	majorCN <- NA
	minorCN <- NA
	if (symmetric){
		if (state==0){
			majorCN = 0; minorCN = 0;
		}else if (state==1){
			majorCN = 0; minorCN = 1;
		}else if(state==2){
			majorCN = 0; minorCN = 2;
		}else if (state==3){
			majorCN = 1; minorCN = 1;
		}else if (state==4){
			majorCN = 0; minorCN = 3;
		}else if (state==5){
			majorCN = 1; minorCN = 2;
		}else if (state==6){
			majorCN = 0; minorCN = 4;
		}else if (state==7){
			majorCN = 1; minorCN = 3;
		}else if (state==8){
			majorCN = 2; minorCN = 2;
		}else if (state==9){
			majorCN = 0; minorCN = 5;
		}else if (state==10){
			majorCN = 1; minorCN = 4;
		}else if (state==11){
			majorCN = 2; minorCN = 3;
		}else if (state==12){
			majorCN = 0; minorCN = 6;
		}else if (state==13){
			majorCN = 1; minorCN = 5;
		}else if (state==14){
			majorCN = 2; minorCN = 4;
		}else if (state==15){
			majorCN = 3; minorCN = 3;
		}else if (state==16){
			majorCN = 0; minorCN = 7;
		}else if (state==17){
			majorCN = 1; minorCN = 6;
		}else if (state==18){
			majorCN = 2; minorCN = 5;
		}else if (state==19){
			majorCN = 3; minorCN = 4;
		}else if (state==20){
			majorCN = 0; minorCN = 8;
		}else if (state==21){
			majorCN = 1; minorCN = 7;
		}else if (state==22){
			majorCN = 2; minorCN = 6;
		}else if (state==23){
			majorCN = 3; minorCN = 5;
		}else if (state==24){
			majorCN = 4; minorCN = 4;
		}
	}else{
		#stop("symmetric=FALSE not yet supported.")
	}
	return(list(majorCN = majorCN, minorCN = minorCN))
}
TITANout<-function(DNAinput,chromosome){
	DNAinput=DNAinput[DNAinput$nfrac>=0.2&DNAinput$nfrac<=0.8,]
	rlength=nchar(as.character(DNAinput$ref))
	alength=nchar(as.character(DNAinput$alt))
	DNAinput=DNAinput[rlength==1&alength==1,]
	titandata=data.frame(chr=DNAinput$chr,pos=DNAinput$pos,ref=DNAinput$ref,refCount=DNAinput$refNumT,Nref=DNAinput$alt,NrefCount=DNAinput$altNumT)
	chro=intersect(chromosome,names(table(titandata$chr))[table(titandata$chr)>1])
	index=match(as.character(titandata$chr),chro)
	titandata=titandata[!is.na(index),]
	DNAinput=DNAinput[!is.na(index),]
	write.table(titandata,"TITAN.input",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
	numClusters <- 2
	params <- loadDefaultParameters(copyNumber=5,numberClonalClusters=numClusters)
	data <- loadAlleleCounts("TITAN.input")
	data$logR=log((DNAinput$refNumT+DNAinput$altNumT)/(DNAinput$refNumN+DNAinput$altNumN),base=10)
	convergeParams <- runEMclonalCN(data,gParams=params$genotypeParams,nParams=params$normalParams,
			pParams=params$ploidyParams,sParams=params$cellPrevParams,maxiter=20,maxiterUpdate=1500,
			txnExpLen=1e15,txnZstrength=1e5,useOutlierState=FALSE,normalEstimateMethod="map",
			estimateS=TRUE,estimatePloidy=TRUE)
	optimalPath <- viterbiClonalCN(data,convergeParams)
	results <- outputTitanResults(data,convergeParams,optimalPath,filename=NULL,posteriorProbs=F)
	results$AllelicRatio=as.numeric(results$AllelicRatio)
	results$LogRatio=as.numeric(results$LogRatio)
	segs <- outputTitanSegments(results, id = "test", convergeParams,filename=NULL)
	ploidy <- tail(convergeParams$phi, 1)
	normal <- tail(convergeParams$n, 1)
	#mean(as.numeric(segs$Cellular_Frequency[!is.na(segs$Cellular_Frequency)])*(1-normal))
	DNAout=list(segment=segs,ploidy=ploidy,alpha=1-normal)
	return(DNAout)
}


FACETSout<-function(DNAinput){
	facetsinput=data.frame(chr=DNAinput$chr,pos=DNAinput$pos,Nsum=DNAinput$Nsum,NAP=DNAinput$altNumN,Tsum=DNAinput$Tsum,TAP=DNAinput$altNumT)
	write.table(facetsinput,"FACETS.input",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
	set.seed(1234)
	xx=preProcSample(file="FACETS.input")
	oo=procSample(xx,cval=150)
	fit<-try(emcncf(oo),silent=TRUE)
	alpha=fit$purity
	if (!is.na(alpha)){
		ploidy=fit$ploidy
		segment=data.frame(chr=fit$cncf$chrom,startpos=fit$start,endpos=fit$end,nMajor=(fit$cncf$tcn.em-fit$cncf$lcn.em),nMinor=fit$cncf$lcn.em)
		segment=segment[!is.na(segment$nMajor)&!is.na(segment$nMinor)&segment$nMajor!=0&segment$nMinor!=0,]
		DNAout=list(alpha=alpha,segment=segment,ploidy=ploidy)
		return (DNAout)
	}
}
ASCATrun<-function(DNAinput,sample,chromosome){
	ASCATdata=ASCATin(data=DNAinput,sample=sample)
	##run ASCAT
	ASCATres<-try(ASCATout(ASCATdata=ASCATdata,sample=sample,chromosome=chromosome),silent=TRUE)
	if (mode(ASCATres)=="list"){
		return(ASCATres)
	}

}


DNACNA<-function(DNAinput,sample,tempth,chromosome,i){
	if (i==1){
		print.noquote("Run germline data")
		ASCATres<-try(ASCATrun(DNAinput,sample,chromosome),silent=TRUE)
		return(ASCATres)
	}else if (i==2){
		registerDoMC(cores = 4)
		TITANres<-try(TITANout(DNAinput=DNAinput,chromosome=chromosome),silent=TRUE)
		return (TITANres)
	}else if (i ==3){
		FACETSres<-try(FACETSout(DNAinput=DNAinput),silent=TRUE)
		return (FACETSres)
	}else{
		sequenzres<-try(sequenzRun(data=DNAinput,sample=sample,chromosome=chromosome),silent=TRUE)
		return (sequenzres)
	}
}

somaticCN<-function(res,i,somaticdata,method){
	if (mode(res[[i]])=="list"){
		if (method[i]=="TITAN" & !is.null(res[[i]])&!is.null(res[[i]]$alpha)){
			DNAalpha=res[[i]]$alpha
			DNAploidy=res[[i]]$ploidy
			segment=res[[i]]$segment
			segment=data.frame(chr=paste("chr",segment$Chromosome,sep=""),startpos=as.numeric(as.character(segment$Start_Position.bp.)),endpos=as.numeric(as.character(segment$End_Position.bp.)),nMajor=segment$MajorCN,nMinor=segment$MinorCN)
			somaticres=MCN(alpha=DNAalpha,segment=segment,somatic=somaticdata)
			return (somaticres)
		}else if (method[i]=="FACETS" & !is.null(res[[i]])&!is.null(res[[i]]$alpha)){
			DNAalpha=res[[i]]$alpha
			DNAploidy=res[[i]]$ploidy
			segment=res[[i]]$segment
			segment=data.frame(chr=paste("chr",segment$chr,sep=""),startpos=segment$startpos,endpos=segment$endpos,nMajor=segment$nMajor,nMinor=segment$nMinor)
			somaticres=MCN(alpha=DNAalpha,segment=segment,somatic=somaticdata)
			return (somaticres)
		}else if (!is.null(res[[i]])&!is.null(res[[i]]$alpha)){
			DNAalpha=res[[i]]$alpha
			DNAploidy=res[[i]]$ploidy
			segment=res[[i]]$segment
			somaticres=MCN(alpha=DNAalpha,segment=segment,somatic=somaticdata)
			return (somaticres)
		}
	}
}

randomCN<-function(res,i,randomdata,somaticres,method){
	if (mode(res[[i]])=="list"){
		if (method[i]=="TITAN" & !is.null(res[[i]])&!is.null(res[[i]]$alpha)&!is.null(somaticres[[i]])){
			if (dim(somaticres[[i]])[1]>0){
				DNAalpha=res[[i]]$alpha
				segment=res[[i]]$segment
				segment=data.frame(chr=paste("chr",segment$Chromosome,sep=""),startpos=as.numeric(as.character(segment$Start_Position.bp.)),endpos=as.numeric(as.character(segment$End_Position.bp.)),nMajor=segment$MajorCN,nMinor=segment$MinorCN)
				registerDoMC(cores = 4)
				randomout=foreach (j = 1:dim(randomdata[[i]])[2], .combine=cbind) %dopar% RMCN(R=randomdata[[i]][,j],somatic=somaticres[[i]],alpha=DNAalpha,segment=segment)
				#randomout=apply(randomdata[[i]],2,RMCN,somatic=somaticres[[i]],alpha=DNAalpha,segment=segment)
				return(randomout)
			}
		}else if (method[i]=="FACETS" & !is.null(res[[i]])&!is.null(res[[i]]$alpha)&!is.null(somaticres[[i]])){
			if (dim(somaticres[[i]])[1]>0){
				DNAalpha=res[[i]]$alpha
				segment=res[[i]]$segment
				segment=data.frame(chr=paste("chr",segment$chr,sep=""),startpos=segment$startpos,endpos=segment$endpos,nMajor=segment$nMajor,nMinor=segment$nMinor)
				registerDoMC(cores = 4)
				randomout=foreach (j = 1:dim(randomdata[[i]])[2], .combine=cbind) %dopar% RMCN(R=randomdata[[i]][,j],somatic=somaticres[[i]],alpha=DNAalpha,segment=segment)
				#randomout=apply(randomdata[[i]],2,RMCN,somatic=somaticres[[i]],alpha=DNAalpha,segment=segment)
				return(randomout)
			}
		}else if (!is.null(res[[i]])&!is.null(res[[i]]$alpha)&!is.null(somaticres[[i]])){
			if (dim(somaticres[[i]])[1]>0){
				DNAalpha=res[[i]]$alpha
				segment=res[[i]]$segment
				registerDoMC(cores = 4)
				randomout=foreach (j = 1:dim(randomdata[[i]])[2], .combine=cbind) %dopar% RMCN(R=randomdata[[i]][,j],somatic=somaticres[[i]],alpha=DNAalpha,segment=segment)
				return(randomout)
			}
		}
	}
}

randomDIF<-function(randomout,somaticres,i){
	if (!is.null(randomout[[i]])){
		if (dim(somaticres[[i]])>1){
			randomDIF=apply(randomout[[i]],2,difCN)
			aveDIF=apply(randomDIF,2,sum)/dim(somaticres[[i]])[1]
			return(aveDIF)
		}
	}
}
realDiff<-function(somaticres,i){
	if (!is.null(somaticres[[i]])){
		realDIF=sum(abs(somaticres[[i]]$SMCN-somaticres[[i]]$SACN))/dim(somaticres[[i]])[1]
		return(realDIF)
	}
}
realSub<-function(somaticres,i){
	if (!is.null(somaticres[[i]])){
		realsub=sum(somaticres[[i]]$SACN==0)/dim(somaticres[[i]])[1]
		return(realsub)
	}
}
subcloneP<-function(x){
	sub0<-function(y){
		return(sum(round(y)==0))
	}
	if (!is.null(x)){
		randomsub=apply(x,2,sub0)
		return(randomsub)
	}
}





DNArun<-function(SNPinput,somaticinput,sample,temppath){
	DNAinput=read.csv(SNPinput,sep="\t",header=TRUE)
	colnames(DNAinput)=c("chr","pos","ref","alt","refNumN","altNumN","refNumT","altNumT")
	chromosome=paste("chr",c(1:22,"X","Y"),sep="")
	if (length(levels(factor(DNAinput$chr)))>5){
		DNAinput$Tsum=DNAinput$refNumT+DNAinput$altNumT
		DNAinput$Nsum=DNAinput$refNumN+DNAinput$altNumN
		DNAinput=DNAinput[DNAinput$Tsum!=0,]
		DNAinput$tfrac=(DNAinput$altNumT)/(DNAinput$Tsum)
		DNAinput$nfrac=(DNAinput$altNumN)/(DNAinput$Nsum)
		DNAinput$TlogR=log2(DNAinput$Tsum/DNAinput$Nsum)
		DNAinput$BAF[DNAinput$tfrac>=0.5]=DNAinput$tfrac[DNAinput$tfrac>=0.5]
		DNAinput$BAF[DNAinput$tfrac<0.5]=1-DNAinput$tfrac[DNAinput$tfrac<0.5]
		DNAinput$Ratio=(DNAinput$altNumT+DNAinput$refNumT)/(DNAinput$altNumN+DNAinput$refNumN)
		DNAinput=DNAinput[DNAinput$Tsum>=10&DNAinput$Nsum>=10,]
		DNAinput=DNAinput[sapply(as.character(DNAinput$chr),nchar)<=5,]
		#rlength=nchar(as.character(DNAinput$ref))
		#alength=nchar(as.character(DNAinput$alt))
		#DNAinput=DNAinput[rlength==1&alength==1,]
		setwd(temppath)
		methods=c("ASCAT","TITAN","FACETS","sequenz")
		somaticdata=read.csv(somaticinput,header=TRUE,sep="\t")
		colnames(somaticdata)=c("chr","pos","ref","alt","refNumN","altNumN","refNumT","altNumT")
		somaticdata$Tsum=somaticdata$altNumT+somaticdata$refNumT
		somaticdata=somaticdata[somaticdata$Tsum>=10,]
		somaticdata$tfrac=somaticdata$altNumT/somaticdata$Tsum
		#rlength=nchar(as.character(somaticdata$ref))
		#alength=nchar(as.character(somaticdata$alt))
		#somaticdata=somaticdata[rlength==1&alength==1,]
		registerDoMC(cores = 4)
		times=1000
		res=foreach(i = 1:4) %dopar% DNACNA(DNAinput=DNAinput,sample=sample,tempth=temppath,chromosome=chromosome,i)
		print.noquote("Somatic mutation copy number")
		somaticres<-foreach(i =1:length(res)) %dopar% somaticCN(res=res,i,somaticdata=somaticdata,method=methods)
		somaticindex=0
		for (i in 1:length(somaticres)){
			if (!is.null(somaticres[[i]])){
				if (dim(somaticres[[i]])[1]>0){
					somaticindex=1
					break
				}
			}
		}
		if (somaticindex!=0){
			randomdata<-foreach(i=1:length(somaticres)) %dopar% randomSample(data=data.frame(Tsum=somaticres[[i]]$Tsum,tfrac=somaticres[[i]]$tfrac),times=times)
			randomout<-foreach(i =1:length(somaticres)) %dopar% randomCN(res,i,randomdata,somaticres,method=methods)
			aveDIF<-foreach(i=1:length(somaticres)) %dopar% randomDIF(randomout,somaticres,i)
			realDIF<-foreach(i=1:length(somaticres)) %dopar% realDiff(somaticres,i)
			realsub<-foreach(i=1:length(somaticres)) %dopar% realSub(somaticres,i)
			randomsub<-foreach(i = 1:length(somaticres)) %dopar% subcloneP(x=randomout[[i]])
			PDl=list()
			PSl=list()
			for (i in 1:length(realDIF)){
				if (!is.null(realDIF[[i]])){
					PDl[[i]]=1-length(aveDIF[[i]][aveDIF[[i]]>=realDIF[[i]]])/length(aveDIF[[i]])
					PSl[[i]]=length(randomsub[[i]][randomsub[[i]]<sum(somaticres[[i]]$SACN==0)])/times
				}
			}
			modelindex=c()
			k=1
			for (i in 1:length(PDl)){
				if (!is.null(PDl[[i]])){
					if(!is.na(PDl[[i]])){
						modelindex[k]=i
						k=k+1
					}
				}
			}
			PD=c()
			PS=c()
			for (k in 1:length(modelindex)){
				PD[k]=PDl[[modelindex[k]]]
				PS[k]=PSl[[modelindex[k]]]
			}
			minimum=min(PD+PS)
			selmeth=methods[modelindex[which(PD+PS==minimum)]]
			if ("ASCAT" %in% selmeth){
				DNAout=res[[1]]
				DNAalpha=res[[1]]$alpha
				DNAploidy=res[[1]]$ploidy
				segment=res[[1]]$segment
				somaticres=somaticres[[1]]
				realDIF=realDIF[[1]]
				realsub=realsub[[1]]
				randomdata=randomdata[[1]]
				PD=PDl[[1]]
				PS=PSl[[1]]
				resmthod="ASCAT"
			}else if ("FACETS" %in% selmeth){
				DNAout=res[[3]]
				DNAalpha=res[[3]]$alpha
				DNAploidy=res[[3]]$ploidy
				segment=res[[3]]$segment
				segment=data.frame(chr=paste("chr",segment$chr,sep=""),startpos=segment$startpos,endpos=segment$endpos,nMajor=segment$nMajor,nMinor=segment$nMinor)
				somaticres=somaticres[[3]]
				realDIF=realDIF[[3]]
				realsub=realsub[[3]]
				randomdata=randomdata[[3]]
				PD=PDl[[3]]
				PS=PSl[[3]]
				resmthod="FACETS"
			}else if ("TITAN" %in% selmeth){
				DNAout=res[[2]]
				DNAalpha=res[[2]]$alpha
				DNAploidy=res[[2]]$ploidy
				segment=res[[2]]$segment
				segment=data.frame(chr=paste("chr",segment$Chromosome,sep=""),startpos=as.numeric(as.character(segment$Start_Position.bp.)),endpos=as.numeric(as.character(segment$End_Position.bp.)),nMajor=segment$MajorCN,nMinor=segment$MinorCN)
				somaticres=somaticres[[2]]
				realDIF=realDIF[[2]]
				realsub=realsub[[2]]
				randomdata=randomdata[[2]]
				PD=PDl[[2]]
				PS=PSl[[2]]
				resmthod="TITAN"
			}else{
				DNAout=res[[4]]
				DNAalpha=res[[4]]$alpha
				DNAploidy=res[[4]]$ploidy
				segment=res[[4]]$segment
				somaticres=somaticres[[4]]
				realDIF=realDIF[[4]]
				realsub=realsub[[4]]
				randomdata=randomdata[[4]]
				PD=PDl[[4]]
				PS=PSl[[4]]
				resmthod="sequenz"
			}
			DNAout=list()
			DNAout$segment=segment
			DNAout$alpha=DNAalpha
			DNAout$ploidy=DNAploidy
			DNAout$somatic=somaticres
			DNAout$method=resmthod
			if (PD >= 0.05 | PS >= 0.05){
				print.noquote("Iter")
				iterout=iterOPT(SNP=DNAinput,somatic=somaticres,segment=segment,realDIF=realDIF,realsub=realsub,randomdata=randomdata,times=times,PD=PD,PS=PS,DNAalpha=DNAalpha,cutoff=50)
				if (length(iterout)==0){
					return(DNAout)
				}else{
					iterout$method=resmthod
					return(iterout)
				}
			}else{
				return(DNAout)
			}

		}else{
			DNAres=list()
			k=1
			for (i in 1:4){
				if (mode(res[[i]])=="list"){
					if (!is.null(res[[i]]$alpha)){
						res[[i]]$method=methods[i]
						DNAres[[k]]=res[[i]]
						k=k+1
					}
				}
			}
			ll=c()
			for (j in 1:length(DNAres)){
				ll[j]=dim(DNAres[[j]]$segment)[1]
			}
			DNAout=DNAres[[which.max(ll)]]
			return(DNAout)
		}
	}else{
			print.noquote(paste(sample,":Germline data is insuccifient",sep=""))
	}

}


DNArun1<-function(SNPinput,somaticinput,sample,temppath){
	DNAinput=read.csv(SNPinput,sep="\t",header=TRUE)
	colnames(DNAinput)=c("chr","pos","ref","alt","refNumN","altNumN","refNumT","altNumT")
	chromosome=paste("chr",c(1:22,"X","Y"),sep="")
	if (length(levels(factor(DNAinput$chr)))>5){
		DNAinput$Tsum=DNAinput$refNumT+DNAinput$altNumT
		DNAinput$Nsum=DNAinput$refNumN+DNAinput$altNumN
		DNAinput$tfrac=(DNAinput$altNumT)/(DNAinput$Tsum)
		DNAinput$nfrac=(DNAinput$altNumN)/(DNAinput$Nsum)
		DNAinput$TlogR=log2(DNAinput$Tsum/DNAinput$Nsum)
		DNAinput=DNAinput[DNAinput$Tsum!=0,]
		DNAinput$BAF[DNAinput$tfrac>=0.5]=DNAinput$tfrac[DNAinput$tfrac>=0.5]
		DNAinput$BAF[DNAinput$tfrac<0.5]=1-DNAinput$tfrac[DNAinput$tfrac<0.5]
		DNAinput$Ratio=(DNAinput$altNumT+DNAinput$refNumT)/(DNAinput$altNumN+DNAinput$refNumN)
		DNAinput=DNAinput[DNAinput$Tsum>=10&DNAinput$Nsum>=10,]
		DNAinput=DNAinput[sapply(as.character(DNAinput$chr),nchar)<=5,]
		#rlength=nchar(as.character(DNAinput$ref))
		#alength=nchar(as.character(DNAinput$alt))
		#DNAinput=DNAinput[rlength==1&alength==1,]
		setwd(temppath)
		methods=c("ASCAT","TITAN","FACETS","sequenz")
		somaticdata=read.csv(somaticinput,header=TRUE,sep="\t")
		colnames(somaticdata)=c("chr","pos","ref","alt","refNumN","altNumN","refNumT","altNumT")
		somaticdata$Tsum=somaticdata$altNumT+somaticdata$refNumT
		somaticdata$tfrac=somaticdata$altNumT/somaticdata$Tsum
		somaticdata=somaticdata[somaticdata$Tsum>=10,]
		#rlength=nchar(as.character(somaticdata$ref))
		#alength=nchar(as.character(somaticdata$alt))
		#somaticdata=somaticdata[rlength==1&alength==1,]
		registerDoMC(cores = 4)
		times=1000
		res=foreach(i = 1:4) %dopar% DNACNA(DNAinput=DNAinput,sample=sample,tempth=temppath,chromosome=chromosome,i)
		newmethod=c()
		newres=list()
		k=1
		DNAres=list()
		for (i in 1:4){
			if (!is.null(res[[i]])){
				if (!is.null(res[[i]]$alpha)){
					if (i ==1 ){
						DNAres$ASCAT=res[[i]]
					}else if (i ==2){
						DNAres$TITAN=res[[i]]
					}else if (i == 3){
						DNAres$FACETS=res[[i]]
					}else{
						DNAres$sequenza=res[[i]]
					}
					newmethod[k]=methods[i]
					newres[[k]]=res[[i]]
					k=k+1
				}
			}
		}
		print.noquote("Somatic mutation copy number")
		somaticres<-foreach(i =1:length(newmethod)) %dopar% somaticCN(res=newres,i,somaticdata=somaticdata,method=newmethod)
		somaticindex=0
		for (i in 1:4){
			if (!is.null(somaticres[[i]])){
				if (dim(somaticres[[i]])[1]>0){
					somaticindex=1
					break
				}
			}
		}
		if (somaticindex!=0){
			randomdata<-foreach(i=1:length(newmethod)) %dopar% randomSample(data=data.frame(Tsum=somaticres[[i]]$Tsum,tfrac=somaticres[[i]]$tfrac),times=times)
			randomout<-foreach(i =1:length(newmethod)) %dopar% randomCN(newres,i,randomdata,somaticres,method=newmethod)
			aveDIF<-foreach(i=1:length(newmethod)) %dopar% randomDIF(randomout,somaticres,i)
			realDIF<-foreach(i=1:length(newmethod)) %dopar% realDiff(somaticres,i)
			realsub<-foreach(i=1:length(newmethod)) %dopar% realSub(somaticres,i)
			randomsub<-foreach(i = 1:length(newmethod)) %dopar% subcloneP(x=randomout[[i]])
			PDl=list()
			PSl=list()
			for (i in 1:length(newmethod)){
				if (!is.null(realDIF[[i]])){
					PDl[[i]]=1-length(aveDIF[[i]][aveDIF[[i]]>=realDIF[[i]]])/length(aveDIF[[i]])
					PSl[[i]]=length(randomsub[[i]][randomsub[[i]]<sum(somaticres[[i]]$SACN==0)])/times
				}
			}
			modelindex=c()
			k=1
			for (i in 1:length(newmethod)){
				if (!is.null(PDl[[i]])){
					if (!is.na(PDl[[i]])){
						modelindex[k]=i
						k=k+1
					}
				}
			}
			PD=c()
			PS=c()
			for (k in 1:length(modelindex)){
				PD[k]=PDl[[modelindex[k]]]
				PS[k]=PSl[[modelindex[k]]]
			}
			minimum=min(PD+PS)
			selmeth=newmethod[modelindex[which(PD+PS==minimum)]]
			if ("ASCAT" %in% selmeth){
				DNAout=newres[[which(newmethod=="ASCAT")]]
				DNAalpha=newres[[which(newmethod=="ASCAT")]]$alpha
				DNAploidy=newres[[which(newmethod=="ASCAT")]]$ploidy
				segment=newres[[which(newmethod=="ASCAT")]]$segment
				somaticres=somaticres[[which(newmethod=="ASCAT")]]
				realDIF=realDIF[[which(newmethod=="ASCAT")]]
				realsub=realsub[[which(newmethod=="ASCAT")]]
				randomdata=randomdata[[which(newmethod=="ASCAT")]]
				PD=PDl[[which(newmethod=="ASCAT")]]
				PS=PSl[[which(newmethod=="ASCAT")]]
				resmthod="ASCAT"
			}else if ("FACETS" %in% selmeth){
				DNAout=newres[[which(newmethod=="FACETS")]]
				DNAalpha=newres[[which(newmethod=="FACETS")]]$alpha
				DNAploidy=newres[[which(newmethod=="FACETS")]]$ploidy
				segment=newres[[which(newmethod=="FACETS")]]$segment
				segment=data.frame(chr=paste("chr",segment$chr,sep=""),startpos=segment$startpos,endpos=segment$endpos,nMajor=segment$nMajor,nMinor=segment$nMinor)
				somaticres=somaticres[[which(newmethod=="FACETS")]]
				realDIF=realDIF[[which(newmethod=="FACETS")]]
				realsub=realsub[[which(newmethod=="FACETS")]]
				randomdata=randomdata[[which(newmethod=="FACETS")]]
				PD=PDl[[which(newmethod=="FACETS")]]
				PS=PSl[[which(newmethod=="FACETS")]]
				resmthod="FACETS"
			}else if ("TITAN" %in% selmeth){
				DNAout=newres[[which(newmethod=="TITAN")]]
				DNAalpha=newres[[which(newmethod=="TITAN")]]$alpha
				DNAploidy=newres[[which(newmethod=="TITAN")]]$ploidy
				segment=newres[[which(newmethod=="TITAN")]]$segment
				segment=data.frame(chr=paste("chr",segment$Chromosome,sep=""),startpos=as.numeric(as.character(segment$Start_Position.bp.)),endpos=as.numeric(as.character(segment$End_Position.bp.)),nMajor=segment$MajorCN,nMinor=segment$MinorCN)
				somaticres=somaticres[[which(newmethod=="TITAN")]]
				realDIF=realDIF[[which(newmethod=="TITAN")]]
				realsub=realsub[[which(newmethod=="TITAN")]]
				randomdata=randomdata[[which(newmethod=="TITAN")]]
				PD=PDl[[which(newmethod=="TITAN")]]
				PS=PSl[[which(newmethod=="TITAN")]]
				resmthod="TITAN"
			}else{
				DNAout=newres[[which(newmethod=="sequenz")]]
				DNAalpha=newres[[which(newmethod=="sequenz")]]$alpha
				DNAploidy=newres[[which(newmethod=="sequenz")]]$ploidy
				segment=newres[[which(newmethod=="sequenz")]]$segment
				somaticres=somaticres[[which(newmethod=="sequenz")]]
				realDIF=realDIF[[which(newmethod=="sequenz")]]
				realsub=realsub[[which(newmethod=="sequenz")]]
				randomdata=randomdata[[which(newmethod=="sequenz")]]
				PD=PDl[[which(newmethod=="sequenz")]]
				PS=PSl[[which(newmethod=="sequenz")]]
				resmthod="sequenz"
			}
			DNAout=list()
			DNAout$segment=segment
			DNAout$alpha=DNAalpha
			DNAout$ploidy=DNAploidy
			DNAout$somatic=somaticres
			DNAout$method=resmthod
			#DNAres$DNAout=DNAout
			return(DNAout)
		}else{
			DNAres=list()
			k=1
			for (i in 1:length(res)){
				if (mode(res[[i]])=="list"){
					if (!is.null(res[[i]]$alpha)){
						res[[i]]$method=methods[i]
						DNAres[[k]]=res[[i]]
						k=k+1
					}
				}
			}
			ll=c()
			for (j in 1:length(DNAres)){
				ll[j]=dim(DNAres[[j]]$segment)[1]
			}
			DNAout=DNAres[[which.max(ll)]]
			return(DNAout)
		}
	}else{
			print.noquote(paste(sample,":Germline data is insuccifient",sep=""))
	}

}

Consistent<-function(GCN,MCN){
	lower=min(c(GCN,MCN))-1
	upper=max(c(GCN,MCN))+1
	da=density(GCN,from=lower,to=upper)
	db=density(MCN,from=lower,to=upper)
	d=data.frame(x=da$x,a=da$y,b=db$y)
	d$w=pmin(d$a,d$b)
	total=integrate.xy(d$x,d$a)+integrate.xy(d$x,d$b)
	intersection=integrate.xy(d$x,d$w)
	overlap=2*intersection/total
	diff=(integrate.xy(d$x,d$b)-intersection)/integrate.xy(d$x,d$b)
	return(list(overlap=overlap,diff=diff))
}

Heterogeneity<-function(DNAout){
	if ("somatic" %in% names(DNAout)){
		segment=DNAout$segment
		somatic=DNAout$somatic
		somatic1=c()
		for (i in 1:dim(segment)[1]){
			subdata=somatic[somatic$chr==as.character(segment$chr[i])&somatic$pos>=segment$startpos[i]&somatic$pos<=segment$endpos[i],]
			if (dim(subdata)[1]>0){
				subdata$Dmajor=segment$nMajor[i]
				subdata$Dminor=segment$nMinor[i]
				somatic1=rbind(somatic1,subdata)
			}
		}
		DNAout$heterogeneity=1-Consistent(c(somatic1$Dmajor,somatic1$Dminor),somatic1$SMCN)$overlap
	}else{
		DNAout$heterogeneity=0
	}
	return(DNAout)
}


###############################################
wildtype<-function(data,segment,type,resout){
  for (i in 1:dim(segment)[1]){
    subdata=data[data$chr==as.character(segment$chr[i])&data$pos>=segment$startpos[i]&data$pos<=segment$endpos[i],]
    Dmajor=segment$nMajor[i]
    Dminor=segment$nMinor[i]
    if (dim(subdata)[1]>0){
    	subres=subdata[,1:6]
    	if (type == "somatic"){
    		index=which(is.na(subdata$BayesP))
    		if (length(index)>0){
	    		subdata$RTCN[index]=(Dmajor+Dminor)*2^(subdata$ratio[index])
	    		subdata$RMCN[index]=subdata$RTCN[index]*subdata$tfrac[index]
	    		subdata$BayesP[index]=1-10^(-abs(subdata$RMCN[index]-subdata$DMCN[index]))
    		}
    		index=which(subdata$RTCN==Inf|subdata$RTCN==-Inf)
    		if (length(index)>0){
    			subdata$RTCN[index]=(Dmajor+Dminor)*2^(subdata$ratio[index])
    			subdata$RMCN[index]=subdata$RTCN[index]*subdata$tfrac[index]
    		}
        }
        if (type=="germline"){
        	altD=subdata$Dminor
        	wildD=subdata$Dmajor
        	altD[subdata$alt==subdata$DMajorAllele]=subdata$Dmajor[subdata$alt==subdata$DMajorAllele]
        	wildD[subdata$alt==subdata$DMajorAllele]=subdata$Dminor[subdata$alt==subdata$DMajorAllele]
            altR=abs(subdata$Rminor)
            wildR=subdata$Rmajor
            altR[subdata$alt==subdata$RMajorAllele]=subdata$Rmajor[subdata$alt==subdata$RMajorAllele]
            wildR[subdata$alt==subdata$RMajorAllele]=abs(subdata$Rminor[subdata$alt==subdata$RMajorAllele])
            subres=cbind(subres,rep("Germline",length=dim(subres)[1]))
            subres=cbind(subres,altD,altR,wildD,wildR)
            subres=cbind(subres,subdata$BayesP)
        }else{
            altD=subdata$DMCN
            altR=subdata$RMCN
            if (abs(altD-Dmajor)<abs(altD-Dminor)){
              wildD=Dminor
            }else{
              wildD=Dmajor
            }
            wildR=subdata$RTCN-subdata$RMCN
            subres=cbind(subres,rep("Somatic",length=dim(subres)))
            subres=cbind(subres,altD,altR,wildD,wildR)
            subres=cbind(subres,subdata$BayesP)
        }
        names(subres)=c("chr","pos","ref","alt","refNum","altNum","type","altD","altR","wildD","wildR","BayesP")
        resout=rbind(resout,subres)
    }
  }
  return(resout)
}

altalle<-function(x){
  dif1=abs(x[2]-x[1])
  dif2=abs(x[2]-x[3])
  if (dif1<=dif2){
    return(data.frame(alt=x[1],ref=x[3]))
  }else{
    return(data.frame(alt=x[3],ref=x[1]))
  }
}

WGD<-function(segment,SNPinput){
	DNAinput=read.csv(SNPinput,sep="\t",header=TRUE)
	colnames(DNAinput)=c("chr","pos","ref","alt","refNumN","altNumN","refNumT","altNumT")
	segment$TCN=segment$Dmajor+segment$Dminor
	if (dim(segment)[1]==1){
		subdata=DNAinput[DNAinput$chr==as.character(segment$chr)&DNAinput$pos>=segment$start&DNAinput$pos<=segment$end,]
		ll=dim(DNAinput)[1]
		gdouble=dim(subdata)[1]
	}else{
		ll=sum(as.numeric(segment$end-segment$start))
		gdouble=sum(as.numeric(segment$end[segment$TCN>2]-segment$start[segment$TCN>2]))
	}

	return(gdouble/ll)

}
DACRE<-function(resout){
	germ=resout[resout$type=="Germline",]
	germ1=do.call(rbind,apply(germ[,8:10],1,altalle))
	germ1=cbind(germ$altR,germ$wildR,germ1,germ$BayesP)
	germ1=as.data.frame(germ1)
	names(germ1)=c("ASELm","ASELw","ASCNm","ASCNw","BayesP")
	germ1$ASEmP=germ1$ASELm-exp(-abs(germ1$ASCNm-1))*germ1$ASCNm
	germ1$ASEwP=germ1$ASELw-exp(-abs(germ1$ASCNw-1))*germ1$ASCNw
	germ1$ASEm=germ1$ASELm-germ1$ASCNm
	germ1$ASEw=germ1$ASELw-germ1$ASCNw
	germ1$ASELm[germ1$ASELm==0]=0.01
	germ1$ASELw[germ1$ASEKw==0]=0.01
	germ1$AEI=log(germ1$ASELm/germ1$ASELw,base=2)
	germ1$DACRE=germ1$ASEmP-germ1$ASEwP
	germres=cbind(germ,germ1$ASEm,germ1$AEI,germ1$DACRE)
	names(germres)=c(names(germ),"eASEL","AEI","DACRE")
	soma=resout[resout$type=="Somatic",]
	if (dim(soma)[1]>0){
		soma$wildR[soma$wildR<0]=soma$altR[soma$wildR<0]*soma$wildD[soma$wildR<0]/soma$altD[soma$wildR<0]
		soma1=soma
		soma1$ASEmP=soma1$altR-exp(-abs(soma1$altD-1))*soma1$altD
		soma1$ASEwP=soma1$wildR-exp(-abs(soma1$wildD-1))*soma1$wildD
		soma1$ASEm=soma1$altR-soma1$altD
		soma1$ASEw=soma1$wildR-soma1$wildD
		soma1$wildR[soma1$wildR==0]=0.01
		soma1$altR[soma1$altR==0]=0.01
		soma1$AEI=log(soma1$altR/soma1$wildR,base=2)
		soma1$DACRE=soma1$ASEmP-soma1$ASEwP
		somares=cbind(soma,soma1$ASEm,soma1$AEI,soma1$DACRE)
		names(somares)=c(names(soma),"eASEL","AEI","DACRE")
		res=rbind(germres,somares)
	}else{
		res=germres
	}
	return(res)
}
