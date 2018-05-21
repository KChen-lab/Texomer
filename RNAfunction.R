# TODO: Add comment
# 
# Author: FWang9
###############################################################################


MajorAlle<-function(nMajor,nMinor,germline){
	germline$RMajorAllele[germline$tfrac>=0.5]=as.character(germline$alt[germline$tfrac>=0.5])
	germline$RMajorAllele[germline$tfrac<0.5]=as.character(germline$ref[germline$tfrac<0.5])
	if (nMajor!=nMinor){
		germline$DMajorAllele[germline$Dfrac>=0.5]=as.character(germline$alt[germline$Dfrac>=0.5])
		germline$DMajorAllele[germline$Dfrac<0.5]=as.character(germline$ref[germline$Dfrac<0.5])
	}else{
		germline$DMajorAllele=germline$RMajorAllele
	}
	germline$y[germline$tfrac>=0.5]=germline$altNum[germline$tfrac>=0.5]
	germline$y[germline$tfrac<0.5]=germline$refNum[germline$tfrac<0.5]
	return(germline)
}
mtmp<-function(data,nMajor,TCN,alpha,NR,theta,pai){
	f1=(1-alpha+alpha*nMajor)/(2*(1-alpha)+alpha*TCN)
	f2=(1-alpha+alpha*NR)/(2*(1-alpha)+alpha*TCN)
	p1=pai*dbetabinom(data$y,f1,data$N,theta,log=FALSE)
	p2=(1-pai)*dbetabinom(data$y,f2,data$N,theta,log=FALSE)
	-sum(log(p1+p2))
}
bbml<-function(model,nMajor,TCN,germline,Dpurity,Rpurity){
		if (Rpurity=="NA"){
			like=matrix(0,nr=1,nc=5)
			colnames(like)=c("likelihood","alpha","Rmajor","pai","theta")
			m0 <- try(mle2(mtmp,start=list(NR=nMajor,theta=9,pai=0.5,alpha=Dpurity),data=list(data=germline,nMajor=nMajor,TCN=TCN),method=model),silent=TRUE)
			if (mode(m0)=="S4"){
				like[1,1]=m0@min
				parameter=m0@coef
				like[1,2]=parameter[1]
				like[1,3]=parameter[2]
				like[1,4]=parameter[4]
				like[1,5]=parameter[3]
			}
		}else{
			like=matrix(0,nr=1,nc=4)
			colnames(like)=c("likelihood","Rmajor","pai","theta")
			m0 <- try(mle2(mtmp,start=list(NR=nMajor,theta=9,pai=0.5),data=list(data=germline,nMajor=nMajor,TCN=TCN,alpha=Rpurity),method=model),silent=TRUE)
			if (mode(m0)=="S4"){
				like[1,1]=m0@min
				parameter=m0@coef
				like[1,2]=parameter[1]
				like[1,3]=parameter[3]
				like[1,4]=parameter[2]
			}
		}
		return(like)
}
Like<-function(nMajor,TCN,germline,Dpurity,Rpurity){
	method=c("Nelder-Mead","BFGS", "CG","SANN")
	like=foreach (j = 1:length(method), .combine=rbind) %dopar% bbml(model=method[j],nMajor,TCN,germline,Dpurity,Rpurity)
	row.names(like)=method
	return(like)
}
IteLike<-function(nMajor,TCN,germline){
	paiti=seq(from=0.1,to=0.9,by=0.1)
	alphati=seq(from=0.1,to=0.9,by=0.1)
	like1=matrix(0,nr=length(paiti)*length(alphati),nc=5)
	like1=t(like1)
	row.names(like1)=c("likelihood","alpha","Rmajor","pai","theta")
	like1=t(like1)
	like1=as.data.frame(like1)
	like2=like1
	s=1
	for (j in 1:length(alphati)){
		for (k in 1:length(paiti)){
			try(m0 <- mle2(mtmp,start=list(NR=nMajor,theta=9),data=list(data=germline,nMajor=nMajor,TCN=TCN,alpha=alphati[j],pai=paiti[k]),method="Nelder-Mead"))
			parameter=m0@coef
			like1[s,1]=m0@min
			like1[s,2]=alphati[j]
			like1[s,3]=parameter[1]
			like1[s,4]=paiti[k]
			like1[s,5]=parameter[2]
			s=s+1
		}
	}
	like1=like1[apply(like1,1,sum)!=0,]
	return(like1[which.min(like1[,1]),])
}
segalpha<-function(segment,RNA,DNA,DNApurity,i,Rpurity){
		nMajor=segment$nMajor[i]
		nMinor=segment$nMinor[i]
		TCN=nMajor+nMinor
		subdata=RNA[RNA$chr==as.character(segment$chr[i])&RNA$pos>=segment$startpos[i]&RNA$pos<=segment$endpos[i],]
		subdata=subdata[subdata$altNum+subdata$refNum>=10,]
		germline=subdata[subdata$type=="germline"&!is.na(subdata$tfrac),]
		subDNA=DNA[DNA$chr==as.character(segment$chr[i])&DNA$pos>=segment$startpos[i]&DNA$pos<=segment$endpos[i],]
		germline=germline[!is.na(match(germline$pos,subDNA$pos)),]
		if (dim(germline)[1] >= 1){
			germline$tfrac=germline$altNum/(germline$refNum+germline$altNum)
			germline$Dfrac=subDNA$altNumT[match(germline$pos,subDNA$pos)]/(subDNA$altNumT[match(germline$pos,subDNA$pos)]+subDNA$refNumT[match(germline$pos,subDNA$pos)])
			germline$N=germline$refNum+germline$altNum
			germline$Dmajor=nMajor
			germline$Dminor=nMinor
			germline=MajorAlle(nMajor=nMajor,nMinor=nMinor,germline=germline)
			like=Like(nMajor=nMajor,TCN=TCN,germline=germline,Dpurity=DNApurity,Rpurity=Rpurity)
			like=as.data.frame(like)
			like=like[like$likelihood!=0,]
			like=like[order(like$likelihood,decreasing=TRUE),]
			for (j in 1:dim(like)[1]){
				if (like$alpha[j]>=0 & like$alpha[j] <= 1 & like$pai[j] >= 0 & like$pai[j] <= 1){
					RNAalpha=like$alpha[j]
					RNARmajor=like$Rmajor[j]
					theta=like$theta[j]
					pai=like$pai[j]
					likelihood=like$likelihood[j]
					break
				}
			}
			if ("RNAalpha" %in% ls()){
				if (abs(RNAalpha-DNApurity)>0.5){
					itelike=IteLike(nMajor=nMajor,TCN=TCN,germline=germline)
					itelikelihood=itelike$likelihood
					if (itelikelihood>likelihood){
						RNAalpha=itelike$alpha
						RNARmajor=itelike$Rmajor
						theta=itelike$theta
						pai=itelike$pai
					}
				}
			}else {
				itelike=IteLike(nMajor=nMajor,TCN=TCN,germline=germline)
				if (sum(itelike)!=0){
					RNAalpha=itelike$alpha
					RNARmajor=itelike$Rmajor
					theta=itelike$theta
					pai=itelike$pai
				}else{
					RNAalpha=like$alpha[which.max(like$likelihood)]
					RNARmajor=like$Rmajor[which.max(like$likelihood)]
					theta=like$theta[which.max(like$likelihood)]
					pai=like$pai[which.max(like$likelihood)]
				}
				
			}
			out=matrix(nr=1,nc=2)
			colnames(out)=c("purity","count")
			out[1,1]=RNAalpha
			out[1,2]=dim(germline)[1]
			return(out)
		}
}
SNPclass<-function(germline,TCN,alpha,nMajor){
	#germline$Rmajor=(germline$y/germline$N*(2*(1-alpha)+alpha*TCN)-(1-alpha))/alpha
	#germline$Rminor=((1-germline$y/germline$N)*(2*(1-alpha)+alpha*TCN)-(1-alpha))/alpha
	Rmajor1=try(normalmixEM(germline$Rmajor, lambda = 0.5, mu = c(min(germline$Rmajor), max(germline$Rmajor)), sigma = sd(germline$Rmajor)),silent=TRUE)
	if (!is.null(names(Rmajor1))){
		res=matrix(nr=3,nc=2)
		row.names(res)=c("mu","pai","sd")
		res[1,]=Rmajor1$mu
		res[2,]=Rmajor1$lambda
		res[3,]=Rmajor1$sigma
		f<-res[2,1]*dnorm(germline$Rmajor,res[1,1],res[3,1])+res[2,2]*dnorm(germline$Rmajor,res[1,2],res[3,2])
		k=which.max(res[2,])
		germline$BayesP=1-res[2,k]*dnorm(germline$Rmajor,res[1,k],res[3,k])/f
	}else{
		f1=dnorm(germline$Rmajor,nMajor,sd(germline$Rmajor))
		f2=dnorm(germline$Rmajor,round(median(germline$Rmajor)),sd(germline$Rmajor))
		germline$BayesP=f2/(f1+f2)
	}
	return(germline)
}
SNPclass2<-function(germline,DTCN,RTCN,alpha,nMajor){
	#germline$Rmajor=(germline$y/germline$N*(2*(1-alpha)+alpha*RTCN)-(1-alpha))/alpha
	#germline$Rminor=((1-germline$y/germline$N)*(2*(1-alpha)+alpha*RTCN)-(1-alpha))/alpha
	f1=dnorm(germline$Rmajor,nMajor,sd(germline$Rmajor))
	f2=dnorm(germline$Rmajor,round(median(germline$Rmajor)),sd(germline$Rmajor))
	germline$BayesP=f2/(f1+f2)
	return(germline)
}
SNVclass<-function(Rsomatic){
	somaticSNV=normalmixEM(Rsomatic$delt, lambda = rep(1, 3)/3, mu = c(min(Rsomatic$delt),0, max(Rsomatic$delt)), sigma = sd(Rsomatic$delt))
	res=matrix(nr=3,nc=3)
	row.names(res)=c("mu","pai","sd")
	res[1,]=somaticSNV$mu
	res[2,]=somaticSNV$lambda
	res[3,]=somaticSNV$sigma
	f<-res[2,1]*dnorm(Rsomatic$delt,res[1,1],res[3,1])+res[2,2]*dnorm(Rsomatic$delt,res[1,2],res[3,2])+res[2,3]*dnorm(Rsomatic$delt,res[1,3],res[3,3])
	k=which.min(abs(res[1,]))
	Rsomatic$BayesP=1-res[2,k]*dnorm(Rsomatic$delt,res[1,k],res[3,k])/f
	return(Rsomatic)
}
SNVclass2<-function(Rsomatic){
	DSNV=normalmixEM(Rsomatic$Ddelt, lambda = 0.5, mu = c(min(Rsomatic$Ddelt), max(Rsomatic$Ddelt)), sigma = sd(Rsomatic$Ddelt))
	resD=matrix(nr=3,nc=2)
	row.names(resD)=c("mu","pai","sd")
	resD[1,]=DSNV$mu
	resD[2,]=DSNV$lambda
	resD[3,]=DSNV$sigma
	fD<-resD[2,1]*dnorm(Rsomatic$delt,resD[1,1],resD[3,1])+resD[2,2]*dnorm(Rsomatic$delt,resD[1,2],resD[3,2])
	RSNV=normalmixEM(Rsomatic$delt, lambda = rep(1, 3)/3, mu = c(min(Rsomatic$delt),0, max(Rsomatic$delt)), sigma = sd(Rsomatic$delt))
	resR=matrix(nr=3,nc=3)
	row.names(resR)=c("mu","pai","sd")
	resR[1,]=RSNV$mu
	resR[2,]=RSNV$lambda
	resR[3,]=RSNV$sigma
	fR<-resR[2,1]*dnorm(Rsomatic$delt,resR[1,1],resR[3,1])+resR[2,2]*dnorm(Rsomatic$delt,resR[1,2],resR[3,2])+resR[2,3]*dnorm(Rsomatic$delt,resR[1,3],resR[3,3])
	Rsomatic$BayesP=fR/(fR+fD)
	return(Rsomatic)
}
RNApurity<-function(Rpurity,Count){
	alpha=c()
	Count=Count[!is.na(Rpurity)]
	Rpurity=Rpurity[!is.na(Rpurity)]
	for (i in 1:length(Count)){
		weight=rep(Rpurity[i],length=Count[i])
		alpha=c(alpha,weight)
	}
	Ralpha=density(alpha)$x[which.max(density(alpha)$y)]
	return(Ralpha)
	
}
#SNPASE<-function(segment,DNA,RNA,DNApurity,Rpurity,ratioS,somatic){
#	for (i in 1:dim(segment)[1]){
#		nMajor=segment$nMajor[i]
#		nMinor=segment$nMinor[i]
#		TCN=nMajor+nMinor
#		subdata=RNA[RNA$chr==as.character(segment$chr[i])&RNA$pos>=segment$startpos[i]&RNA$pos<=segment$endpos[i],]
#		subdata=subdata[subdata$altNum+subdata$refNum>=10,]
#		germline=subdata[subdata$type=="germline"&!is.na(subdata$tfrac),]
#		subDNA=DNA[DNA$chr==as.character(segment$chr[i])&DNA$pos>=segment$startpos[i]&DNA$pos<=segment$endpos[i],]
#		index=match(germline$pos,subDNA$pos)
#		germline=germline[!is.na(index),]
#		germline$ND=subDNA$refNumT[index[!is.na(index)]]+subDNA$altNumT[index[!is.na(index)]]
#		Rsomatic=subdata[subdata$type=="somatic"&!is.na(subdata$tfrac),]
#		index=match(paste(Rsomatic$chr,Rsomatic$pos),paste(somatic$chr,somatic$pos))
#		Rsomatic=Rsomatic[!is.na(index),]
#		Rsomatic$DMCN=somatic$SMCN[index[!is.na(index)]]
#		Rsomatic$Dfrac=somatic$tfrac[index[!is.na(index)]]
#		if (dim(germline)[1]>0){
#			germline$tfrac=germline$altNum/(germline$refNum+germline$altNum)
#			germline$Dfrac=subDNA$altNumT[match(germline$pos,subDNA$pos)]/(subDNA$altNumT[match(germline$pos,subDNA$pos)]+subDNA$refNumT[match(germline$pos,subDNA$pos)])
#			germline$N=germline$refNum+germline$altNum
#			germline$Dmajor=nMajor
#			germline$Dminor=nMinor
#			germline=MajorAlle(nMajor=nMajor,nMinor=nMinor,germline=germline)
#			if (segment$BayesP[i] <= 0.5){
#				RTCN=TCN
#				segment$RTCN[i]=RTCN
#				segment$RMajor[i]=nMajor
#				segment$RMinor[i]=nMinor
#			}else{
#				TCNRNA=(germline$N/germline$ND/ratioS*(2*(1-DNApurity)+DNApurity*TCN)-2*(1-Rpurity))/Rpurity
#				if (length(TCNRNA[TCNRNA>0])==0){
#					RTCN=round(median(-TCNRNA))
#				}else{
#					RTCN=round(median(TCNRNA[TCNRNA>0]))
#				}
#				RMajor=(germline$y/germline$N*(2*(1-Rpurity)+Rpurity*RTCN)-(1-Rpurity))/Rpurity
#				RMinor=((1-germline$y/germline$N)*(2*(1-Rpurity)+Rpurity*RTCN)-(1-Rpurity))/Rpurity
#				if (round(median(RMinor))<0){
#					segment$RTCN[i]=RTCN+abs(round(median(RMinor)))
#					segment$RMajor[i]=round(median(RMajor))+abs(round(median(RMinor)))
#					segment$RMinor[i]=abs(round(median(RMinor)))
#				}else{
#					segment$RTCN[i]=RTCN
#					segment$RMajor[i]=round(median(RMajor))
#					segment$RMinor[i]=round(median(RMinor))
#				}
#			}
#			if (dim(germline)[1]>=20){
#				if (segment$BayesP[i] <= 0.5){
#					germline=SNPclass(germline=germline,TCN=RTCN,alpha=Rpurity,nMajor=nMajor)
#				}else{
#					germline=SNPclass2(germline,DTCN=TCN,RTCN=RTCN,alpha=Rpurity,nMajor=nMajor)
#				}
#			}else{
#				germline$Rmajor=(germline$y/germline$N*(2*(1-Rpurity)+Rpurity*TCN)-(1-Rpurity))/Rpurity
#				germline$Rminor=((1-germline$y/germline$N)*(2*(1-Rpurity)+Rpurity*TCN)-(1-Rpurity))/Rpurity
#				germline$BayesP=1-10^(-abs(round(germline$Rmajor)-nMajor))
#			}
#			if (!("Rgermlineres" %in% ls())){
#				Rgermlineres = germline
#			}else{
#				Rgermlineres=rbind(Rgermlineres,germline)
#			}
#	
#		}
#		if (dim(Rsomatic)[1]>=1){
#			Rsomatic$RMCN=Rsomatic$tfrac*(2*(1-Rpurity)+Rpurity*RTCN)/Rpurity
#			Rsomatic$Ddelt=Rsomatic$DMCN-round(Rsomatic$DMCN)
#			Rsomatic$delt=Rsomatic$RMCN-Rsomatic$DMCN
#			if (dim(Rsomatic)[1]>=10){
#				if (RTCN==TCN){
#					Rsomatic=SNVclass(Rsomatic)
#				}else{
#					Rsomatic=SNVclass2(Rsomatic)
#				}
#			}else{
#				Rsomatic$BayesP=1-10^(-abs(round(Rsomatic$RMCN)-round(Rsomatic$DMCN)))
#			}
#			if (!("Rsomaticres" %in% ls())){
#				Rsomaticres = Rsomatic
#			}else{
#				Rsomaticres=rbind(Rsomaticres,Rsomatic)
#			}	 
#		}
#	}
#	out=list(germline=Rgermlineres,segment=segment)
#	if ("Rsomaticres" %in% ls()){
#		out$somatic=Rsomaticres
#	}
#	return(out)
#}



RNAalpha<-function(segment,RNA,DNA,DNApurity,Rpurity){
	segout=foreach (i = 1:dim(segment)[1], .combine=rbind) %dopar% segalpha(segment,RNA,DNA,DNApurity,i,Rpurity)
	alpha=RNApurity(segout[,1],segout[,2])
	return(alpha)
}


segmentASE<-function(segment,alphaD,alphaR,DNA,RNA,ratioS){
	for (i in 1:dim(segment)[1]){
		subRNA=RNA[RNA$chr==as.character(segment$chr[i])&RNA$pos>=segment$startpos[i]&RNA$pos<=segment$endpos[i],]
		index=match(paste(subRNA$chr,subRNA$pos),paste(DNA$chr,DNA$pos))
		subRNA=subRNA[!is.na(index),]
		subDNA=DNA[index[!is.na(index)],]
		segment$NumR[i]=median(subRNA$refNum+subRNA$altNum)
		segment$NumD[i]=median(subDNA$refNumT+subDNA$altNumT)
	}
	subseg=segment[!is.na(segment$NumR)&!is.na(segment$NumD),]
	subseg$ratioR=subseg$NumR/subseg$NumD
	subseg$ratioE=(2*(1-alphaR)+alphaR*(subseg$nMajor+subseg$nMinor))/(2*(1-alphaD)+alphaD*(subseg$nMajor+subseg$nMinor))
	subseg$ratioEW=ratioS*subseg$ratioE
	subseg$ratioRE=subseg$ratioR/subseg$ratioEW
	ratioRE1 <- normalmixEM(subseg$ratioRE, lambda = rep(1, 3)/3, mu = c(min(subseg$ratioRE),1, max(subseg$ratioRE)), sigma = sd(subseg$ratioRE))
	res=matrix(nr=3,nc=3)
	row.names(res)=c("mu","pai","sd")
	res[1,]=ratioRE1$mu
	res[2,]=ratioRE1$lambda
	res[3,]=ratioRE1$sigma
	k=which.max(res[2,])
	f<-res[2,1]*dnorm(subseg$ratioRE,res[1,1],res[3,1])+res[2,2]*dnorm(subseg$ratioRE,res[1,2],res[3,2])+res[2,3]*dnorm(subseg$ratioRE,res[1,3],res[3,3])
	subseg$BayesP=1-res[2,k]*dnorm(subseg$ratioRE,res[1,k],res[3,k])/f
	return(subseg)
}
SNPASE<-function(segment,DNA,RNA,DNApurity,Rpurity,ratioS,somatic){
	for (i in 1:dim(segment)[1]){
		nMajor=segment$nMajor[i]
		nMinor=segment$nMinor[i]
		TCN=nMajor+nMinor
		subdata=RNA[RNA$chr==as.character(segment$chr[i])&RNA$pos>=segment$startpos[i]&RNA$pos<=segment$endpos[i],]
		subdata=subdata[subdata$altNum+subdata$refNum>=10,]
		germline=subdata[subdata$type=="germline"&!is.na(subdata$tfrac),]
		subDNA=DNA[DNA$chr==as.character(segment$chr[i])&DNA$pos>=segment$startpos[i]&DNA$pos<=segment$endpos[i],]
		index=match(germline$pos,subDNA$pos)
		germline=germline[!is.na(index),]
		germline$ND=subDNA$refNumT[index[!is.na(index)]]+subDNA$altNumT[index[!is.na(index)]]
		Rsomatic=subdata[subdata$type=="somatic"&!is.na(subdata$tfrac),]
		index=match(paste(Rsomatic$chr,Rsomatic$pos),paste(somatic$chr,somatic$pos))
		Rsomatic=Rsomatic[!is.na(index),]
		Rsomatic$DMCN=somatic$SMCN[index[!is.na(index)]]
		Rsomatic$Dfrac=somatic$tfrac[index[!is.na(index)]]
		if (dim(germline)[1]>=1){
			germline$tfrac=germline$altNum/(germline$refNum+germline$altNum)
			germline$Dfrac=subDNA$altNumT[match(germline$pos,subDNA$pos)]/(subDNA$altNumT[match(germline$pos,subDNA$pos)]+subDNA$refNumT[match(germline$pos,subDNA$pos)])
			germline$N=germline$refNum+germline$altNum
			germline$Dmajor=nMajor
			germline$Dminor=nMinor
			germline=MajorAlle(nMajor=nMajor,nMinor=nMinor,germline=germline)
			TCNRNA=(germline$N/germline$ND/ratioS*(2*(1-DNApurity)+DNApurity*TCN)-2*(1-Rpurity))/Rpurity
			germline$ratio=log(germline$N/germline$ND,base=2)
			if (dim(germline)[1]>1){
				peak=density(germline$ratio)$x[which.max(density(germline$ratio)$y)]
			}else{
				peak=germline$ratio
			}
			RTCN1=median(TCNRNA)
			RRTCN=germline$ratio-peak+RTCN1
			delt=0
			if (sum(RRTCN<0)>0){
				delt=min(RRTCN)
				RRTCN=RRTCN-min(RRTCN)
			}
			germline$RTCN=RRTCN
			DRTCN=germline$RTCN*TCN/round(median(RRTCN))
			germline$Rmajor=(germline$y/germline$N*(2*(1-Rpurity)+Rpurity*germline$RTCN)-(1-Rpurity))/Rpurity
			germline$Rminor=((1-germline$y/germline$N)*(2*(1-Rpurity)+Rpurity*germline$RTCN)-(1-Rpurity))/Rpurity
			if (segment$BayesP[i] <= 0.5){
				RTCN=round(median(DRTCN))
				germline$RTCN=DRTCN
				germline$Rmajor=(germline$y/germline$N*(2*(1-Rpurity)+Rpurity*germline$RTCN)-(1-Rpurity))/Rpurity
				germline$Rminor=((1-germline$y/germline$N)*(2*(1-Rpurity)+Rpurity*germline$RTCN)-(1-Rpurity))/Rpurity
				segment$RTCN[i]=RTCN
				segment$Rmajor[i]=round(median(germline$Rmajor))
				segment$Rminor[i]=round(median(germline$Rminor))
				#germline$RMajor=round(median(germline$Rmajor))
				#germline$RMinor=round(median(germline$Rminor))
			}else{
				RTCN=median(germline$RTCN)
				segment$RTCN[i]=RTCN
				segment$Rmajor[i]=round(median(germline$Rmajor))
				segment$Rminor[i]=round(median(germline$Rminor))
				#germline$RMajor=round(median(germline$Rmajor))
				#germline$RMinor=round(median(germline$Rminor))
			}
			if (dim(germline)[1]>=20){
				if (segment$BayesP[i] <= 0.5){
					germline=SNPclass(germline=germline,TCN=RTCN,alpha=Rpurity,nMajor=nMajor)
				}else{
					germline=SNPclass2(germline,DTCN=TCN,RTCN=RTCN,alpha=Rpurity,nMajor=nMajor)
				}
			}else{
				germline$Rmajor=(germline$y/germline$N*(2*(1-Rpurity)+Rpurity*TCN)-(1-Rpurity))/Rpurity
				germline$Rminor=((1-germline$y/germline$N)*(2*(1-Rpurity)+Rpurity*TCN)-(1-Rpurity))/Rpurity
				germline$Rminor[germline$tfrac==0]=0
				k=which(germline$tfrac!=0&germline$Rminor<0)
				germline$Rminor[k]=abs(germline$Rminor[k])
				germline$RTCN=germline$Rmajor+germline$Rminor
				germline$BayesP=1-exp(-abs(germline$RTCN-TCN))
			}
			germline$Rminor[germline$tfrac==0]=0
			k=which(germline$tfrac!=0&germline$Rminor<0)
			germline$Rminor[k]=abs(germline$Rminor[k])
			germline$RTCN=germline$Rmajor+germline$Rminor
			if (!("Rgermlineres" %in% ls())){
				Rgermlineres = germline
			}else{
				Rgermlineres=rbind(Rgermlineres,germline)
			}
		}
		if (dim(Rsomatic)[1]>=1){
			index=match(paste(Rsomatic$chr,Rsomatic$pos,sep=":"),paste(somatic$chr,somatic$pos,sep=":"))
			Rsomatic$ND=somatic$Tsum[index]
			Rsomatic$ratio=(Rsomatic$refNum+Rsomatic$altNum)/Rsomatic$ND
			if (is.na(RTCN)){
				RTCN=median(germline$RTCN)
			}
			Rsomatic$RTCN=Rsomatic$ratio-median(log(germline$N/germline$ND,base=2))+RTCN
			Rsomatic$RMCN=Rsomatic$tfrac*(2*(1-Rpurity)+Rpurity*Rsomatic$RTCN)/Rpurity
			Rsomatic$Ddelt=Rsomatic$DMCN-round(Rsomatic$DMCN)
			Rsomatic$delt=Rsomatic$RMCN-Rsomatic$DMCN
			if (dim(Rsomatic)[1]>=10){
				if (RTCN==TCN){
					Rsomatic=SNVclass(Rsomatic)
				}else{
					Rsomatic=SNVclass2(Rsomatic)
				}
			}else{
				Rsomatic$BayesP=1-exp(-abs(Rsomatic$RMCN-Rsomatic$DMCN))
			}
			if (!("Rsomaticres" %in% ls())){
				Rsomaticres = Rsomatic
			}else{
				Rsomaticres=rbind(Rsomaticres,Rsomatic)
			}
		}
	}
	out=list(germline=Rgermlineres,segment=segment)
	if ("Rsomaticres" %in% ls()){
		out$somatic=Rsomaticres
	}
	return(out)
}

RNArun<-function(SNPinput,RNAinput,DNAout){
	print.noquote("estimation at RNA-level")
	DNA=read.csv(SNPinput,,header=TRUE,sep="\t")
	colnames(DNA)=c("chr","pos","ref","alt","refNumN","altNumN","refNumT","altNumT")
	DNA=DNA[(DNA$refNumN+DNA$altNumN)>=10&(DNA$refNumT+DNA$altNumT)>=10,]
	rlength=nchar(as.character(DNA$ref))
	alength=nchar(as.character(DNA$alt))
	DNA=DNA[rlength==1&alength==1,]
	DNA=DNA[!is.na(DNA$pos),]
	RNA=read.csv(RNAinput,header=TRUE,sep="\t")
	colnames(RNA)=c("chr","pos","ref","alt","refNum","altNum","type")
	RNA=RNA[(RNA$refNum+RNA$altNum)>=10,]
	rlength=nchar(as.character(RNA$ref))
	alength=nchar(as.character(RNA$alt))
	RNA=RNA[rlength==1&alength==1,]
	if (dim(RNA)[1]>10){
		RNA$tfrac=RNA$altNum/(RNA$altNum+RNA$refNum)
		segment=DNAout$segment
		somatic=DNAout$somatic
		DNApurity=DNAout$alpha
		ratioS=median(RNA$altNum+RNA$refNum)/median(DNA$altNumT+DNA$refNumT)
		registerDoMC(cores = 4)
		print.noquote("Inferring tumor purity of RNA data")
		Rpurity=RNAalpha(segment=segment,RNA=RNA,DNA=DNA,DNApurity=DNApurity,Rpurity="NA")
		print.noquote("Inferring allelic imbalance of segments")
		segment=segmentASE(segment,alphaD=DNApurity,alphaR=Rpurity,DNA=DNA,RNA=RNA,ratioS=ratioS)
		print.noquote("Inferring allelic imbalance of SNPs/SNVs")
		RNAout=SNPASE(segment,DNA,RNA,DNApurity,Rpurity,ratioS,somatic)
		RNAout$alpha=Rpurity
		return(RNAout)
	}else{
		print.noquote("RNA data is insufficient!")
	}
	
}