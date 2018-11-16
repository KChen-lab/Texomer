'''
Created on May 23, 2017

@author: FWang9
'''
import sys,os,time
import subprocess
from optparse import OptionParser
from tempfile import gettempdir
def VarscanTrans(varscan,case):
    data=open(varscan)
    germlineout=open(case+".germline","w")
    somaticout=open(case+".somatic","w")
    out="chr"+"\t"+"position"+"\t"+"ref"+"\t"+"alt"+"\t"+"refNumN"+'\t'+"altNumN"+"\t"+"refNumT"+"\t"+"altNumT"
    print >> germlineout,out
    print >> somaticout,out
    next(data)
    for line in data:
        line=line[0:-1].split("\t")
        if "chr" not in line[0]:
            chr="chr"+line[0]
        else:
            chr=line[0]
        position=line[1]
        ref=line[2]
        var=line[3]
        refNumN=line[4]
        altNumN=line[5]
        refNumT=line[8]
        altNumT=line[9]
        gt=line[7]
        type=line[12]
        pgerm=float(line[13])
        psoma=float(line[14])
        if type == "Germline" and gt not in ["A","C","G","T"] and pgerm < 0.01:
            out=chr+"\t"+position+"\t"+ref+"\t"+var+"\t"+refNumN+"\t"+altNumN+"\t"+refNumT+"\t"+altNumT
            print >> germlineout,out
        if type == "LOH" and psoma < 0.01:
            out=chr+"\t"+position+"\t"+ref+"\t"+var+"\t"+refNumN+"\t"+altNumN+"\t"+refNumT+"\t"+altNumT
            print >> germlineout,out
        if type == "Somatic" and psoma < 0.01:
            out=chr+"\t"+position+"\t"+ref+"\t"+var+"\t"+refNumN+"\t"+altNumN+"\t"+refNumT+"\t"+altNumT
            print >> somaticout,out
def SNVexp(germline,somatic,samfile,outputname,case,bedtools):
    chro=range(23)
    chro.remove(0)
    chrom=[]
    for ele in chro:
        chrom.append("chr"+str(ele))
    chrom.append("chrX")
    chrom.append("chrY")
    SNV={}
    chrindex=0
    SNPdata=open(germline)
    next(SNPdata)
    for line in SNPdata:
        line=line[0:-1].split("\t")
        if "chr" not in line[0]:
            chr="chr"+line[0]
        else:
            chr=line[0]
        pos=int(line[1])
        ref=line[2]
        alt=line[3]
        SNV.setdefault(chr,{})[pos]="germline:"+ref+":"+alt
    SNPdata.close()
    somaticdata=open(somatic)
    next(somaticdata)
    for line in somaticdata:
        line=line[0:-1].split("\t")
        if "chr" not in line[0]:
            chr="chr"+line[0]
        else:
            chr=line[0]
        pos=int(line[1])
        ref=line[2]
        alt=line[3]
        SNV.setdefault(chr,{})[pos]="somatic:"+ref+":"+alt
    somaticdata.close()
    sam=open(samfile)
    tempout=open(case+".bed","w")
    for line in sam:
        line=line[0:-1].split("\t")
        readid=line[0]
        if "@" not in readid:
            if "chr" not in line[2]:
                chr="chr"+line[2]
            else:
                chr=line[2]
            start=line[3]
            seq=line[9]
            out=chr+"\t"+start+"\t"+str(int(start)+len(seq)-1)+"\t"+readid+"\t"+seq
            print >> tempout,out
    sam.close()
    tempout.close()
    output=open(outputname,'w')
    out="chr\tpos\tref\talt\trefNum\taltNum\ttype"
    print >> output,out
    snvpos=open(case+".snptemp.bed","w")
    chrom1=list(set(chrom).intersection(SNV.keys()))
    for chromosome in chrom1:
        temp=SNV[chromosome]
        for key in sorted(temp.keys()):
            type=temp[key].split(":")[0]
            ref=temp[key].split(":")[1]
            alt=temp[key].split(":")[2]
            out=chromosome+"\t"+str(key)+"\t"+str(key)+"\t"+ref+"\t"+alt+"\t"+type
            print >> snvpos,out
    snvpos.close()
    os.system(bedtools+" intersect -a "+case+".bed -b "+case+".snptemp.bed -wa -wb > "+case+".snp.read.bed")
    snpdata=open(case+".snp.read.bed")
    read={}
    SNPpos={}
    for line in snpdata:
        line=line[0:-1].split("\t")
        chr=line[0]
        start=int(line[1])
        end=int(line[2])
        readid=line[3]
        sequen=line[4]
        pos=int(line[6])
        ref=line[8]
        alt=line[9]
        if sequen[pos-start] == ref:
            SNPpos.setdefault(chr,[]).append(pos)
            read.setdefault(chr+":"+str(pos),[]).append("REF\t"+readid)
        elif sequen[pos-start] == alt:
            SNPpos.setdefault(chr,[]).append(pos)
            read.setdefault(chr+":"+str(pos),[]).append("ALT\t"+readid)
    snpdata.close()
    for chromosome in chrom:
        temppos=sorted(list(set(SNPpos[str(chromosome)])))
        temp=SNV[str(chromosome)]
        for key in temppos:
            type=temp[key].split(":")[0]
            ref=temp[key].split(":")[1]
            alt=temp[key].split(":")[2]
            nref=0
            nalt=0
            for id in list(set(read[str(chromosome)+":"+str(key)])):
                if "REF\t" in id:
                    nref=nref+1
                elif "ALT\t" in id:
                    nalt=nalt+1
            out=str(chromosome)+"\t"+str(key)+"\t"+ref+"\t"+alt+"\t"+str(nref)+"\t"+str(nalt)+"\t"+type
            print >> output,out
    output.close()
def datatrans(input,output):
    data=open(input)
    outputdata=open(output,"w")
    line=next(data)[0:-1]
    print >> outputdata,line
    for line in data:
        print >> outputdata,"chr"+line[0:-1]
    data.close()
    outputdata.close()
def SNPcalling(Tumor,Normal,samtools,picard,ref,Varscan):
    os.system(samtools+" index "+Tumor)
    os.system(samtools+" index "+Normal)
    tumorpile="tumor.pileup"
    normalpile="normal.pileup"
    commond=samtools+" mpileup -f "+ref+" "+Tumor+" > "+tumorpile
    os.system(commond)
    commond=samtools+" mpileup -f "+ref+" "+Normal+" > "+normalpile
    os.system(commond)
    commond="java -jar "+Varscan+" somatic "+normalpile+" "+tumorpile+" "+"case.varscan"
    os.system(commond)
    VarscanTrans(varscan="case.varscan.snp",case="case")

def getPath(path):
    path1=path.split("/")
    if path1[0] == ".":
        if (len(path1)==1):
            newpath=os.getcwd()
        else:
            newpath=os.getcwd()
            for ele in path1:
                if ele !=".":
                    newpath=newpath+"/"+ele
    elif path1[0]=="..":
        i = 0
        for ele in path1:
            if ele == "..":
                i=i+1
        path2=os.getcwd()
        path2=path2.split("/")
        newpath="/"+path2[0]
        if len(path2)-i > 1:
            for j in range(1,len(path2)-i):
                newpath=newpath+"/"+path2[j]
    else:
        newpath=path
    return newpath


def main():
    usage = "usage: python %prog [-t <tumor bam file>] [-n <normal bam file>] [-r <RNA bam file>] [-v <Varscan mutation calling output>] [-g <Defiend germline mutation input file>] [-s <Defined somatic mutation input file>] [-o <output path>] [-u <Optimization>] [-e <Defined expression file of mutation>] [-f location of reference sequence] -p <Texomer path> -I <input form> "
    description = "[options] is optional. The path of Texomer and the form of input file are required in Texomer. Three forms of input file are allowed: BAM, Varscan and Defined. If -I BAM is selected, -t and -n are required. Texomer requires bam file is alligned based on GRCH38. If -I Varscan is selected, -v is required. If -I Defined is selected, -g and -s are required. RNA-seq bam file is optional for Texomer. Texomer would do estimation only at DNA level if without RNA data input."
    op = OptionParser(version="%prog 2.0",description=description,usage=usage,add_help_option=False)
    op.add_option("-h","--help",action="help",
                  help="Show this help message and exit.")
    op.add_option("-p","--Texomer",dest="Texomer",type="str",
                  help="the path of Texomer")
    op.add_option("-I","--Input",dest="Input",type="str",
                  help="the form of input file: BAM, Varscan, and Defined")
    op.add_option("-t","--Tumor",dest="Tumor",type="str",
                  help="Tumor WES bam file.")
    op.add_option("-n","--Normal",dest="Normal",type="str",
                  help="Normal WES bam file.")
    op.add_option("-o","--outpath",dest="outpath",type="str",
                  help="the output path.")
    op.add_option("-r","--RNA",dest="RNA",type="str",
                  help="RNA-seq bam file")
    op.add_option("-g","--germline",dest="germline",type="str",
                  help="You can input your own germline mutation file together with somatic mutation file (-s). The file include 8 columns: chromosome, position, RefAllele, Altallele, read counts of RefAllele in normal, read counts of Altallele in normal, read counts of RefAllele in tumor and read counts of Altallele in tumor with header seperated by tab. \t")
    op.add_option("-s","--somatic",dest="somatic",type="str",
                  help="You can input your own somatic mutation file together with germline mutation file (-g). The file include 8 columns: chromosome, position, RefAllele, Altallele, read counts of RefAllele in normal, read counts of Altallele in normal, read counts of RefAllele in tumor and read counts of Altallele in tumor with header seperated by tab. \t")
    op.add_option("-v","--varscan",dest="varscan",type="str",
                  help="the output of Varscan2 based on somatic calling")
    op.add_option("-e","--snvexpress",dest="snvexpress",type="str",
                  help="The allelic read count of mutation from RNA-seq including 7 columns: chromosome, position, ref, alt, refnum, altnum and type (germline or somatic)")
    op.add_option("-u","--iter",dest="iter",type="int",
                  help="Optimization using somatic mutation, 0 corresponding to no optimization and 1 is optimization. The default = 1")
    op.add_option("-f","--reference",dest="reference",type="str",
                  help="The location of reference sequence.")
    (options,args) = op.parse_args()
    if not options.Texomer or not options.Input:
        op.print_help()
        sys.exit(1)
    if options.Input not in ["BAM","Varscan","Defined"]:
        print "Please select input form from BAM, varscan and Defined!"
        op.print_help()
        sys.exit(1)
    currentPath=os.getcwd()
    Texomerpath=options.Texomer
    Texomerpath=getPath(Texomerpath)
    samtools=Texomerpath+"/samtools"
    bedtools=Texomerpath+"/bedtools"
    picard=Texomerpath+"/picard.jar"
    Varscan=Texomerpath+"/VarScan.v2.3.9.jar"
    if not options.outpath:
        outpath=os.getcwd()
    else:
        outpath=options.outpath
        outpath=getPath(outpath)
        outpath1=outpath.split('/')
        path1="/"
        for i in range(1,(len(outpath1)-1)):
            path1=path1+outpath1[i]+"/"
        if outpath1[len(outpath1)-1] not in os.listdir(path1):
            os.chdir(path1)
            os.system("mkdir "+outpath1[len(outpath1)-1])
    os.chdir(outpath)
    os.system("mkdir temp")
    os.chdir(currentPath)
    #Rscript = options.Rscript
    if options.Input == "Varscan":
        if not options.varscan:
            print "Please input mutation file matched with Varscan output through -v."
            sys.exit(1)
        else:
            varscan=options.varscan
            print "Transfer varscan2 output"
            varscan=getPath(varscan)
            case=varscan.split("/")[len(varscan.split("/"))-1]
            os.chdir(outpath+"/temp")
            VarscanTrans(varscan,case)
            germlinename = os.getcwd()+"/"+case+".germline"
            somaticname = os.getcwd()+"/"+case+".somatic"
    elif options.Input == "BAM":
        if options.Tumor and options.Normal:
            if not options.reference:
                print "Please input the location of refereence sequence through -f."
                sys.exit(1)
            else:
                ref=options.reference
                ref=getPath(ref)
                Tumor=options.Tumor
                Tumor=getPath(Tumor)
                Normal=options.Normal
                Normal=getPath(Normal)
                case="case"
                os.chdir(outpath+"/temp")
                SNPcalling(Tumor=Tumor,Normal=Normal,samtools=samtools,picard=picard,ref=ref,Varscan=Varscan)
                germlinename = os.getcwd()+"/"+case+".germline"
                somaticname = os.getcwd()+"/"+case+".somatic"
                os.system("rm "+Tumor+".bai")
                os.system("rm "+Normal+".bai")
        else:
            print "Please input bam file through -t and -n."
            sys.exit(1)
    elif options.Input == "Defined":
        if options.germline and options.somatic:
            germlinename = options.germline
            germlinename=getPath(germlinename)
            somaticname = options.somatic
            somaticname=getPath(somaticname)
        else:
            print "Please input defined germline and somatic mutation files through -g and -s."
            sys.exit(1)
    os.chdir(outpath+"/temp")
    data=open(somaticname)
    next(data)
    line=next(data)
    data.close()
    if "chr" not in line:
        case=somaticname.split("/")[len(somaticname.split("/"))-1]
        somaticoutname = os.getcwd()+"/"+case+".somatic"
        datatrans(input=somaticname,output=somaticoutname)
        somaticname=somaticoutname
    data=open(germlinename)
    next(data)
    line=next(data)
    data.close()
    if "chr" not in line:
        case=germlinename.split("/")[len(germlinename.split("/"))-1]
        germlineoutname = os.getcwd()+"/"+case+".germline"
        datatrans(input=germlinename,output=germlineoutname)
        germlinename=germlineoutname
    os.chdir(currentPath)
    if not options.RNA:
        if not options.snvexpress:
            print "No RNA data! Only estimate at DNA level!"
            RNAname="NA"
        else:
            RNAname=options.snvexpress
            RNAname=getPath(RNAname)
            data=open(RNAname)
            next(data)
            line=next(data)[0:-1].split("\t")
            data.close()
            if "chr" not in line[0]:
                case=outputname.split("/")[len(outputname.split("/"))-1]
                os.chdir(outpath+"/temp")
                datatrans(input=outputname,output=case)
                #outputname=case
    else:
        RNAname = options.RNA
        RNAname=getPath(RNAname)
        case=RNAname.split("/")[len(RNAname.split("/"))-1]
        print "Transfer bam to sam"
        os.chdir(outpath+"/temp")
        os.system(samtools+ " index "+RNAname)
        os.system(samtools+" view "+RNAname+" > "+RNAname+".sam")
        samfilename=RNAname+".sam"
        outputname=case+".SNV"
        if os.path.exists(outputname):
            os.remove(outputname)
        print "Align mutations to RNA data"
        SNVexp(germline=germlinename,somatic=somaticname,samfile=samfilename,outputname=outputname,case=case,bedtools=bedtools)
        os.system("rm "+RNAname+".sam")
        os.system("rm "+RNAname+".bai")
        RNAname=outputname
    print "Infering at DNA level"
    if not options.iter:
        optindex=1
    else:
        optindex=options.iter
    os.system("Rscript "+Texomerpath+"/Texomer.R "+germlinename+" "+somaticname+" "+RNAname+" "+Texomerpath+" "+outpath+" "+str(optindex))
    os.chdir(outpath)
    os.system("rm -r temp")




if __name__ == "__main__":

    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me, see you!\n")
        sys.exit(0)
