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
        if type == "Germline" and gt not in ["A","C","G","T"]:
            out=chr+"\t"+position+"\t"+ref+"\t"+var+"\t"+refNumN+"\t"+altNumN+"\t"+refNumT+"\t"+altNumT
            print >> germlineout,out
        if type == "Somatic":
            out=chr+"\t"+position+"\t"+ref+"\t"+var+"\t"+refNumN+"\t"+altNumN+"\t"+refNumT+"\t"+altNumT
            print >> somaticout,out
def SNVexp(germline,somatic,samfile,outputname,case):
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
    os.system("intersectBed -a "+case+".bed -b "+case+".snptemp.bed -wa -wb > "+case+".snp.read.bed")
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
    chro=list(set(chrom).intersection(SNPpos.keys()))
    for chromosome in chro:
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

def main():
    usage = "usage: python %prog <-g germline> <-s somatic> <-r RNA.bam> <-p script path> <-v varscan output>[...]"
    description = "Input files require RNA-seq data, germline and somatic mutations. We also can convert the output of the somatic programs of Varscan2 through -v if you don't have germline and somatic mutation files. For example: python %prog -g germline -s somatic -R RNA.bam -i RNA.bam.bai -o output. -i and -o are not necessary."
    op = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    op.add_option("-h","--help",action="help",
                  help="Show this help message and exit.")
    op.add_option("-p","--Rscript",dest="Rscript",type="str",
                  help="the path of Texomer")
    op.add_option("-g","--germline",dest="germline",type="str",
                  help="You can input your own germline mutation file if no -v. The file name of germline input file include 8 columns: chromosome, position, RefAllele, Altallele, read counts of RefAllele in normal, read counts of Altallele in normal, read counts of RefAllele in tumor and read counts of Altallele in tumor with header seperated by tab. \t")
    op.add_option("-s","--somatic",dest="somatic",type="str",
                  help="You can input your own somatic mutation file if no -v. The file name of somatic input file include 8 columns: chromosome, position, RefAllele, Altallele, read counts of RefAllele in normal, read counts of Altallele in normal, read counts of RefAllele in tumor and read counts of Altallele in tumor with header seperated by tab. \t")
    op.add_option("-v","--varscan",dest="varscan",type="str",
                  help="the output of Varscan2 based on somatic calling")
    op.add_option("-r","--RNA",dest="RNA",type="str",
                  help="bam file of RNA-seq")
    op.add_option("-i","--bai",dest="bai",type="str",
                  help="index bam file of RNA-seq. this is optional.")
    op.add_option("-e","--snvexpress",dest="snvexpress",type="str",
                  help="The allelic read count of mutation from RNA-seq including 7 columns: chromosome, position, ref, alt, refnum, altnum and type (germline or somatic)")
    op.add_option("-t","--iter",dest="iter",type="int",
                  help="optimal through somatic mutation, 0 corresponding to no optimation and 1 is optimation. The default = 1")
    op.add_option("-o","--outpath",dest="outpath",type="str",
                  help="the output path. This is optional")
    (options,args) = op.parse_args()
    if not options.Rscript:
        op.print_help()
        sys.exit(1)
    if not options.outpath:
        outpath=os.getcwd()
    else:
        outpath=options.outpath
    os.chdir(outpath)
    os.system("mkdir temp")
    os.chdir("temp")
    Rscript = options.Rscript
    if not options.germline or not options.somatic:
        if not options.varscan:
            op.print_help()
            sys.exit(1)
        else:
            varscan=options.varscan
            print "Transfer varscan2 output"
            case=varscan.split("/")[len(varscan.split("/"))-1]
            VarscanTrans(varscan,case)
            germlinename = os.getcwd()+"/"+case+".germline"
            if not options.somatic:
                somaticname = os.getcwd()+"/"+case+".somatic"
            else:
                somaticname= options.somatic
    else:
        germlinename = options.germline
        somaticname = options.somatic    
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
    if not options.snvexpress:
        if not options.RNA:
            print "No RNA data! Only estimate at DNA level!"
            outputname="NA"
            #op.print_help()
            #sys.exit(1)
        else:
            RNAname = options.RNA
            case=RNAname.split("/")[len(RNAname.split("/"))-1]
            if RNAname.split(".")[len(RNAname.split("."))-1] == "bam":
                if not options.bai:
                    print "Generate index"
                    os.system("samtools index "+RNAname)
                print "Transfer bam to sam"
                os.system("samtools view "+RNAname+" > "+RNAname+".sam")
                samfilename=RNAname+".sam"
            else:
                samfilename=RNAname
            outputname=case+".SNV"
            if os.path.exists(outputname):
                os.remove(outputname)
            print "Align mutations to RNA data"
            SNVexp(germline=germlinename,somatic=somaticname,samfile=samfilename,outputname=outputname,case=case)
    else:
        outputname=options.snvexpress
        data=open(outputname)
        next(data)
        line=next(data)[0:-1].split("\t")
        data.close()
        if "chr" not in line[0]:
            case=outputname.split("/")[len(outputname.split("/"))-1]
            datatrans(input=outputname,output=case)
            outputname=case
    print "Infering at DNA level"
    if not options.iter:
        optindex=1
    else:
        optindex=options.iter
    os.system("Rscript "+Rscript+"/Texomer.R "+germlinename+" "+somaticname+" "+outputname+" "+Rscript+" "+outpath+" "+str(optindex))
    os.chdir(outpath)
    os.system("rm -r temp")

    


if __name__ == "__main__":
    
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me, see you!\n")
        sys.exit(0)    




