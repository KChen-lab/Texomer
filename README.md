DACRE-scan

Author: Fang Wang
Email: fwang9@mdanderson.org or kchen3@mdanderson.org
Draft date: Feb. 12, 2018

Description
===========
DACRE-scan is a tool used in cancer genomic studies to discover single nucleotide variants (both somatic and germline) 
with differential allelic cis-regulatory effects(DACRE) through integration whole- exome sequencing (WES) data and whole transcriptome sequencing (WTS) data.
It can report estimations of purity at both DNA and RNA levels, heterogeneity at DNA level, allelic specific copy number and expression level corresponding
to tumor component and DACRE score of mutated allele. If the input only include information of DNA level, we may only report the estimates at DNA level and 
focused on somatic mutaions. 

System requirements and dependency
==================================
DACRE-scan runs on a x86_64 Linux system. It depends on samtools and bedtools to extract reads of variants from WTS data, R(version >= 3.2)
for DACRE running which depends on R package bbmle, emdbook, copynumber,TitanCNA, facets, mixtools, ASCAT and Sequenza. I have already put R packeages 
in the release, You can run DACRE-scan directly.

Installation
============
Please download and copy the distribution to your specific location. For example, the downloaded distribuition is DACRE_scan.tar.gz.
	Type 'tar zxvf DACRE_scan.tar.gz'
Then, run SNVexpress.py in DACRE_scan folder.

Usage
=====
Input files: 

DNA input files: 
	two kinds of input files were allowed in DACRE_scan:

	(1) input the output of varscan including germline and somatic mutations through -v;

	(2) two seperate files one is germline mutation through -g and another is somatic mutation through -s. 
	The format of these two files are the same, including 8 columns: chromosome (start from "chr"), 
	position, Reference allele, alternative allele, ref allelic read counts in normal, alt allelic 
	read counts in normal, ref and alt allelic read counts in tumor with header and tab separate.

RNA input files:
	two kinds of input files were allowed in ARDE:

	(1) bam file from RNA-seq. Load samtools and bedtools first if you input bam file.

	(2) allelic read counts of mutation from RNA-seq inlcuding 7 cloumns: chromosome (start from "chr"),
	positive, reference allele, alternative allele, reference and alternative allelic read count from
	RNA-seq, type of mutation (germline or somatic) with header and tab separate.    

Output file: Multiple files would be output. If you input RNA file, ARDE would output both at DNA and RNA level. otherwise, it 
	would output estimation at DNA level. 

	(1) .segment file: position of segmentation; allele-specific copy number (Dmajor and Dminor), 
	allele-specific expression levels(Rmajor and Rminor), Posterior probability of discordant expression corresponding to each segment.

	(2) .mutation file: position of mutation, reference(ref) and alternative(alt) allele, allelic read count from RNA-seq(refNum and altNum),
	mutation type(germline or somatic), allele specific copy number (altD corresponding to alternative allele and wildD corresponding to reference allele),
	allele specific expression level (altR corresponding to alternative allele and wildR corresponding to reference allele),
	posterior probability of discordant expression and DACRE score for each mutation.

	(3) .summary file: purity at DNA and/or RNA level; ploidy and heterogeneity at DNA level

run DACRE-scan:
***The python script for easy (hopefully) run of DACRE_scan is in the release directory. You can tune the
parameters as you wish.***

Python DACRE_scan.py –p Rscript/path –g germline.input –s somatic.input –r RNA.bam <-t iter.optimal.index> <-e mutation.expression> <-o output.path>

About the default parameters
========================
DACRE_scan optimizes estimation of purity and allele specific copy number through combing somatic mutation. So default -t is 1 corresponding to do the optimized iteration.

You can set -t 0 if you don't want to do optimization.

Thank you very much for testing and using DACRE_scan. We appreciate so much for
your feedback!


Example
=====
python SNVexpress.py -g TCGA-3C-AAAU-01A-11D-A41F-09.snp.germline.input -s TCGA-3C-AAAU-10A-01D-A41F-09.mutect.vcf -e TCGA-3C-AAAU.RNA.SNV -p segACN
