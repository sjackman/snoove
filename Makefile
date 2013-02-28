# Map reads to a reference, call variants and create a dendrogram
# Written by Shaun Jackman <sjackman@gmail.com>

# The data prefix
p=/data/MiSEQ_data/DATA_2013

# The reference sequence
r=ref/NC_018521.fa

all: Campylobacter_jejuni_NCTC11168-stamp

.PHONY: all
.DELETE_ON_ERROR:
.SECONDARY:

lanes:
	find $p -name *.fastq.gz |cut -d/ -f5 |sort -u >$@

%.samples:
	grep -h $* $p/*/SampleSheet.csv |cut -d, -f1 |sort -u >$@

%/files:
	mkdir -p $*
	ls $p/*/Data/Intensities/BaseCalls/*.fastq.gz |grep -h $* >$@

%-stamp: %.samples
	make `sed 's/$$/\/bwa.vcf.gz.tbi/' $^`
	touch $@

%.fa.bwt: %.fa
	bwa index $<

%.fa.fai: %.fa
	samtools faidx $<

%/bwa.bam: %/files
	bwa aln $r `grep _R1_ $<` >$*/R1.sai
	bwa aln $r `grep _R2_ $<` >$*/R2.sai
	bwa sampe $r $*/R1.sai $*/R2.sai `<$<` |samtools view -Su - |samtools sort - $*/bwa
	rm -f $*/R1.sai $*/R2.sai

%.bam.bai: %.bam
	samtools index $<

%.bcf: %.bam
	samtools mpileup -uf $r $< >$@

%.vcf: %.bcf
	bcftools view -vcg $< >$*.vcf

%.vcf.gz: %.vcf
	bgzip $<

%.vcf.gz.tbi: %.vcf.gz
	tabix -pvcf $<

%.fa: %.vcf.gz
	zcat $< |bin/vcftofa >$@

%.tree: %.fa
	FastTree $< >$@
