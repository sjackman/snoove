#!/usr/bin/make -rRf
# Map reads to a reference, call variants and create a dendrogram
# Written by Shaun Jackman <sjackman@gmail.com>

# The data prefix
p=/data/MiSEQ_data/DATA_2013

all: Campylobacter_jejuni_NCTC11168.tree \
	Clostridium_difficile_630.tree \
	MRSA_USA300.tree \
	Salmonella_Enteritidis.tree

.PHONY: all
.DELETE_ON_ERROR:
.SECONDARY:

SampleSheet.csv: $p/*/SampleSheet.csv
	cat $^ |tr -d '\r' >$@

# Actions per species

%.fa.bwt: %.fa
	bwa index $<

%.fa.fai: %.fa
	samtools faidx $<

%.samples: SampleSheet.csv
	awk -F, '$$11=="$*" {print $$1}' $< |sort -u >$@

%-stamp: %.samples ref/%.fa.bwt
	make r=ref/$*.fa `sed 's/$$/\/bwa.bam.bai/' $<`
	touch $@

%.bcf: %.samples %-stamp
	samtools mpileup -Igf ref/$*.fa `sed 's/$$/\/bwa.bam/' $<` >$@

%.vcf: %.bcf
	bcftools view -vcg $< >$*.vcf

%.vcf.gz: %.vcf
	bgzip $<

%.vcf.gz.tbi: %.vcf.gz
	tabix -pvcf $<

%.fa: %.vcf.gz
	zcat $< |bin/vcftofa >$@

%.tree: %.fa
	FastTree -nt $< >$@

# Actions per sample

fastq-files: $p/*/Data/Intensities/BaseCalls/*.fastq.gz
	ls $^ >$@

%/fastq-files: fastq-files
	mkdir -p $*
	grep -m2 /$*_ $< >$@

%/bwa.bam: %/fastq-files
	bwa aln $r `grep _R1_ $<` >$*/R1.sai
	bwa aln $r `grep _R2_ $<` >$*/R2.sai
	bwa sampe $r $*/R1.sai $*/R2.sai `<$<` |samtools view -Su - |samtools sort - $*/bwa
	rm -f $*/R1.sai $*/R2.sai

%.bam.bai: %.bam
	samtools index $<

%/fastqc/stamp: %/fastq-files
	mkdir -p $*/fastqc
	fastqc -o $*/fastqc `<$<`
	touch $@
