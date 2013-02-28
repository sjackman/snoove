Pathogen Variant Calling Pipeline
=================================

This pipeline aligns reads to a reference and calls variants using BWA and samtools.

Written by Shaun Jackman <sjackman@gmail.com>

Required Software
=================

 * BWA 0.6.2
 * samtools 0.1.18
 * tabix 0.2.6
 * PhyML 20120412

Input
=====

* FASTQ files of short read sequencing data
* A reference genome sequence

Output
======

 * a BAM file of reads aligned to the reference
 * a VCF file of variant calls
 * a FASTA file of the genotype of those variants
 * a SNP tree based on those genotypes
