Pathogen Variant Calling Pipeline
=================================

This pipeline aligns reads to a reference and calls variants using BWA and samtools.

Written by Shaun Jackman <sjackman@gmail.com>

Required Software
=================

 * BWA
 * samtools
 * tabix
 * FastTree

Input
=====

* FASTQ files of short read sequencing data
* A reference genome sequence

Output
======

 * a BAM file of reads aligned to the reference
 * a VCF file of variant calls
 * a FASTA file of the genotype of those variants
 * a dendrogram derived from those genotypes
