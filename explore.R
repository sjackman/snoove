# Explore a VCF file

#source('http://bioconductor.org/biocLite.R')
#biocLite('VariantAnnotation')

library('lattice')
library('reshape2')
library('VariantAnnotation')

vcf <- readVcf('Salmonella_Enteritidis.vcf.gz', 'NC_018521')

# Histogram of the quality values.
qual <- fixed(vcf)$QUAL
histogram(qual, main='Histogram of quality')

# Violin plot of the quality vs. genotype.
gt <- data.frame(geno(vcf)$GT)
qual.gt <- melt(cbind(qual, gt), id.vars='qual',
	variable.name='sample', value.name='gt')
qual.gt$gt <- factor(qual.gt$gt)
bwplot(qual ~ gt, qual.gt, panel = panel.violin,
	main='Quality vs. genotype',
	ylab='Quality', xlab='Genotype')
bwplot(qual ~ gt | sample, qual.gt, panel = panel.violin,
	main='Quality vs. genotype',
	ylab='Quality', xlab='Genotype')

# Examine genotype qualities.
gt.gq <- data.frame(gt = geno(vcf)$GT, gq=geno(vcf)$GQ)
gt <- melt(t(geno(vcf)$GT), value.name='gt')
gq <- melt(t(geno(vcf)$GQ), value.name='gq')
gt.gq <- cbind(gt[1], gt$gt, gq$gq)
colnames(gt.gq) <- c('sample', 'gt', 'gq')
histogram(gt.gq$gq, main='Histogram of genotype quality')
bwplot(gq ~ gt, gt.gq, panel = panel.violin,
	main='Genotype quality vs. genotype',
	ylab='Genotype quality', xlab='Genotype')
bwplot(gq ~ gt | sample, gt.gq, panel = panel.violin,
	main='Genotype quality vs. genotype',
	ylab='Genotype quality', xlab='Genotype')

# Examine a FASTA file for 'N' characters

library(Biostrings)

fa <- readDNAStringSet('Salmonella_all.fa')
fa <- fa[!names(fa) %in% c('Reference', 'Consensus'),]
n <- max(width(fa))
h <- alphabetFrequency(fa)
pn <- h[,'N'] / n
hist(pn, 20)
table(pn < 0.5)
names(fa)[pn > 0.5]
