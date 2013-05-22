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

# Identify failed samples
len <- max(width(fa))
h <- alphabetFrequency(fa)
pn <- h[,'N'] / len
hist(pn, 20)
boxplot(pn[pn<0.5])
table(pn < 0.5)

# Remove failed samples
names(fa)[pn >= 0.5]
fa <- fa[pn < 0.5]

# Identify loci with many N
c <- t(consensusMatrix(fa))
n <- length(fa)
cn <- c[,'N'] / n
hist(cn, 20)
boxplot(cn[cn < 0.2])
table(cn < 0.2)
table(cn < 0.2) / len

fa.good <- DNAStringSet(apply(
	as.matrix(fa)[,cn < 0.2], 1, function(x) paste(x, collapse='')))
writeXStringSet(DNAStringSet(fa.good), 'Salmonella_good.fa',
	width=max(width(fa.good)))
